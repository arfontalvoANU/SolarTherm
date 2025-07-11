within examples;
model FieldAnnualOpticalInterpolation "Field annual interpolation"
	extends Modelica.Icons.Example;
	// Imports
	import nSI = Modelica.SIunits.Conversions.NonSIunits;
	import SI = Modelica.SIunits;
	import FI = SolarTherm.Models.Analysis.Finances;
	import SolarTherm.Utilities.*;
	import metadata = SolarTherm.Utilities.Metadata_Optics;
	import SolarTherm.Models.Sources.SolarFunctions.*;
	import Modelica.SIunits.Conversions.*;

	// Parameters
	parameter String wea_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Weather/seville_spain_tmy_2005_2014_shifted.motab");
	parameter String pri_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Prices/aemo_vic_2014.motab");
	parameter Real pi = Modelica.Constants.pi;
	parameter Real sigma = Modelica.Constants.sigma;

	parameter nSI.Angle_deg lon = -5.326 "Longitude (+ve East)";
	parameter nSI.Angle_deg lat = 37.558 "Latitude (+ve North)";
	parameter nSI.Time_hour t_zone = 1 "Local time zone (UCT=0)";
	parameter Integer year = 2006 "Meteorological year";

	parameter SolarTherm.Types.Solar_angles angles = SolarTherm.Types.Solar_angles.elo_hra "Angles used in the lookup table file";
	parameter String opt_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Optics/OELTs_Solstice.motab");
	parameter String fileW1 = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Optics/OELTs_Solstice_W1.motab");
	parameter String fileW2 = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Optics/OELTs_Solstice_W2.motab");
	parameter String fileW3 = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Optics/OELTs_Solstice_W3.motab");
	parameter String fileW4 = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Optics/OELTs_Solstice_W4.motab");
	parameter String fileW5 = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Optics/OELTs_Solstice_W5.motab");
	parameter String fileW6 = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Optics/OELTs_Solstice_W6.motab");
	parameter String fileW7 = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Optics/OELTs_Solstice_W7.motab");
	parameter String fileW8 = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Optics/OELTs_Solstice_W8.motab");
	parameter String fileW9 = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Optics/OELTs_Solstice_W9.motab");
	parameter Real metadata_list[23] = metadata(opt_file);
	parameter Integer n_heliostat = integer(metadata_list[1]) "Number of heliostats";
	parameter SI.Area A_heliostat = metadata_list[2] "Heliostat module reflective area";
	parameter SI.Area A_field = n_heliostat * A_heliostat "Heliostat field reflective area";
	parameter Real SM = metadata_list[23] "Solar multiple";

	parameter SI.Efficiency ab_rec = 0.93 "Receiver coating absorptance";
	parameter SI.Efficiency em_rec = 0.87 "Receiver coating emissivity";
	parameter SI.Diameter D_receiver = metadata_list[4] "Receiver diameter";
	parameter SI.Length H_receiver = metadata_list[5] "Receiver height";
	parameter SI.HeatFlowRate Q_rec_out_des = Q_flow_des * SM "Heat from receiver at design";
	parameter Real[4] CL = {metadata_list[8],metadata_list[9],metadata_list[10],metadata_list[11]};
	parameter Real[4] C4L = {metadata_list[12],metadata_list[13],metadata_list[14],metadata_list[15]};
	parameter Real[5] CH = {metadata_list[16],metadata_list[17],metadata_list[18],metadata_list[19],metadata_list[20]};
	parameter SI.Temperature T_cold_set_CS = from_degC(500) "Cold tank target temperature";

	parameter SI.HeatFlowRate Q_flow_des = P_gross/eff_blk "Heat to power block at design, By product of PB initialisation, regardless which PB model is chosen e.g CEA or SAM";
	parameter Real par_fr = 0.1 "Parasitics fraction of power block rating at design point";
	parameter SI.Power P_gross = P_name / (1 - par_fr);

	parameter Integer horizon = 12 "Forecast horizon of the receiver dispatch algorithm";
	parameter Real dt = 300 "Forecast time step, in seconds";
	parameter Real const_t = -dt;
	parameter SI.Angle ele_min = from_deg(8) "Heliostat stow deploy angle";
	parameter Real nu_min = 0.25;
	parameter SI.HeatFlowRate Q_start = nu_min*Q_rec_out_des/eff_rec;
	parameter SI.HeatFlowRate Q_stop = nu_min*Q_rec_out_des/eff_rec;

	parameter SI.Area A_rec = pi*D_receiver*H_receiver "Area of receiver aperture";
	parameter SI.Efficiency eff_rec = metadata_list[7] "Receiver efficiency";
	parameter SI.Efficiency eff_blk = 0.3774 "Power block efficiency at design point";
	parameter SI.Power P_name = 19.9e6 "Nameplate rating of power block";
	parameter Real t_storage(unit="h") = 12 "Hours of storage";
	parameter SI.Energy E_max = P_name*t_storage*3600/eff_blk "Max stored energy";

	parameter SI.Energy E_up_u = 0.95*E_max "Upper energy limit";
	parameter SI.Energy E_up_l = 0.93*E_max "Upper energy limit";
	parameter SI.Energy E_low_u = 0.07*E_max "Lower energy limit";
	parameter SI.Energy E_low_l = 0.05*E_max "Lower energy limit";
	parameter SI.Irradiance dni_stop = 100 "DNI at which concentrator stops";
	parameter SI.Irradiance dni_start = 200 "DNI at which concentrator starts";

	parameter SI.Time t_con_up_min = 39*60 "Minimum operation time after concentrator starts";
	parameter SI.Time t_con_on_delay = 20*60 "Delay until concentrator starts";
	parameter SI.Time t_con_off_delay = 15*60 "Delay until concentrator shuts off";
	parameter SI.Time t_blk_on_delay = 15*60 "Delay until power block starts";
	parameter SI.Time t_blk_off_delay = 10*60 "Delay until power block shuts off";

	parameter Integer ramp_order = 1 "ramping filter order";

	parameter Integer n_sched_states = 1 "Number of schedule states";
	parameter Integer sch_state_start(min=1, max=n_sched_states) = 1 "Starting schedule state";
	parameter SI.Time t_sch_next_start = 0 "Time to next schedule change";
	parameter SI.HeatFlowRate Q_flow_sched_val[n_sched_states] = {
			P_name/eff_blk
			} "Heat flow at schedule states";
	parameter SI.Time t_delta[n_sched_states] = {
			24*3600
			} "Time differences between schedule states";

	parameter FI.Money C_field = 120*A_field "Field cost";
	parameter FI.Money C_site = 0 "Site improvements cost";
	parameter FI.Money C_tower = 0 "Tower cost";
	parameter FI.Money C_receiver = 135*A_rec "Receiver cost";
	parameter FI.Money C_storage = (30/(1e3*3600))*E_max "Storage cost";
	parameter FI.Money C_block = (1440/1e3)*P_name "Power block cost";
	parameter FI.Money C_bop = 0 "Balance of plant cost";
	parameter FI.Money C_land = 0 "Land cost";
	parameter FI.Money C_cap = C_field + C_site + C_tower + C_receiver + C_storage + C_block + C_bop + C_land "Capital costs";

	parameter FI.MoneyPerYear C_year =
			10*A_field // field cleaning/maintenance
			"Cost per year";
	parameter Real C_prod(unit="$/J/year") = 0 "Cost per production per year";
	parameter Real r_disc = 0.05 "Discount rate";
	parameter Integer t_life = 20 "Lifetime of plant";
	parameter Integer t_cons = 1 "Years of construction";

	parameter Boolean constrained = false "Constraint is present in optimisation if true";
		// If there is a constraint, then "constrained" must be a 'variable' Boolean
		// whose value is determined through an if statement with a constraint being the condition.
		// Note in this case the if block must be moved to the equation section.
	parameter Real distance = 0 "Distance to be added to a constant offset as added penalty when a constraint is not respected";
		// e.g. for euclidean distance: if T > T0 then constrained=true & distance=sqrt((T-T0)^2)
		// e.g. for quadratic distance: if T > T0 then constrained=true & distance=(T-T0)^2
		// e.g. for a constraint like T1 < T < T2, then T0 = (T1 + T2)/2
	parameter Real eps = Modelica.Constants.small;

	parameter String opticsTableNames[5] = {"optics_" + String(i-1) for i in 1:5};

	// Models
	SolarTherm.Models.CSP.CRS.HeliostatsField.Optical.TableMDBA oelts(
		angles = angles,
		file = opt_file,
		hra=sun.hra,
		dec=sun.dec,
		lat=lat,
		dni=sun.dni,
		ele=ele);

	// Models
	SolarTherm.Models.CSP.CRS.HeliostatsField.Optical.TableCosine W1(
		angles = angles,
		file = fileW1,
		hra=sun.hra,
		dec=sun.dec,
		lat=lat,
		dni=sun.dni,
		ele=ele);

	SolarTherm.Models.CSP.CRS.HeliostatsField.Optical.TableCosine W2(
		angles = angles,
		file = fileW2,
		hra=sun.hra,
		dec=sun.dec,
		lat=lat,
		dni=sun.dni,
		ele=ele);

	SolarTherm.Models.CSP.CRS.HeliostatsField.Optical.TableCosine W3(
		angles = angles,
		file = fileW3,
		hra=sun.hra,
		dec=sun.dec,
		lat=lat,
		dni=sun.dni,
		ele=ele);

	SolarTherm.Models.CSP.CRS.HeliostatsField.Optical.TableCosine W4(
		angles = angles,
		file = fileW4,
		hra=sun.hra,
		dec=sun.dec,
		lat=lat,
		dni=sun.dni,
		ele=ele);

	SolarTherm.Models.CSP.CRS.HeliostatsField.Optical.TableCosine W5(
		angles = angles,
		file = fileW5,
		hra=sun.hra,
		dec=sun.dec,
		lat=lat,
		dni=sun.dni,
		ele=ele);

	SolarTherm.Models.CSP.CRS.HeliostatsField.Optical.TableCosine W6(
		angles = angles,
		file = fileW6,
		hra=sun.hra,
		dec=sun.dec,
		lat=lat,
		dni=sun.dni,
		ele=ele);

	SolarTherm.Models.CSP.CRS.HeliostatsField.Optical.TableCosine W7(
		angles = angles,
		file = fileW7,
		hra=sun.hra,
		dec=sun.dec,
		lat=lat,
		dni=sun.dni,
		ele=ele);

	SolarTherm.Models.CSP.CRS.HeliostatsField.Optical.TableCosine W8(
		angles = angles,
		file = fileW8,
		hra=sun.hra,
		dec=sun.dec,
		lat=lat,
		dni=sun.dni,
		ele=ele);

	SolarTherm.Models.CSP.CRS.HeliostatsField.Optical.TableCosine W9(
		angles = angles,
		file = fileW9,
		hra=sun.hra,
		dec=sun.dec,
		lat=lat,
		dni=sun.dni,
		ele=ele);

	SolarTherm.Models.Sources.DataTable.DataTable data(
		lon = lon,
		lat = lat,
		t_zone = t_zone,
		year = year,
		file = wea_file) annotation(
			Placement(visible = true, transformation(extent = {{-132, -56}, {-102, -28}}, rotation = 0)));

	SolarTherm.Models.Sources.SolarModel.Sun sun(
		lon = lon,
		lat = lat,
		t_zone = t_zone,
		year = year,
		redeclare function solarPosition = SolarTherm.Models.Sources.SolarFunctions.PSA_Algorithm) annotation(
			Placement(transformation(extent = {{-82, 60}, {-62, 80}})));
			
	Modelica.Blocks.Sources.RealExpression dni_input(y = data.DNI);

	// Variables
	SI.HeatFlowRate Q_raw "Raw field output";
	SI.Angle ele "Elevation angle";

	FI.SpotPriceTable pri(file=pri_file);

	SI.HeatFlowRate Q_flow_chg "Heat flow into tank";
	SI.HeatFlowRate Q_flow_dis "Heat flow out of tank";
	SI.HeatFlowRate Q_ref;
	SI.HeatFlowRate Q_conv;
	SI.HeatFlowRate Q_rad;
	SI.HeatFlowRate Q_rec_out_raw;
	SI.Power P_elec "Output power of power block";

	Real fr_dfc(min=0, max=1) "Target energy fraction of the heliostat fistateld at the defocused state";
	Boolean full "True if the storage tank is full";

	Real fr_ramp_con (min=0.0, max=1.0) "ramping transition rate for the concentrator";
	SolarTherm.Utilities.Transition.Ramp ramp_up_blk(ramp_order=ramp_order, t_dur= t_blk_on_delay, up=true);
	SolarTherm.Utilities.Transition.Ramp ramp_down_blk(ramp_order=ramp_order, t_dur= t_blk_off_delay, up=false);
	Real fr_ramp_blk (min=0.0, max=1.0) "ramping transition rate for the power block";

	SI.Energy E(min=0, max=E_max) "Stored energy";
	SI.Energy E_rec "Energy dispatched from the receiver";
	SI.Energy E_rec_raw "Raw Energy dispatched from the receiver";

	SI.HeatFlowRate Q_flow_sched "Discharge schedule";

	Integer con_state(min=1, max=5) "Concentrator state";
	Integer blk_state(min=1, max=4) "Power block state";
	Integer sch_state(min=1, max=n_sched_states) "Schedule state";

	SI.HeatFlowRate Q_flow_rec_raw "Heat flow into receiver before curtailment";
	SI.HeatFlowRate Q_flow_rec "Heat flow into receiver";
	FI.Money R_spot(start=0, fixed=true) "Spot market revenue";
	SI.Energy E_elec(start=0, fixed=true) "Generate electricity";

	SI.Angle elo "Ecliptic longitude";
	SI.Angle azi "Ecliptic longitude";
	SI.Temperature T_ext_linear "Space average temperature for convective loss";
	SI.Temperature T_ext_4_linear "Fourth degree space average temperature for radiative loss";
	SI.CoefficientOfHeatTransfer h_conv "Heat transfer coefficient";

protected
	SI.Time  t_con_w_now "Time of concentrator current warm-up event";
	SI.Time  t_con_w_next "Time of concentrator next warm-up event";
	SI.Time  t_con_c_now "Time of concentrator current cool-down event";
	SI.Time  t_con_c_next "Time of concentrator next cool-down event";
	SI.Time  t_blk_w_now "Time of power block current warm-up event";
	SI.Time  t_blk_w_next "Time of power block next warm-up event";
	SI.Time  t_blk_c_now "Time of power block current cool-down event";
	SI.Time  t_blk_c_next "Time of power block next cool-down event";
	SI.Time  t_sch_next "Time of next schedule change";

function Phis "Compute coordinates along a line"
	input Real x     "Independent variable";
	output Real y    "Value of y at the specified x";
protected
	Real eps = Modelica.Constants.small;
algorithm
	if x<1 then
		y := exp(-1./(x+eps))/(exp(-1./(x+eps)) + exp(-1./(1-(x+eps))));
	else
		y := 1;
	end if;
end Phis;

initial equation
	E = E_low_l;
	Q_flow_sched = Q_flow_sched_val[sch_state_start];
	con_state = 1;
	blk_state = 1;
	sch_state = sch_state_start;
	t_con_w_now = 0;
	t_con_w_next = 0;
	t_con_c_now = 0;
	t_con_c_next = 0;
	t_blk_w_now = 0;
	t_blk_w_next = 0;
	t_blk_c_now = 0;
	t_blk_c_next = 0;
	t_sch_next = t_sch_next_start;

	if E > E_up_u then
		full = true;
	elseif E < E_up_l then
		full = false;
	else
		full = true;
	end if;

algorithm
	// Discrete equation system not yet supported (even though correct)
	// Putting in algorithm section instead
	when con_state == 2 and E >= E_up_u then
		con_state := 1; // off sun
	elsewhen con_state == 3 and (Q_raw <= Q_stop) and t_con_off_delay > 0 then
		con_state := 5; // ramp down
	elsewhen con_state == 3 and (Q_raw <= Q_stop) and t_con_off_delay <= 0 then
		con_state := 1; // off sun(no ramp-down)
	elsewhen con_state == 3 and full then
		con_state := 4; // on sun at part load
	elsewhen con_state == 4 and not full then
		con_state := 3; // on sun at full load
	elsewhen con_state == 4 and (Q_raw <= Q_stop) and t_con_off_delay > 0 then
		con_state := 5; // ramp down
	elsewhen con_state == 4 and (Q_raw <= Q_stop) and t_con_off_delay <= 0 then
		con_state := 1; // off sun (no ramp-down)
	elsewhen con_state == 1 and data.DNI >= dni_start and E <= E_up_l and t_con_on_delay > 0 then
		con_state := 2; // start onsteering (i.e. ramp up)
	elsewhen con_state == 1 and data.DNI >= dni_start and E <= E_up_l and t_con_on_delay <= 0 then
		con_state := 3; // on sun at full (no ramp-up)
	elsewhen con_state == 2 and time >= t_con_w_next then
		con_state := 3; // on sun at full load
	elsewhen con_state == 5 and time >= t_con_c_next then
		con_state := 1; // off sun
	end when;

	when blk_state == 2 and Q_flow_sched <= 0 then
		blk_state := 1; // turn off (or stop ramping) due to no demand
	elsewhen blk_state == 2 and E <= E_low_l then
		blk_state := 1; // turn off (or stop ramping) due to empty tank
	elsewhen blk_state == 3 and Q_flow_sched <= 0 and t_blk_off_delay > 0 then
		blk_state := 4; // ramp down due to no demand
	elsewhen blk_state == 3 and Q_flow_sched <= 0 and t_blk_off_delay <= 0 then
		blk_state := 1; // turn off (no ramp-down) due to no demand
	elsewhen blk_state == 3 and E <= E_low_l and t_blk_off_delay > 0 then
		blk_state := 4; // ramp down due to empty tank
	elsewhen blk_state == 3 and E <= E_low_l and t_blk_off_delay <= 0 then
		blk_state := 1; // turn off (no ramp down) due to empty tank
	elsewhen blk_state == 2 and time >= t_blk_w_next then
		blk_state := 3; // operational, ramp-up completed
	elsewhen blk_state == 1 and Q_flow_sched > 0 and E >= E_low_u  and t_blk_on_delay > 0 then
		blk_state := 2; // ramp up, demand and tank has capacity
	elsewhen blk_state == 1 and Q_flow_sched > 0 and E >= E_low_u  and t_blk_on_delay <= 0 then
		blk_state := 3; // operational (no ramp-up)
	elsewhen blk_state == 4 and time >= t_blk_c_next then
		blk_state := 1; // turn off after the ramp-down is complete
	end when;

	when time >= t_sch_next then
		sch_state := mod(pre(sch_state), n_sched_states) + 1;
	end when;

	when con_state == 2 then
		t_con_w_now := time;
		t_con_w_next := time + t_con_on_delay;
	end when;

	when con_state == 5 then
		t_con_c_now := time;
		t_con_c_next := time + t_con_off_delay;
	end when;

	when blk_state == 2 then
		t_blk_w_now := time;
		t_blk_w_next := time + t_blk_on_delay;
	end when;

	when blk_state == 4 then
		t_blk_c_now := time;
		t_blk_c_next := time + t_blk_off_delay;
	end when;

	when sch_state == 1 then
		Q_flow_sched := Q_flow_sched_val[1];
		t_sch_next := time + t_delta[1];
	end when;

	when E > E_up_u then
		full := true;
	elsewhen E < E_up_l then
		full := false;
	end when;

	if con_state == 2 then
		fr_ramp_con := Phis(((time - t_con_w_now) / t_con_on_delay));
	elseif con_state == 5 then
		fr_ramp_con := 1 - Phis(((time - t_con_c_now) / t_con_off_delay));
	else
		fr_ramp_con := 0;
	end if;


	if blk_state == 2 then
		fr_ramp_blk := if ramp_order == 0 then 0.0 else abs(ramp_up_blk.y);
	elseif blk_state == 4 then
		fr_ramp_blk := if ramp_order == 0 then 0.0 else abs(ramp_down_blk.y);
	else
	 fr_ramp_blk := 0;
	end if;

equation
	connect(sun.dni, dni_input.y);
	elo=SolarTherm.Models.Sources.SolarFunctions.eclipticLongitude(sun.dec);
	ele=SolarTherm.Models.Sources.SolarFunctions.elevationAngle(sun.dec,sun.hra,lat);
	azi=SolarTherm.Models.Sources.SolarFunctions.solarAzimuth(sun.dec,sun.hra,lat);
	Q_raw = data.DNI*oelts.nu*A_field;
	h_conv = CH[5] + CH[4]*data.Wspd + CH[3]*data.Wspd^2 + CH[2]*data.Wspd^3 + CH[1]*data.Wspd^4;
	T_ext_linear = CL[1] + CL[2]*(max(1,Q_flow_rec_raw))/1e6 + CL[3]*(max(1,data.Tdry)) + CL[4]*(max(1,data.Wspd));
	T_ext_4_linear = C4L[1] + C4L[2]*(max(1,Q_flow_rec_raw))/1e6 + C4L[3]*(max(1,data.Tdry)) + C4L[4]*(max(1,data.Wspd));
	Q_ref = (1-ab_rec)*(max(1,Q_flow_rec_raw));
	Q_rad = em_rec*sigma*H_receiver*H_receiver*pi*(T_ext_4_linear^4 - (max(1,data.Tdry))^4);
	Q_conv = h_conv*H_receiver*H_receiver*pi*(T_ext_linear-(max(1,data.Tdry)));
	Q_rec_out_raw = max(0, Q_flow_rec_raw - Q_ref - Q_rad - Q_conv);

	ramp_up_blk.x = t_blk_w_now;
	ramp_down_blk.x = t_blk_c_now;

	Q_flow_chg = eff_rec*Q_flow_rec;

	der(E) = Q_flow_chg - Q_flow_dis;
	der(E_rec) = Q_flow_rec;
	der(E_rec_raw) = Q_flow_rec_raw;

	if con_state <= 1 then
		Q_flow_rec = 0;
		Q_flow_rec_raw = 0;
		fr_dfc = 0;
	elseif con_state == 2 then
		Q_flow_rec = fr_ramp_con * oelts.nu*data.DNI*A_field;
		Q_flow_rec_raw = fr_ramp_con * oelts.nu*data.DNI*A_field;
		fr_dfc = if ramp_order == 0 then 0 else 1;
	elseif con_state == 5 then
		Q_flow_rec = fr_ramp_con * oelts.nu*data.DNI*A_field;
		Q_flow_rec_raw = fr_ramp_con * oelts.nu*data.DNI*A_field;
		fr_dfc = if ramp_order == 0 then 0 else 1;
	else
		if full then
			if eff_rec*oelts.nu*data.DNI*A_field > Q_flow_dis then
				Q_flow_rec = min(Q_flow_dis/eff_rec, oelts.nu*data.DNI*A_field);
				Q_flow_rec_raw = oelts.nu*data.DNI*A_field;
				fr_dfc = Q_flow_dis / (oelts.nu*data.DNI*A_field + eps);
			else
				Q_flow_rec = oelts.nu*data.DNI*A_field;
				Q_flow_rec_raw = oelts.nu*data.DNI*A_field;
				fr_dfc = 1;
			end if;
		else
			Q_flow_rec = oelts.nu*data.DNI*A_field;
			Q_flow_rec_raw = oelts.nu*data.DNI*A_field;
			fr_dfc = 1;
		end if;
	end if;

	if blk_state <=1 then
		Q_flow_dis = 0;
		P_elec = 0;
	elseif blk_state == 2 then
		Q_flow_dis = if ramp_order == 0 then Q_flow_sched else fr_ramp_blk * Q_flow_sched;
		P_elec = eff_blk*Q_flow_dis;
	elseif blk_state == 4 then
		Q_flow_dis = fr_ramp_blk * Q_flow_sched;
		P_elec = eff_blk*Q_flow_dis;
	else
		Q_flow_dis = Q_flow_sched;
		P_elec = eff_blk*Q_flow_dis;
	end if;

	der(E_elec) = P_elec;
	der(R_spot) = P_elec*pri.price;
	annotation(experiment(StartTime=0.0, StopTime=31536000.0, Interval=300, Tolerance=1e-06));
end FieldAnnualOpticalInterpolation;
within examples;
model SimpleSystemOperation "High temperature Sodium-sCO2 system"
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

	parameter nSI.Angle_deg lon = -5.326 "Longitude (+ve East)";
	parameter nSI.Angle_deg lat = 37.558 "Latitude (+ve North)";
	parameter nSI.Time_hour t_zone = 1 "Local time zone (UCT=0)";
	parameter Integer year = 2006 "Meteorological year";

	parameter SolarTherm.Types.Solar_angles angles = SolarTherm.Types.Solar_angles.elo_hra "Angles used in the lookup table file";
	parameter String opt_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Data/salt_N06230_OD22.40_WT1.20_565/OELTs_Solstice.motab");
	parameter Real metadata_list[23] = metadata(opt_file);
	parameter Integer n_heliostat = integer(metadata_list[1]) "Number of heliostats";
	parameter SI.Area A_heliostat = metadata_list[2] "Heliostat module reflective area";
	parameter SI.Area A_field = n_heliostat * A_heliostat "Heliostat field reflective area";
	parameter Real SM = metadata_list[23] "Solar multiple";

	parameter SI.Efficiency ab_rec = 0.98 "Receiver coating absorptance";
	parameter SI.Efficiency em_rec = 0.91 "Receiver coating emissivity";
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

	parameter FI.MoneyPerYear C_year = 10*A_field "Cost per year";
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

	// Flux interpolation parameters
	parameter Integer N = 450 "Number of tube segments in flowpath";
	parameter String opticsTableNames[5] = {"optics_" + String(i-1) for i in 1:5};
	parameter String tableNames[N] = {"flux_" + String(i) for i in 1:N};
	parameter String tablemflowNames[5] = {"mflow_" + String(i) for i in 1:5}; //1:0.56 2:0.87 3:1.0 4:1.20 5:1.39
	parameter String file_dni3 = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Data/salt_N06230_OD22.40_WT1.20_565/FLUX_fp1_d1.0.motab");
	parameter String file_mflow = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Data/salt_N06230_OD22.40_WT1.20_565/MFLOW_Solstice_fp1.motab");

	SolarTherm.Models.Sources.DataTable.DataTable data(
		lon = lon,
		lat = lat,
		t_zone = t_zone,
		year = year,
		file = wea_file) annotation(
			Placement(visible = true, transformation(extent = {{-132, -56}, {-102, -28}}, rotation = 0)));

	Modelica.Blocks.Sources.RealExpression u1(y = to_deg(elo));
	Modelica.Blocks.Sources.RealExpression u2(y = to_deg(sun.hra));

	SolarTherm.Models.Sources.SolarModel.Sun sun(lat = lat, lon = lon, t_zone = t_zone, year = year);
	SI.CoefficientOfHeatTransfer h_conv "Heat transfer coefficient";
	Integer con_state(min=1, max=5) "Concentrator state";

	FI.SpotPriceTable pri(file=pri_file);
	SI.Angle elo "Ecliptic longitude";
	SI.Angle ele "Elevation angle";
	SI.Angle azi "Azimuth angle";
	SI.Power P_elec "Output power of power block";

	FI.Money R_spot(start=0, fixed=true) "Spot market revenue";
	SI.Energy E_elec(start=0, fixed=true) "Generate electricity";

	SI.HeatFlux CG[N] "Interpolated incident flux";
	SI.MassFlowRate m_flow_tb "Interpolated mass flow rate";

	SI.HeatFlowRate Q_raw "Raw field output";

protected
	SolarTherm.Models.CSP.CRS.HeliostatsField.Optical.TableMDBA oelts(
		angles = angles,
		file = opt_file,
		hra=sun.hra,
		dec=sun.dec,
		lat=lat,
		dni=sun.dni,
		ele=ele);

	Modelica.Blocks.Tables.CombiTable2D m_flow[5](
		each fileName = file_mflow, 
		each tableOnFile = true, 
		each smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
		tableName = tablemflowNames);

	Modelica.Blocks.Tables.CombiTable2D flux_dni3[N](
		each fileName = file_dni3, 
		each tableOnFile = true, 
		each smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
		tableName = tableNames);

equation
	Q_raw = sun.dni*oelts.nu*A_field;
	con_state = 1;
	sun.dni=950.2*(1-exp(-0.075*to_deg(ele)));
	h_conv = CH[5] + CH[4]*data.Wspd + CH[3]*data.Wspd^2 + CH[2]*data.Wspd^3 + CH[1]*data.Wspd^4;
	elo=SolarTherm.Models.Sources.SolarFunctions.eclipticLongitude(sun.dec);
	// Flux interpolation
	for i in 1:N loop
		connect(u1.y, flux_dni3[i].u1);
		connect(u2.y, flux_dni3[i].u2);
	end for;
	// Mass flow rate interpolation
	for i in 1:5 loop
		connect(u1.y, m_flow[i].u1);
		connect(u2.y, m_flow[i].u2);
	end for;
	ele = SolarTherm.Models.Sources.SolarFunctions.elevationAngle(sun.dec,sun.hra,lat);
	azi = SolarTherm.Models.Sources.SolarFunctions.solarAzimuth(sun.dec,sun.hra,lat);

	for i in 1:N loop
		CG[i] = max(0, flux_dni3[i].y);
	end for;
	m_flow_tb = max(0, m_flow[3].y);

	P_elec = P_name;

	der(E_elec) = P_elec;
	der(R_spot) = P_elec*pri.price;
	annotation(experiment(StartTime=0.0, StopTime=31536000.0, Interval=300, Tolerance=1e-06));
end SimpleSystemOperation;

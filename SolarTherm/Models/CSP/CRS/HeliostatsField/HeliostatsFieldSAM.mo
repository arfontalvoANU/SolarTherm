within SolarTherm.Models.CSP.CRS.HeliostatsField;
model HeliostatsFieldSAM
  extends Interfaces.Models.Heliostats;
  parameter nSI.Angle_deg lon=133.889 "Longitude (+ve East)" annotation(Dialog(group="System location"));
  parameter nSI.Angle_deg lat=-23.795 "Latitude (+ve North)" annotation(Dialog(group="System location"));
  parameter Integer n_h=1 "Number of heliostats" annotation(Dialog(group="Technical data"));
  parameter SI.Area A_h=4 "Heliostat's Area" annotation(Dialog(group="Technical data"));
  parameter Real he_av=0.99 "Heliostat availability" annotation(Dialog(group="Technical data"));
  replaceable model Optical =
      SolarTherm.Models.CSP.CRS.HeliostatsField.Optical.Constant
  constrainedby
    SolarTherm.Models.CSP.CRS.HeliostatsField.Optical.OpticalEfficiency
    "Total optical efficency"
    annotation (Dialog(group="Technical data",__Dymola_label="nu"), choicesAllMatching=true);
  parameter Boolean use_on = false
    "= true to display when solar field is connected"
      annotation (Dialog(group="Operating strategy"), Evaluate=true, HideResult=true, choices(checkBox=true));
  parameter Boolean use_defocus = false "= true to use defocus strategy"
      annotation (Dialog(group="Operating strategy"), Evaluate=true, HideResult=true, choices(checkBox=true));
  parameter Boolean use_wind = false "= true to use Wind-stop strategy"
      annotation (Dialog(group="Operating strategy"), Evaluate=true, HideResult=true, choices(checkBox=true));
  parameter SI.Angle ele_min=from_deg(8) "Heliostat stow deploy angle" annotation(min=0,Dialog(group="Operating strategy"));
  parameter SI.HeatFlowRate Q_design=529.412 "Receiver design thermal power (with heat losses)" annotation(min=0,Dialog(group="Operating strategy"));
  parameter Real nu_start=0.60 "Receiver energy start-up fraction" annotation(min=0,Dialog(group="Operating strategy"));
  parameter Real nu_min=0.25 "Minimum receiver turndown energy fraction" annotation(min=0,Dialog(group="Operating strategy"));
  parameter Real nu_defocus=1 "Receiver limiter energy fraction at defocus state" annotation(Dialog(group="Operating strategy",enable=use_defocus));
  parameter SI.Velocity Wspd_max=15 "Wind stow speed" annotation(min=0,Dialog(group="Operating strategy",enable=use_wind));

  parameter SI.Energy E_start=90e3 "Start-up energy of a single heliostat" annotation(Dialog(group="Parasitic loads"));
  parameter SI.Power W_track=0.055e3 "Tracking power for a single heliostat" annotation(Dialog(group="Parasitic loads"));
  parameter SI.Time t_startup=720 "Receiver startup delay time";

  Optical optical(hra=solar.hra, dec=solar.dec, lat=lat);
  SI.HeatFlowRate Q_raw;
  SI.HeatFlowRate Q_net;

  Modelica.Blocks.Interfaces.BooleanOutput on if use_on annotation (Placement(
        transformation(extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={0,-100}),                           iconTransformation(extent={{-12,-12},
            {12,12}},
        rotation=-90,
        origin={0,-100})));
  Modelica.Blocks.Interfaces.BooleanInput defocus if use_defocus annotation (Placement(
        transformation(extent={{-126,-88},{-86,-48}}),iconTransformation(extent={{-110,
            -72},{-86,-48}})));

  Modelica.Blocks.Interfaces.RealInput Wspd if use_wind annotation (Placement(
        transformation(extent={{-126,50},{-86,90}}), iconTransformation(extent={
            {-110,50},{-86,74}})));

  SI.Time t_ramp_start(start = 0.0) "current time of simulation when ramping starts";
  SI.Time t_ramp_end(start = t_startup) "current time of simulation when ramping starts";
  Real rec_ramp_fraction;

  SI.Angle elo;
  SI.Angle ele;
  SI.Angle zen;
  SI.Angle zen2;
  SI.Angle azi;
  SI.Energy E_raw;
  SI.Energy E_dni;
  SI.Energy E_field;

  SI.Power W_loss;
  Real damping;
//protected
  Boolean on_hf;
protected
  SI.Power W_loss1;
  SI.Power W_loss2;
  //SI.Time t_start=30*60;
  parameter SI.Time t_start=3600 "Start-up traking delay";
  discrete Modelica.SIunits.Time t_on(start=0, fixed=true) "Sunrise time instant";
  Modelica.Blocks.Interfaces.BooleanInput on_internal
    "Needed to connect to conditional connector";
  Modelica.Blocks.Interfaces.BooleanInput defocus_internal
    "Needed to connect to conditional connector";
  Modelica.Blocks.Interfaces.RealInput Wspd_internal
    "Needed to connect to conditional connector";
  parameter SI.HeatFlowRate Q_start=nu_start*Q_design "Heliostat field start power" annotation(min=0,Dialog(group="Operating strategy"));
  parameter SI.HeatFlowRate Q_min=nu_min*Q_design "Heliostat field turndown power" annotation(min=0,Dialog(group="Operating strategy"));
  parameter SI.HeatFlowRate Q_defocus=nu_defocus*Q_design "Heat flow rate limiter at defocus state" annotation(Dialog(group="Operating strategy",enable=use_defocus));
  parameter SI.Energy E_startup=nu_start*t_startup*Q_design "Minimum energy required to startup the receiver during the delay time";
algorithm

	when E_raw>E_startup then//when Q_raw>Q_start then
		if on_internal == false then
		//on_internal:=true;
			t_ramp_start := time;
			t_ramp_end := time + t_startup;
		end if;
		on_internal:=true;
	elsewhen Q_raw<Q_min then
	on_internal:=false;
	end when;
initial equation
   on_internal=E_raw>E_startup;
//   on_internal=Q_raw>Q_start;
equation
  if use_on then
    connect(on,on_internal);
  end if;
  if use_defocus then
    connect(defocus,defocus_internal);
  else
    defocus_internal = false;
  end if;
    if use_wind then
    connect(Wspd,Wspd_internal);
  else
    Wspd_internal = -1;
  end if;

  on_hf=(ele>ele_min) and (Wspd_internal<Wspd_max);

  E_raw= if on_hf then max(he_av*n_h*A_h*solar.dni*optical.nu,0)*t_startup*0.5 else 0;

  Q_raw= if on_hf then max(he_av*n_h*A_h*solar.dni*optical.nu,0) else 0;

  if rec_ramp_fraction <1.0-1e-6 then
    Q_net= if on_internal then (if defocus_internal then min(Q_defocus,(Q_raw - Q_start)) else (Q_raw - Q_start)) else 0;
  else
    Q_net= if on_internal then (if defocus_internal then min(Q_defocus,Q_raw) else Q_raw) else 0;
  end if;

  rec_ramp_fraction = min(1.0, (time - t_ramp_start) / (t_ramp_end - t_ramp_start));

  heat.Q_flow= -Q_net;
  elo=SolarTherm.Models.Sources.SolarFunctions.eclipticLongitude(solar.dec);
//   optical.hra=solar.hra;
//   optical.dec=solar.dec;

  ele=SolarTherm.Models.Sources.SolarFunctions.elevationAngle(
    solar.dec,
    solar.hra,
    lat);
  zen=pi/2-ele;
  zen2=SolarTherm.Models.Sources.SolarFunctions.aparentSolarZenith(
    solar.dec,
    solar.hra,
    lat);
  azi=SolarTherm.Models.Sources.SolarFunctions.solarAzimuth(
    solar.dec,
    solar.hra,
    lat);

  der(E_field) = Q_net;
  der(E_dni) = he_av*n_h*A_h*solar.dni;
  damping= if on_internal then Q_net/(Q_raw+1e-3) else 1;
  W_loss1=if ele>1e-2 then n_h*he_av*damping*W_track else 0;
  when ele>1e-2 then
    t_on=time;
  end when;
  W_loss2= if time<t_on+t_start then n_h*he_av*damping*E_start/t_start else 0;
  W_loss=W_loss1+W_loss2;
  annotation (Documentation(info = "<html>
</html>", revisions = "<html>
<ul>
<li>Alberto de la Calle:<br>Released first version. </li>
</ul>
</html>"),
    Icon(graphics = {Text(origin = {0, -140}, lineColor = {0, 0, 255}, extent = {{-140, 20}, {140, -20}}, textString = "%name")}));
end HeliostatsFieldSAM;

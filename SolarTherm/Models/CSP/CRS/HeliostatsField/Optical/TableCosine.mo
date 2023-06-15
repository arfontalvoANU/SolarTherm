within SolarTherm.Models.CSP.CRS.HeliostatsField.Optical;
model TableCosine "From table"
  extends OpticalEfficiency;
  import SolarTherm.Utilities.FluxInterpolation;
  parameter String file "File where optical data matrix is stored" annotation (Dialog(
      group="Technical data",
      enable=tableOnFile,
      loadSelector(filter="TMY3 custom-built files (*.motab);;MATLAB MAT-files (*.mat)",
                   caption="Open file in which optical data is present")));
  parameter SolarTherm.Types.Solar_angles angles=SolarTherm.Types.Solar_angles.elo_hra
    "Table angles"
    annotation (Dialog(group="Table data interpretation"));

  parameter String opticsTableNames[5] = {"optics_" + String(i-1) for i in 1:5};

  parameter SI.Angle ele_min=Modelica.SIunits.Conversions.from_deg(8) "Heliostat stow deploy angle"
    annotation(min=0,Dialog(group="Operating strategy"));

  SI.Angle angle1;
  SI.Angle angle2;
  input SI.Angle ele;
  input SI.Irradiance dni;

  Modelica.Blocks.Tables.CombiTable2D nu_table[5](
    each fileName = file, 
    each tableOnFile = true, 
    each smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
    tableName = opticsTableNames);

  Modelica.Blocks.Sources.RealExpression angle2_input(y=to_deg(angle2));
  Modelica.Blocks.Sources.RealExpression angle1_input(y=to_deg(angle1));
  Real loss;

equation

  if angles==SolarTherm.Types.Solar_angles.elo_hra then
    angle1=SolarTherm.Models.Sources.SolarFunctions.eclipticLongitude(dec);
    angle2=hra;
  elseif angles==SolarTherm.Types.Solar_angles.dec_hra then
    angle1=dec;
    angle2=hra;
  elseif angles==SolarTherm.Types.Solar_angles.ele_azi then
    angle1=SolarTherm.Models.Sources.SolarFunctions.elevationAngle(dec,hra,lat);
    angle2=SolarTherm.Models.Sources.SolarFunctions.solarAzimuth(dec,hra,lat);
  else
    angle1=SolarTherm.Models.Sources.SolarFunctions.solarZenith(dec,hra,lat);
    angle2=SolarTherm.Models.Sources.SolarFunctions.solarAzimuth(dec,hra,lat);
  end if;
  nu = 1;

  loss = min(1,max(0, FluxInterpolation(nu_table[1].y, nu_table[2].y, nu_table[3].y, nu_table[4].y, nu_table[5].y, ele, dni, ele_min)));
//  loss = max(0, nu_table[3].y);

  for i in 1:5 loop
    connect(angle2_input.y, nu_table[i].u2);
    connect(angle1_input.y, nu_table[i].u1);
  end for;

end TableCosine;
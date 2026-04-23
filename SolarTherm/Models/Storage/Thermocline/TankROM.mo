within SolarTherm.Models.Storage.Thermocline;

model TankROM "TES Component model of a single thermocline tank with spherical fillers"
  extends SolarTherm.Interfaces.Models.StorageFluid_Thermocline;
  extends SolarTherm.Icons.PackedBedSpheres;
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  //Initialise Material Packages
  replaceable package Medium = SolarTherm.Media.Sodium.Sodium_pT;
  replaceable package Fluid_Package = SolarTherm.Materials.PartialMaterial;
  replaceable package Filler_Package = SolarTherm.Materials.PartialMaterial;
  // Parameters
  parameter Integer Correlation = 1 "Interfacial convection correlation {1 = WakaoKaguei, 2 = MelissariArgyropoulos, 3 = Conservative}";
  parameter SI.Energy E_max = 144.0e9 "Maximum storage capacity";
  parameter SI.Height H_tank = 1.2;
  parameter SI.Diameter D_tank = 0.148;
  parameter Real epsilon = 0.4 "Volume of tank occupied by the heat transfer fluid";
  parameter SI.Length ds = 0.02 "Filler sphere diameter";
  parameter Integer Nz = 100 "Number of vertical elements";
  parameter SI.CoefficientOfHeatTransfer U_wall = 0.1 "Overall heat loss coefficient through tank walls (W/m2/K)";
  parameter SI.Temperature T_min = 293 "HTF discharging temperature at design";
  parameter SI.Temperature T_max = 823 "HTF charging temperature at design";
  parameter SI.Temperature T_start = 293 "Initial (uniform) temperature of all components, in Kelvin, defaults to T_min";
  final parameter Real ZDH[Nz] = Relative_Tank_Axes(H_tank, Nz) "Non-dimensional tank vertical axis";
  //Input and Output Ports
  Modelica.Blocks.Interfaces.RealOutput T_top_measured "Temperature at the top of the tank as an output signal (K)";
  Modelica.Blocks.Interfaces.RealOutput T_bot_measured "Temperature at the bottom of the tank as an output signal (K)";
  Modelica.Blocks.Interfaces.RealOutput T_p_top_measured "Temperature of the innermost solid element at the the hot-end of the TES (K)";
  Modelica.Blocks.Interfaces.RealOutput T_p_bot_measured "Temperature of the innermost solid element at the the cold-end of the TES (K)";
  Modelica.Blocks.Interfaces.RealOutput h_bot_outlet "Enthaply at the bottom of the tank as an output signal (J/kg)";
  Modelica.Blocks.Interfaces.RealOutput h_top_outlet "Enthaply at the top of the tank as an output signal (J/kg)";
  Modelica.Blocks.Interfaces.RealInput T_amb "Ambient Temperature" annotation(
    Placement(visible = true, transformation(origin = {-50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-46, 40}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput p_amb "Ambient Pressure" annotation(
    Placement(visible = true, transformation(origin = {48, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0), iconTransformation(origin = {-46, -40}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput Level "Theoretical Tank Level" annotation(
    Placement(visible = true, transformation(extent = {{40, 16}, {60, 36}}, rotation = 0), iconTransformation(origin = {45, 37.5}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  //Analysis of fluid entering and exiting storage
  Fluid_Package.State fluid_top "Fluid entering/exiting top";
  Fluid_Package.State fluid_bot "Fluid entering/exiting bottom";
  SI.Pressure p_drop_total "Sum of all pressure drops (Pa)";
  SI.Power W_dot_loss_pump "losses due to pressure drop (W)";
  SI.Energy E;
  SI.HeatFlowRate Q_TES_in, Q_loss;
  Real e_top;
  Real e_bot;
  //SI.HeatFlowRate Q_loss_TES = (f_min + Level*(f_max - f_min))*E_max;
equation
  e_bot = 0.741704 + 0.261359/(1 + exp(+15.317555*(Level - 0.646759)));
  e_top = 0.723180 + 0.280802/(1 + exp(-14.901512*(Level - 0.352074))); 
  fluid_a.m_flow + fluid_b.m_flow = 0;
  if fluid_a.m_flow > 0 then //HTF flowing from top to bottom
    fluid_top.h = inStream(fluid_a.h_outflow);
    h_top_outlet = inStream(fluid_a.h_outflow);
    T_top_measured = fluid_top.T;
    T_bot_measured = T_max - e_bot*(T_max-T_min);
    fluid_bot.T = T_bot_measured;
    fluid_bot.h = h_bot_outlet;
    fluid_a.h_outflow = h_bot_outlet;
    fluid_b.h_outflow = h_bot_outlet;
    Q_loss = 16523.4053782364*Level+9663.97237617296;
  else //discharging
    fluid_bot.h = inStream(fluid_b.h_outflow);
    h_bot_outlet = inStream(fluid_b.h_outflow);
    T_bot_measured = fluid_bot.T;
    T_top_measured = T_min + e_top*(T_max-T_min);
    fluid_top.T = T_top_measured;
    fluid_top.h = h_top_outlet;
    fluid_a.h_outflow = h_top_outlet;
    fluid_b.h_outflow = h_top_outlet;
    Q_loss = 16671.4488706001*Level+9199.98246747617;
  end if;
  der(E) = Q_TES_in - Q_loss;
  Level = E/E_max;
  fluid_a.p = p_amb;
  fluid_b.p = p_amb;
  p_drop_total = 0.0;
  W_dot_loss_pump = 0.0;
  T_p_top_measured = fluid_top.T;
  T_p_bot_measured = fluid_bot.T;
  annotation(
    Documentation(info = "<html>
		<p>This model contains the fluid_a (top) and fluid_b (bottom) ports, basically a complete CSP component. 
		This model simply connects the Thermocline_Spheres_Section models to the correct ports.</p>
		</html>", revisions = "<html>
        <ul>
		<li><i>Dec 2020</i> by Zebedee Kee:<br>
		Released first version.</li>
		<li><i>Feb 2026</i> by Armando Fontalvo:<br>
		Simplification for verification purposes.</li>
		</ul>
		</html>"));
end TankROM;
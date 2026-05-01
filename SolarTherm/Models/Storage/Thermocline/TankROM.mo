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
  //************************* NEW STUFF ***************************
  Modelica.SIunits.SpecificEnthalpy h_in "Enthalpy at inlet";
  Modelica.SIunits.SpecificEnthalpy h_out "Enthalpy at outlet";
  Modelica.SIunits.SpecificEnthalpy h_top "Enthalpy at TES top";
  Modelica.SIunits.SpecificEnthalpy h_bot "Enthalpy at TES bottom";
  final parameter SI.Density rhof_avg = (Fluid_Package.rho_Tf(T_min,0) + Fluid_Package.rho_Tf(T_max,0)) / 2;
  final parameter Modelica.SIunits.Area A_tank = 0.25*Modelica.Constants.pi*D_tank^2;
  final parameter Modelica.SIunits.Mass m_tank = rhof_avg*H_tank*A_tank;
  final parameter Modelica.SIunits.Mass m_section = m_tank/Nz;
  // Charging e_bot
  parameter Real Cc = 644641.7434;
  parameter Real Lc = 182774.3657;
  parameter Real kc = 15.1010;
  parameter Real tc = 0.6568;
  // Discharging e_top
  parameter Real Cd = 1167096.5487;
  parameter Real Ld = 207697.2662;
  parameter Real kd = 15.0865;
  parameter Real td = 0.3507;
  parameter Modelica.SIunits.SpecificEnthalpy hf_start = Fluid_Package.h_Tf(T_start, 0.0);
  Fluid_Package.State state_top;
  Fluid_Package.State state_bot;

initial equation
  state_top.T = T_start;
  state_bot.T = T_start;
equation
  if fluid_a.m_flow > 1e-6 then
    fluid_top.h = inStream(fluid_a.h_outflow);
    fluid_bot.h = fluid_b.h_outflow;
  elseif fluid_a.m_flow < -1e-6 then
    fluid_top.h = fluid_a.h_outflow;
    fluid_bot.h = inStream(fluid_b.h_outflow);
  else
    fluid_top.T = 298.15;
    fluid_bot.T = 298.15;
  end if;
  
  //Calculate tank energy level
  der(E) = Q_TES_in - Q_loss;
  Level = E/E_max;
  
  //Determine tank outlet enthalpy used by external control system
  h_bot_outlet = h_bot;
  h_top_outlet = h_top;
  
  //Mass balance
  fluid_a.m_flow + fluid_b.m_flow = 0; //always true for a steady state component
  if fluid_a.m_flow > 0 then //mass is flowing into the top, direction is downwards so Tank_A.m_flow is (negative), charging
    h_in = inStream(fluid_a.h_outflow);
    m_section*der(h_top) = fluid_a.m_flow*(h_in - h_top);
    der(h_bot) = Lc*kc*exp(-kc*(Level - tc))/(1 + exp(-kc*(Level - tc)))^2*der(Level);
    h_out = h_bot;
    fluid_a.h_outflow = h_in;
    fluid_b.h_outflow = h_out;
    Q_loss = 16523.4053782364*Level+9663.97237617296;
  else //discharging
    h_in = inStream(fluid_b.h_outflow);
    m_section*der(h_bot) = fluid_b.m_flow*(h_in - h_bot);
    der(h_top) = Ld*kd*exp(-kd*(Level - td))/(1 + exp(-kd*(Level - td)))^2*der(Level);
    h_out = h_top;
    fluid_a.h_outflow = h_out;
    fluid_b.h_outflow = h_in;
    Q_loss = 16671.4488706001*Level+9199.98246747617;
  end if;

  fluid_a.p = p_amb;
  fluid_b.p = p_amb;
  //T_amb = Tank_A.T_amb;
  state_top.h = h_top;
  state_bot.h = h_bot;
  T_top_measured = state_top.T;
  T_bot_measured = state_bot.T;
  T_p_top_measured = state_top.T;
  T_p_bot_measured = state_bot.T;

  p_drop_total = 0.0;
  W_dot_loss_pump = 0.0;
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

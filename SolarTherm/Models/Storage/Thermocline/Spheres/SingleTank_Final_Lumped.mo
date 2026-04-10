within SolarTherm.Models.Storage.Thermocline.Spheres;

model SingleTank_Final_Lumped "TES Component model of a single thermocline tank with spherical fillers"
  extends SolarTherm.Interfaces.Models.StorageFluid_Thermocline;
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  //Initialise Material Packages
  replaceable package Medium = SolarTherm.Media.Sodium.Sodium_pT;
  replaceable package Fluid_Package = SolarTherm.Materials.PartialMaterial;
  replaceable package Filler_Package = SolarTherm.Materials.PartialMaterial;
  
  //Storage Parameter Settings
  parameter Integer Correlation = 1 "Interfacial convection correlation {1 = WakaoKaguei, 2 = MelissariArgyropoulos, 3 = Conservative}";

  //Storage CApacity
  parameter SI.Energy E_max = 144.0e9 "Maximum storage capacity";
  
  //Aspect ratios (H/D) of tank
  parameter SI.Height H_tank = 1.2;
  parameter SI.Diameter D_tank = 0.148;
  
  //Porosity of tank filler materials
  parameter Real epsilon = 0.4 "Porosity";
  
  //Filler diameter of materials
  parameter SI.Length ds = 0.02 "Filler sphere diameter";
  
  //Discretization Settings
  parameter Integer Nz = 10;
  
  //Heat loss coefficient of tanks
  parameter SI.CoefficientOfHeatTransfer U_wall = 0.1 "W/m2K";
  
  //Temperature Settings
  parameter SI.Temperature T_min = 293 "Minimum temperature (design) also starting T";
  parameter SI.Temperature T_max = 823 "Maximum design temperature (design)";
  parameter SI.Temperature T_start = 293 "Initial (uniform) temperature of all components (K), defaults to T_min";

  //Input and Output Ports
  Modelica.Blocks.Interfaces.RealOutput T_top_measured "Temperature at the top of the tank as an output signal (K)"
    annotation (Placement(visible = true,
      transformation(extent = {{40, 50}, {60, 70}}, rotation = 0),
      iconTransformation(origin = {45, 55}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));

  Modelica.Blocks.Interfaces.RealOutput T_bot_measured "Temperature at the bottom of the tank as an output signal (K)"
    annotation (Placement(visible = true,
      transformation(extent = {{40, -70}, {60, -50}}, rotation = 0), 
      iconTransformation(origin = {45, -57}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
          
  Modelica.Blocks.Interfaces.RealOutput T_p_top_measured = Tank_A.Ts[Nz] "Temperature of the innermost solid element at the the hot-end of the TES (K)"
    annotation (Placement(visible = true,
      transformation(extent = {{40, 36}, {60, 56}}, rotation = 0), 
      iconTransformation(origin = {45, 43}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));

  Modelica.Blocks.Interfaces.RealOutput T_p_bot_measured = Tank_A.Ts[1] "Temperature of the innermost solid element at the the cold-end of the TES (K)"
    annotation (Placement(visible = true,
      transformation(extent = {{40, -54}, {60, -34}}, rotation = 0), 
      iconTransformation(origin = {45, -45}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  
  Modelica.Blocks.Interfaces.RealOutput h_bot_outlet "Enthaply at the bottom of the tank as an output signal (J/kg)"
    annotation (Placement(visible = true,
      transformation(origin = {-40, -70},extent = {{-10, -10}, {10, 10}}, rotation = -90),
      iconTransformation(origin = {-27, -69}, extent = {{-5, -5}, {5, 5}}, rotation = -90)));
          
  Modelica.Blocks.Interfaces.RealOutput h_top_outlet "Enthaply at the top of the tank as an output signal (J/kg)"
    annotation (Placement(visible = true,
      transformation(origin = {-40, 56},extent = {{10, -10}, {-10, 10}}, rotation = -90),
      iconTransformation(origin = {-27, 65}, extent = {{5, -5}, {-5, 5}}, rotation = -90)));

  Modelica.Blocks.Interfaces.RealInput T_amb "Ambient Temperature"
    annotation (
        Placement(
            visible = true,
            transformation(
                origin={-50, 0},
                extent={{-10, -10}, {10, 10}},
                rotation=0), 
            iconTransformation(
                origin={-46, 40},
                extent={{-6, -6}, {6, 6}},
                rotation=0)));
  Modelica.Blocks.Interfaces.RealInput p_amb "Ambient Pressure" 
    annotation (
        Placement(
            visible = true,
            transformation(
                origin={48, 0},
                extent={{10, -10}, {-10, 10}},
                rotation=0),
            iconTransformation(
                origin={-46, -40},
                extent={{-6, -6}, {6, 6}},
                rotation=0)));
  
  //Initialize Tank
  SolarTherm.Models.Storage.Thermocline.Spheres.Section_Final_Lumped Tank_A(
    redeclare replaceable package Fluid_Package = Fluid_Package,
    redeclare replaceable package Filler_Package = Filler_Package,
    Correlation = Correlation,
    E_max = E_max,
    epsilon = epsilon,
    ds = ds,
    T_min = T_min,
    T_max = T_max,
    T_start = T_start,
    Nz = Nz,
    U_wall = U_wall,
    H_tank = H_tank,
    D_tank = D_tank) 
    annotation (
        Placement(
            visible = true,
            transformation(
                origin={0, 0},
                extent={{10, -10}, {-10, 10}},
                rotation=0),
            iconTransformation(
                origin={0, 0},
                extent={{-10, -10}, {10, 10}},
                rotation=0)));

  Modelica.Blocks.Interfaces.RealOutput Level "Theoretical Tank Level"
    annotation (Placement(visible = true,
      transformation(extent = {{40, 16}, {60, 36}}, rotation = 0),
      iconTransformation(origin = {45, 21}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));

  //Tank Non-dimensionalized vertical axis
  final parameter Real ZDH[Nz] = Tank_A.ZDH;
  
  //Plotting Temperature degC
  Real Tf_degC[Nz](start=fill(T_min,Nz));
  
  //Analysis of fluid entering and exiting storage
  Fluid_Package.State fluid_top "Fluid entering/exiting top";
  Fluid_Package.State fluid_bot "Fluid entering/exiting bottom";
  
  SI.Pressure p_drop_total "Sum of all pressure drops (Pa)";
  SI.Power W_dot_loss_pump "losses due to pressure drop (W)";

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
  //Convert from Kelvin to degC for easier plotting
  Tf_degC = (Tank_A.Tf).-273.15;
  
  //Calculate tank energy level
  Level = Tank_A.Level;
  
  //Determine tank outlet enthalpy used by external control system
  h_bot_outlet = Tank_A.hf[1];
  h_top_outlet = Tank_A.hf[Nz];
  
  //Mass balance
  fluid_a.m_flow = -1.0*fluid_b.m_flow; //always true for a steady state component
  if fluid_a.m_flow > 0.0 then //mass is flowing into the top, direction is downwards so Tank_A.m_flow is (negative), charging
    Tank_A.m_flow = -1.0*fluid_a.m_flow;
    Tank_A.h_in = inStream(fluid_a.h_outflow);
    fluid_a.h_outflow = Tank_A.h_in;
    fluid_b.h_outflow = Tank_A.h_out;
  else //discharging
    Tank_A.m_flow = -1.0*fluid_a.m_flow;
    Tank_A.h_in = inStream(fluid_b.h_outflow);
    fluid_a.h_outflow = Tank_A.h_out;
    fluid_b.h_outflow = Tank_A.h_in;
  end if;

  fluid_a.p = p_amb;
  fluid_b.p = p_amb;
  T_amb = Tank_A.T_amb;
  T_top_measured = Tank_A.Tf[Nz];
  T_bot_measured = Tank_A.Tf[1];

  p_drop_total = Tank_A.p_drop_total;
  W_dot_loss_pump = Tank_A.W_loss_pump;

annotation(
    Icon(graphics = {
    Text(origin = {7, -16}, extent = {{-49, 42}, {35, -10}}, textString = "TC"), 
    Ellipse(origin = {-36, 56}, fillColor = {240, 50, 0}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}), 
    Ellipse(origin = {-16, 56}, fillColor = {240, 50, 0}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}), 
    Ellipse(origin = {-26, 56}, fillColor = {240, 50, 0}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}), 
    Ellipse(origin = {24, 56}, fillColor = {240, 50, 0}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-6, 56}, fillColor = {240, 50, 0}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {4, 56}, fillColor = {240, 50, 0}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {14, 56}, fillColor = {240, 50, 0}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {34, 56}, fillColor = {240, 50, 0}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-30, 46}, fillColor = {220, 50, 20}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-20, 46}, fillColor = {220, 50, 20}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-10, 46}, fillColor = {220, 50, 20}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {30, 46}, fillColor = {220, 50, 20}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {0, 46}, fillColor = {220, 50, 20}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {10, 46}, fillColor = {220, 50, 20}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {20, 46}, fillColor = {220, 50, 20}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-36, 36}, fillColor = {200, 50, 40}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {34, 36}, fillColor = {200, 50, 40}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {24, 36}, fillColor = {200, 50, 40}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {14, 36}, fillColor = {200, 50, 40}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {4, 36}, fillColor = {200, 50, 40}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-6, 36}, fillColor = {200, 50, 40}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-6, 36}, fillColor = {200, 50, 40}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-16, 36}, fillColor = {200, 50, 40}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-26, 36}, fillColor = {200, 50, 40}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-30, 26}, fillColor = {180, 50, 60}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {30, 26}, fillColor = {180, 50, 60}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {20, 26}, fillColor = {180, 50, 60}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {10, 26}, fillColor = {180, 50, 60}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {10, 26}, fillColor = {180, 50, 60}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {0, 26}, fillColor = {180, 50, 60}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-10, 26}, fillColor = {180, 50, 60}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-20, 26}, fillColor = {180, 50, 60}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-36, 16}, fillColor = {160, 50, 80}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {34, 16}, fillColor = {160, 50, 80}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-30, 6}, fillColor = {140, 50, 100}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {30, 6}, fillColor = {140, 50, 100}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-36, -4}, fillColor = {120, 50, 120}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {34, -4}, fillColor = {120, 50, 120}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-30, -14}, fillColor = {100, 50, 140}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {30, -14}, fillColor = {100, 50, 140}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-36, -24}, fillColor = {80, 50, 160}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {34, -24}, fillColor = {80, 50, 160}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {24, -24}, fillColor = {80, 50, 160}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {14, -24}, fillColor = {80, 50, 160}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {4, -24}, fillColor = {80, 50, 160}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-6, -24}, fillColor = {80, 50, 160}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-16, -24}, fillColor = {80, 50, 160}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-26, -24}, fillColor = {80, 50, 160}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-30, -34}, fillColor = {60, 50, 180}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {30, -34}, fillColor = {60, 50, 180}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {20, -34}, fillColor = {60, 50, 180}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {10, -34}, fillColor = {60, 50, 180}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {0, -34}, fillColor = {60, 50, 180}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-10, -34}, fillColor = {60, 50, 180}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-20, -34}, fillColor = {60, 50, 180}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-36, -44}, fillColor = {40, 50, 200}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {34, -44}, fillColor = {40, 50, 200}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {24, -44}, fillColor = {40, 50, 200}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {14, -44}, fillColor = {40, 50, 200}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {4, -44}, fillColor = {40, 50, 200}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-6, -44}, fillColor = {40, 50, 200}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-16, -44}, fillColor = {40, 50, 200}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-26, -44}, fillColor = {40, 50, 200}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-30, -54}, fillColor = {40, 50, 200}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {30, -54}, fillColor = {40, 50, 200}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {20, -54}, fillColor = {40, 50, 200}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {10, -54}, fillColor = {40, 50, 200}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {0, -54}, fillColor = {40, 50, 200}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-10, -54}, fillColor = {40, 50, 200}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Ellipse(origin = {-20, -54}, fillColor = {40, 50, 200}, fillPattern = FillPattern.Sphere, extent = {{-4, 4}, {6, -6}}),
    Text(origin = {-50, 50}, extent = {{-16, 7}, {6, -3}}, textString = "T_amb"), 
    Text(origin = {-50,-55}, extent = {{-16, 10}, {8, -6}}, textString = "p_amb"), 
    Text(origin = { 75,-128}, lineColor = {0, 0, 255}, extent = {{-150, 150}, {150, 110}}, textString = "%name")}, 
    coordinateSystem(initialScale = 0.1)), 
    Documentation(
		info="<html>
		<p>This model contains the fluid_a (top) and fluid_b (bottom) ports, basically a complete CSP component. 
		This model simply connects the Thermocline_Spheres_Section models to the correct ports.</p>
		</html>",
        revisions ="<html>
        <ul>
		<li><i>Dec 2020</i> by Zebedee Kee:<br>
		Released first version.</li>
		<li><i>Feb 2026</i> by Armando Fontalvo:<br>
		Simplification for verification purposes.</li>
		</ul>
		</html>"));
end SingleTank_Final_Lumped;

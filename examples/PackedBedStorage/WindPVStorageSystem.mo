within PackedBedStorage;
model WindPVStorageSystem
  extends Modelica.Icons.Example;
  import Modelica.SIunits.Conversions.*;
  import Modelica.Constants.*;
  parameter String PV_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Renewable/PV_Latrobe_CF.motab");
  parameter String Wind_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Renewable/Wind_Latrobe_CF.motab");
  parameter String schd_input = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Schedules/Latrobe_Qflow.motab");
  parameter String wea_input = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Weather/Latrobe_38.21S_146.47E_2024.motab");

  replaceable package Medium = SolarTherm.Media.Air.Air_CoolProp_1bar;
  replaceable package Fluid = SolarTherm.Materials.Air_CoolProp_Table_1bar;
  replaceable package Filler = SolarTherm.Materials.Steatite;
  
  //Parameter Inputs
  parameter Real RM = 2.0 "Renewable Multiple (pre-transmission oversizing)";
  parameter Real HM = 2.0 "Heater Multiple";
  parameter Real PV_fraction = 0.5 "PV_fraction";
  parameter Real t_storage = 4 "Hours of storage (hours)";
  parameter Real util_storage_des = 0.2683; //Utilisation determined via component-level analysis
  parameter Real level_storage_mid = 0.5081; //Midpoint of minimum and maximum storage levels determine via component-level analysis
  
  //Heater Parameters
  parameter Real eff_heater = 0.99 "Electrical-to-heat conversion efficiency of the heater";
  
  //Renewable Parameters
  parameter Modelica.SIunits.Power P_renewable_des = RM*P_heater_des;
  parameter Modelica.SIunits.HeatFlowRate Q_heater_des = HM*Q_boiler_des;
  parameter Modelica.SIunits.Power P_heater_des = Q_heater_des/eff_heater;
  parameter Modelica.SIunits.Power PV_ref_size = 1;
  parameter Modelica.SIunits.Power Wind_ref_size = 1;
  //Results
  Modelica.SIunits.Energy E_supplied(start=0) "Energy supplied by the boiler to the industrial process (J)";
  Modelica.SIunits.Energy E_demand(start=0) "Energy demanded by the industrial process (J)";
  Real Capacity_Factor(start=0) "Capacity factor of the system";
  
  //Discretisation and geometry
  parameter Integer Nz = 30;
  
  //Misc Parameters
  parameter Integer Correlation = 1;
  
  //Temperature Controls
  parameter Modelica.SIunits.Temperature T_min = 613 "Minimum system temperature (K)";
  parameter Modelica.SIunits.Temperature T_max = 1173 "Maximum system temperature (K)";
  parameter Modelica.SIunits.Temperature T_boiler_start = T_min + 0.85*(T_max-T_min) "Temperature above-which TES can start discharge (K)";
  parameter Modelica.SIunits.Temperature T_boiler_min = T_boiler_start - 50 "Temperature below-which TES stops discharge (K)";
  parameter Modelica.SIunits.Temperature T_heater_max = T_max - 0.85*(T_max-T_min) "Temperature above-which TES stops charging (K)";
  parameter Modelica.SIunits.Temperature T_heater_start = T_heater_max - 50 "Temperature below-which TES can start charging (K)";

  //Level-Controls
  parameter Modelica.SIunits.Time t_stor_start_dis = 0.1*t_storage "Number of effective storage seconds stored before TES can start discharging (1 hour)";  
  
  //Calculated Parameters
  parameter Modelica.SIunits.Energy E_max = t_storage * 3600.0 * Q_boiler_des "Maximum tank stored energy (J)";
  parameter Modelica.SIunits.HeatFlowRate Q_boiler_des = 5e6 "Heat-rate to boiler at design (W)";
  parameter Modelica.SIunits.MassFlowRate m_boiler_des = Q_boiler_des/(h_air_hot_set-h_air_cold_set) "Design boiler input mass flow rate";

  parameter Modelica.SIunits.Temperature T_hot_set = T_max "Ideal hot temperature of the system";
  parameter Modelica.SIunits.Temperature T_cold_set = T_min "Ideal cold temperature of the system";
  parameter Medium.ThermodynamicState state_air_cold_set = Medium.setState_pTX(Medium.p_default, T_cold_set) "Cold air thermodynamic state at design";
  parameter Medium.ThermodynamicState state_air_hot_set = Medium.setState_pTX(Medium.p_default, T_hot_set) "Hold air thermodynamic state at design";
  parameter Modelica.SIunits.SpecificEnthalpy h_air_cold_set = Medium.specificEnthalpy(state_air_cold_set) "Cold air specific enthalpy at design";
  parameter Modelica.SIunits.SpecificEnthalpy h_air_hot_set = Medium.specificEnthalpy(state_air_hot_set) "Hot air specific enthalpy at design";
  
  parameter Modelica.SIunits.Diameter D_tank = 3.44;
  parameter Modelica.SIunits.Height H_tank = 8;
  parameter Real epsilon = 0.4;
  parameter Modelica.SIunits.CoefficientOfHeatTransfer U_wall = 0.339 "W/m2K";
  parameter Modelica.SIunits.Diameter ds = 0.02 "Filler sphere diameter";
  parameter SI.Temperature T_start = T_cold_set "Initial (uniform) temperature of all components (K), defaults to T_min";

  SolarTherm.Models.Storage.Thermocline.Spheres.SingleTank_Final_Lumped TES(
    redeclare package Medium = Medium, 
    redeclare package Fluid_Package = Fluid, 
    redeclare package Filler_Package = Filler, 
    D_tank = D_tank,
    H_tank = H_tank,
    epsilon = epsilon,
    ds = ds,
    Nz = Nz, 
    T_max = T_hot_set, 
    T_min = T_cold_set, 
    T_start = T_start,
    Correlation = Correlation,
    U_wall = U_wall,
    E_max = E_max) 
    annotation(Placement(visible = true, transformation(origin = {32, 2}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));

  SolarTherm.Models.Fluid.Pumps.PumpSimple_EqualPressure pumpCold(redeclare package Medium = Medium) 
    annotation(Placement(visible = true, transformation(origin = {-18, -78}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));

  SolarTherm.Models.Fluid.Valves.PBS_TeeJunction_LoopBreaker Splitter_Top(redeclare package Medium = Medium) 
    annotation(Placement(visible = true, transformation(origin = {32, 67.6249}, extent = {{-18, -13.9366}, {18, 13.9366}}, rotation = 0)));

  SolarTherm.Models.Fluid.Valves.PBS_TeeJunction Splitter_Bot(redeclare package Medium = Medium) 
    annotation(Placement(visible = true, transformation(origin = {32, -58}, extent = {{17, 0}, {-17, -22.039}}, rotation = 0)));

  SolarTherm.Models.Fluid.Pumps.PumpSimple_EqualPressure pumpHot(redeclare package Medium = Medium) 
    annotation(Placement(visible = true, transformation(origin = {91, 79}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));

  Modelica.Blocks.Sources.RealExpression Tamb(y = weather.y[1]) 
    annotation(Placement(visible = true, transformation(origin = {-5, 17}, extent = {{-9, -12}, {9, 12}}, rotation = 0)));

  Modelica.Blocks.Sources.RealExpression p_amb(y = 101325) 
    annotation(Placement(visible = true, transformation(origin = {-5, -13}, extent = {{-9, -12}, {9, 12}}, rotation = 0)));
  
  SolarTherm.Models.Control.WindPV_Thermocline_Control Control(
    redeclare package HTF = Medium, 
    E_max = E_max, 
    Q_flow_boiler_des = Q_boiler_des, 
    T_boiler_min = T_boiler_min, 
    T_boiler_start = T_boiler_start, 
    T_heater_max = T_heater_max, 
    T_heater_start = T_heater_start, 
    T_target = T_max, 
    util_storage_des = util_storage_des, 
    h_target = h_air_hot_set, 
    level_mid = level_storage_mid, 
    m_flow_0 = 1e-8, 
    m_flow_boiler_des = m_boiler_des, 
    m_flow_min = 1e-8, 
    m_flow_tol = 0.001 * m_boiler_des, 
    t_stor_start_dis = t_stor_start_dis, 
    t_wait = 1.0 * 3600.0)
      annotation(Placement(visible = true,transformation(origin = {114, 16},extent = {{-10, -10}, {10, 10}},rotation = 0)));

  Modelica.Blocks.Sources.CombiTimeTable weather(
    fileName = wea_input,
    columns = {2,3},
    smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
    tableName = "weather",
    tableOnFile = true);
      //annotation(Placement(visible = true,transformation(origin = {-100, 75},extent = {{-10, -10}, {10, 10}},rotation = 0)));

  Modelica.Blocks.Sources.CombiTimeTable Q_schd(
    fileName = schd_input, 
    smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative, 
    tableName = "Q_flow", 
    tableOnFile = true) 
      annotation(Placement(visible = true,transformation(origin = {140, 48},extent = {{10, -10}, {-10, 10}},rotation = 0)));

  SolarTherm.Models.Fluid.HeatExchangers.Boiler_Basic Boiler(
    redeclare package Medium = Medium,
    T_cold_set = T_cold_set,
    T_hot_set = T_hot_set)
      annotation(Placement(visible = true,transformation(origin = {158, 0},extent = {{-10, -10}, {10, 10}},rotation = 0)));

  SolarTherm.Models.CSP.CRS.Receivers.Basic_Heater heater(
    redeclare package Medium = Medium,
    P_heater_des = P_heater_des,
    Q_flow_heater_des = Q_heater_des, 
    eff_heater = eff_heater, 
    T_cold_set = T_cold_set, 
    T_hot_set = T_hot_set) 
      annotation(Placement(visible = true,transformation(origin = {-46, 10},extent = {{-10, -10}, {10, 10}},rotation = 90)));

  Modelica.Blocks.Sources.CombiTimeTable PV_input(
    fileName = PV_file, 
    tableName = "Power", 
    tableOnFile = true, 
    smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative) 
      annotation(Placement(visible = true,transformation(origin = {-124, 34},extent = {{-10, -10}, {10, 10}},rotation = 0)));

  Modelica.Blocks.Math.Add Grid_Sum(
    k1 =  P_renewable_des *PV_fraction / PV_ref_size, 
    k2 =  P_renewable_des *(1.0 - PV_fraction) / Wind_ref_size) 
      annotation(Placement(visible = true,transformation(origin = {-85, 10},extent = {{-10, -10}, {10, 10}},rotation = 0)));

  Modelica.Blocks.Sources.CombiTimeTable Wind_input(
    fileName = Wind_file,
    smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
    tableName = "Power",
    tableOnFile = true)
      annotation(Placement(visible = true,transformation(origin = {-124, 4},extent = {{-10, -10}, {10, 10}},rotation = 0)));

equation
  der(E_supplied) = Boiler.Q_flow;
  der(E_demand) = Control.Q_flow_demand;
  if time > 86400.0 then
    Capacity_Factor = E_supplied/E_demand;
  else
    Capacity_Factor = 0.0;
  end if;

  connect(TES.T_top_measured, Control.T_top_tank);
  connect(TES.T_bot_measured, Control.T_bot_tank);
  connect(TES.h_top_outlet, Control.h_tank_top);
  connect(TES.h_bot_outlet, Control.h_tank_bot);
  connect(TES.fluid_b, Splitter_Bot.fluid_c) annotation(
    Line(points = {{32, -28}, {32, -68}}, color = {0, 127, 255}, thickness = 0.5));
  connect(Tamb.y, TES.T_amb) annotation(
    Line(points = {{5, 17}, {13, 17}}, color = {0, 0, 127}));
  connect(p_amb.y, TES.p_amb) annotation(
    Line(points = {{5, -13}, {13, -13}}, color = {0, 0, 127}));
  connect(Splitter_Bot.fluid_b, pumpCold.fluid_a) annotation(
    Line(points = {{17, -77}, {-8, -77}}, color = {0, 127, 255}, thickness = 0.5));
  connect(TES.Level, Control.Level) annotation(
    Line(points = {{49, 16}, {103, 16}}, color = {0, 0, 127}));
  connect(Splitter_Top.fluid_c, TES.fluid_a) annotation(
    Line(points = {{32, 67}, {32, 32}}, color = {0, 127, 255}, thickness = 0.5));
  connect(Splitter_Top.fluid_b, pumpHot.fluid_a) annotation(
    Line(points = {{47, 79}, {82, 79}}, color = {0, 127, 255}, thickness = 0.5));
  connect(Boiler.fluid_b, Splitter_Bot.fluid_a) annotation(
    Line(points = {{158, -10}, {158, -77}, {45, -77}}, color = {0, 127, 255}, thickness = 0.5));
  connect(pumpHot.fluid_b, Boiler.fluid_a) annotation(
    Line(points = {{100, 79}, {158, 79}, {158, 10}}, color = {0, 127, 255}, thickness = 0.5));
  connect(pumpCold.fluid_b, heater.fluid_a) annotation(
    Line(points = {{-28, -78}, {-46, -78}, {-46, 1}}, color = {0, 127, 255}, thickness = 0.5));
  connect(heater.fluid_b, Splitter_Top.fluid_a) annotation(
    Line(points = {{-46, 19}, {-46, 79}, {18, 79}}, color = {0, 127, 255}, thickness = 0.5));
  connect(Control.curtail, heater.curtail) annotation(
    Line(points = {{114, 5}, {114, -45}, {-55.5, -45}, {-55.5, -2}}, color = {255, 0, 255}, pattern = LinePattern.Dash));
  connect(Control.Q_flow_curtail, heater.Q_flow_curtail) annotation(
    Line(points = {{108, 5}, {108, -35}, {-50, -35}, {-50, -2}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
  connect(Q_schd.y[1], Control.Q_flow_demand) annotation(
    Line(points = {{129, 48}, {120, 48}, {120, 27}}, color = {0, 0, 127}));
  connect(PV_input.y[1], Grid_Sum.u1) annotation(
    Line(points = {{-113, 34}, {-104, 34}, {-104, 16}, {-96, 16}}, color = {0, 0, 127}));
  connect(Wind_input.y[1], Grid_Sum.u2) annotation(
    Line(points = {{-113, 4}, {-96, 4}}, color = {0, 0, 127}));
  connect(Grid_Sum.y, heater.P_supply) annotation(
    Line(points = {{-80, 10}, {-57, 10}}, color = {0, 0, 127}));
  connect(Boiler.h_out_signal, Control.h_boiler_outlet) annotation(
    Line(points = {{149, 7}, {149, 16}, {124, 16}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
  connect(heater.Q_flow_heater_raw, Control.Q_flow_heater_raw) annotation(
    Line(points = {{-53, 21}, {-53, 48}, {108, 48}, {108, 27}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
  connect(Control.m_flow_heater_signal, pumpCold.m_flow) annotation(
    Line(points = {{120, 5}, {120, -58}, {-18, -58}, {-18, -70}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
  connect(Control.m_flow_boiler_signal, pumpHot.m_flow) annotation(
    Line(points = {{114, 27}, {114, 92}, {91, 92}, {91, 86}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
  annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-200, -100}, {200, 100}}, initialScale = 0.1)),
    Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false)), 
    experiment(StopTime = 3.1536e+07, StartTime = 0, Tolerance = 1.0e-6, Interval = 300, maxStepSize = 60, initialStepSize = 60));

end WindPVStorageSystem;

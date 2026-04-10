within SolarTherm.Systems;
model PBS_1T_Spheres_Verification "Packed-bed storage with air and spheres os steatite"
    import SI = Modelica.SIunits;
    import CN = Modelica.Constants;
    import CV = Modelica.SIunits.Conversions;
    extends Modelica.Icons.Example;
    package Medium = SolarTherm.Media.Air.Air_CoolProp_1bar;
    package Fluid_Package = SolarTherm.Materials.Air_CoolProp_Table_1bar;
    package Filler_Package = SolarTherm.Materials.Steatite;

    // Case selection
    parameter Verification.ScenarioBank data;
    parameter Verification.CaseData currentCase = data.case1;

    //Heat Transfer Convection Coefficient
    parameter Integer Correlation = 1 "Wakao & Kaguei";

    //Numerical Discretisation
    parameter Integer Nz = 200 "Number of fluid CVs in each tank";//360

    //Design Parameters
    parameter SI.Energy E_max = 20 * 3.6e9 "Storage capacity (J)";
    parameter SI.Temperature T_min = 613 "Minimum temperature";
    parameter SI.Temperature T_max = 1173 "Maximum temperature";
    parameter SI.Temperature T_start = 293 "Packed-bed initial temperature";
    final parameter SI.Temperature T_stop_charging = T_max - 0.85*(T_max - T_min);
    final parameter SI.Temperature T_stop_discharging = T_min + 0.85*(T_max - T_min);
    parameter Real epsilon = 0.4 "Packed-bed porosity";
    parameter SI.Length ds = 0.03 "Filler diameter";
    parameter SI.Length H_tank = 8 "Tank height";
    parameter SI.Diameter D_tank = 3.44 "Tank diameter";
    parameter SI.CoefficientOfHeatTransfer U_wall = 0.678 "W/m2K";

    // Calculated parameters
    parameter SI.SpecificEnthalpy h_f_min = Fluid_Package.h_Tf(T_min, 0);
    parameter SI.SpecificEnthalpy h_f_max = Fluid_Package.h_Tf(T_max, 1);
    parameter SI.Power P_charging = 10e6 "Charging rate (W)";
    parameter SI.Power P_discharging = 5e6 "Discharging rate (W)";
    parameter SI.MassFlowRate m_flow_charge = P_charging / (h_f_max - h_f_min);
    parameter SI.MassFlowRate m_flow_discharge = P_discharging / (h_f_max - h_f_min);

    //Models
    SolarTherm.Models.Storage.Thermocline.Spheres.SingleTank_Final_Lumped TES(
        redeclare package Medium = Medium,
        redeclare package Fluid_Package = Fluid_Package,
        redeclare package Filler_Package = Filler_Package,
        Nz = Nz,
        T_max = T_max,
        T_min = T_min,
        T_start = T_start,
        E_max = E_max,
        H_tank = H_tank,
        D_tank = D_tank,
        epsilon = epsilon,
        ds = ds,
        U_wall = U_wall,
        Correlation = Correlation) 
        annotation(Placement(
            visible = true, 
            transformation(
                origin = {0, 0}, 
                extent = {{-20, -20}, {20, 20}}, 
                rotation = 0)));
    
    Modelica.Fluid.Sources.FixedBoundary Source(
        redeclare package Medium = Medium, 
        T = T_max,
        p = 101325,
        nPorts = 1)
        annotation(Placement(
            visible = true, 
            transformation(
                origin = {-85, 50}, 
                extent = {{-10, -10}, {10, 10}}, 
                rotation = 0)));

    Modelica.Fluid.Sources.FixedBoundary Source2(
        redeclare package Medium = Medium, 
        T = T_min,
        nPorts = 1,
        p = 101325)
        annotation(Placement(
            visible = true, 
            transformation(
                origin = {85, -50}, 
                extent = {{10, -10}, {-10, 10}}, 
                rotation = 0)));

    SolarTherm.Models.Fluid.Sources.FluidSink2 Sink(
        redeclare package Medium = Medium)
        annotation(Placement(
            visible = true, 
            transformation(
                origin = {-85, -50}, 
                extent = {{10, -10}, {-10, 10}}, 
                rotation = 0)));

    SolarTherm.Models.Fluid.Sources.FluidSink2 Sink2(
        redeclare package Medium = Medium)
        annotation(Placement(
            visible = true, 
            transformation(
                origin = {85, 50}, 
                extent = {{-10, -10}, {10, 10}}, 
                rotation = 0)));

    Modelica.Blocks.Sources.RealExpression m_flow_chg(y = m_char) 
        annotation(Placement(
            visible = true, 
            transformation(
                origin = {-80, 80},
                extent = {{-10, -10}, {10, 10}},
                rotation = 0)));

    Modelica.Blocks.Sources.RealExpression m_flow_dis(y = m_disc) 
        annotation(Placement(
            visible = true, 
            transformation(
                origin = {80, 80},
                extent = {{10, -10}, {-10, 10}},
                rotation = 0)));

    SolarTherm.Models.Fluid.Pumps.PumpSimple pump(
        redeclare package Medium = Medium)
        annotation(Placement(
            visible = true,
            transformation(
                origin = {-30, 50},
                extent = {{-10, -10}, {10, 10}},
                rotation = 0)));

    SolarTherm.Models.Fluid.Pumps.PumpSimple_EqualPressure pump2(
        redeclare package Medium = Medium)
        annotation(Placement(
            visible = true, 
            transformation(
                origin = {30, 50}, 
                extent = {{-10, -10}, {10, 10}}, 
                rotation = 0)));

    SolarTherm.Models.Fluid.Pumps.PumpSimple_EqualPressure pump3(
        redeclare package Medium = Medium)
        annotation(Placement(
            visible = true, 
            transformation(
                origin = {-50, -50}, 
                extent = {{10, -10}, {-10, 10}}, 
                rotation = 0)));

    SolarTherm.Models.Fluid.Pumps.PumpSimple pump4(
        redeclare package Medium = Medium)
        annotation(Placement(
            visible = true, 
            transformation(
                origin = {50, -50}, 
                extent = {{10, -10}, {-10, 10}}, 
                rotation = 0)));

    //Ambient conditions
    Modelica.Blocks.Sources.RealExpression T_amb(y = 298.15) 
        annotation(Placement(
            visible = true, 
            transformation(
                origin = {-30, 8},
                extent = {{-10, -10}, {10, 10}},
                rotation = 0)));
    Modelica.Blocks.Sources.RealExpression P_amb(y = 101325) 
        annotation(Placement(
            visible = true,
            transformation(
                origin = {-30,-8},
                extent = {{-10, -10}, {10, 10}},
                rotation = 0)));

    SolarTherm.Models.Fluid.HeatExchangers.mass_loop_breaker breaker(
        redeclare package Medium = Medium) 
        annotation(Placement(
            visible = true, 
            transformation(
                origin = {0, 25}, 
                extent = {{-10, -10}, {10, 10}}, 
                rotation = -90)));

    SolarTherm.Models.Fluid.Valves.PBS_TeeJunction split1(
        redeclare package Medium = Medium)
        annotation(Placement(
            visible = true,
            transformation(
                origin = {0, 28.5},
                extent = {{-10, 0}, {10, 25}},
                rotation = 0)));

    SolarTherm.Models.Fluid.Valves.PBS_TeeJunction split2(
        redeclare package Medium = Medium)
        annotation(Placement(
            visible = true,
            transformation(
                origin = {0, -28.5},
                extent = {{-10, 0}, {10, 25}},
                rotation = 180)));

    //Mass flow Signals starts in charging state
    SI.MassFlowRate m_char(start = m_flow_charge);
    SI.MassFlowRate m_disc(start = 0);

    // Control
    parameter Real h_standby = 12;
    final parameter SI.Time t_standby = 3600*h_standby;

    // Control Variables
    Modelica.Blocks.Continuous.LimPID pid_chg(
      Ti = 60,
      k = 1,
      yMin = 0,
      yMax = m_flow_charge,
      y_start = m_flow_charge,
      initType = Modelica.Blocks.Types.InitPID.InitialOutput,
      limitsAtInit = true);

    Integer state;
    SI.Time t_next_event;

algorithm
    when TES.Tank_A.Tf[1] > T_stop_charging then
        t_next_event := time + t_standby;
        state := 1;
    end when;
    when time > t_next_event and state < 2 then
        t_next_event := time + t_standby;
        state := 2;
    end when;
    when TES.Tank_A.Tf[Nz] < T_stop_discharging and state > 0 then
        t_next_event := time + t_standby;
        state := 0;
    end when;

equation
    // Controlled
    pid_chg.u_m = TES.Tank_A.Tf[1];
    pid_chg.u_s = T_stop_charging;

    // State Logic
    if state == 2 then
        m_char = 0;
        m_disc = m_flow_discharge;
    elseif state == 1 then
        m_char = 1e-12;
        m_disc = 0;
    else
        m_char = pid_chg.y;
        m_disc = 0;
    end if;

    connect(T_amb.y, TES.T_amb) annotation(
    Line(points = {{-19, 8}, {-10, 8}}, color = {0, 0, 127}));
    connect(P_amb.y, TES.p_amb) annotation(
    Line(points = {{-19, -8}, {-10, -8}}, color = {0, 0, 127}));
    connect(Source.ports[1], pump.fluid_a) annotation(
    Line(points = {{-74, 50}, {-40, 50}}, color = {0, 127, 255}));
    connect(pump.fluid_b, split1.fluid_a) annotation(
    Line(points = {{-20, 50}, {-8, 50}}));
    connect(split1.fluid_b, pump2.fluid_a) annotation(
    Line(points = {{8, 50}, {20, 50}}, color = {0, 127, 255}));
    connect(split1.fluid_c, breaker.port_a) annotation(
    Line(points = {{0, 41}, {0, 32}}, color = {0, 127, 255}));
    connect(breaker.port_b, TES.fluid_a) annotation(
    Line(points = {{0, 20}, {0, 16}}, color = {0, 127, 255}));
    connect(pump2.fluid_b, Sink2.port_a) annotation(
    Line(points = {{40, 50}, {76, 50}}, color = {0, 127, 255}));
    connect(pump3.fluid_b, Sink.port_a) annotation(
    Line(points = {{-60, -50}, {-74, -50}}, color = {0, 127, 255}));
    connect(pump4.fluid_a, Source2.ports[1]) annotation(
    Line(points = {{60, -50}, {76, -50}}, color = {0, 127, 255}));
    connect(m_flow_chg.y, pump.m_flow) annotation(
    Line(points = {{-68, 80}, {-30, 80}, {-30, 59}}, color = {0, 0, 127}));
    connect(m_flow_chg.y, pump3.m_flow) annotation(
    Line(points = {{-68, 80}, {-50, 80}, {-50, -42}}, color = {0, 0, 127}));
    connect(split2.fluid_c, TES.fluid_b) annotation(
    Line(points = {{0, -40}, {0, -16}}, color = {0, 127, 255}));
    connect(pump4.fluid_b, split2.fluid_a) annotation(
    Line(points = {{40, -50}, {8, -50}}, color = {0, 127, 255}));
    connect(split2.fluid_b, pump3.fluid_a) annotation(
    Line(points = {{-8, -50}, {-40, -50}}));
    connect(m_flow_dis.y, pump2.m_flow) annotation(
    Line(points = {{70, 80}, {30, 80}, {30, 59}}, color = {0, 0, 127}));
    connect(m_flow_dis.y, pump4.m_flow) annotation(
    Line(points = {{70, 80}, {50, 80}, {50, -42}}, color = {0, 0, 127}));

annotation(
    experiment(StopTime = 864000, StartTime = 0, Tolerance = 1e-6, Interval = 60),
    Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false)),
    Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false)),
    Documentation(info =
        "<html>
        <ul>
        <li> <i>Dec 2020</i> by Z. Kee:<br> Resealed first version. </li>
        <li> <i>Feb 2026</i> by A. Fontalvo:<br> Simplification for verification purposes. </li>
        </ul>
        </html>"));
end PBS_1T_Spheres_Verification;

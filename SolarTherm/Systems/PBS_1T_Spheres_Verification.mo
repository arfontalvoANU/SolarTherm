within SolarTherm.Systems;

model PBS_1T_Spheres_Verification "Packed-bed storage with air and spheres os steatite"
    import SI = Modelica.SIunits;
    import CN = Modelica.Constants;
    import CV = Modelica.SIunits.Conversions;
    extends Modelica.Icons.Example;
    package Medium = SolarTherm.Media.Air.Air_CoolProp_1bar;
    package Fluid_Package = SolarTherm.Materials.Air_CoolProp_Table_1bar;
    package Filler_Package = SolarTherm.Materials.Steatite;

    //Heat Transfer Convection Coefficient
    parameter Integer Correlation = 1 "Wakao & Kaguei";

    //Numerical Discretisation
    parameter Integer N_f = 100 "Number of fluid CVs in each tank";//360
    parameter Integer N_p = 10 "Number of filler CVs  in main tank";

    //Design Parameters
    parameter SI.Energy E_max = 2.2691652e7 "Storage capacity (J)";
    parameter SI.Temperature T_max = 823 "Maximum temperature";
    parameter SI.Temperature T_min = 293 "Minimum temperature";
    parameter Real eta = 0.4 "Packed-bed porosity"; 
    parameter Real ar = 1.2/0.148 "Tank aspect ratio";
    parameter SI.CoefficientOfHeatTransfer U_loss_tank = 0.678 "W/m2K";
    parameter SI.Length d_p = 0.02 "Filler diameter";
    parameter SI.MassFlowRate m_flow = 0.255 * 0.25*Modelica.Constants.pi*TES.Tank_A.D_tank^2;

    //Models
    SolarTherm.Models.Storage.Thermocline.Spheres.SingleTank_Final_Lumped TES(
        redeclare package Medium = Medium,
        redeclare package Fluid_Package = Fluid_Package,
        redeclare package Filler_Package = Filler_Package,
        N_f = N_f,
        N_p = N_p,
        T_max = T_max,
        T_min = T_min,
        E_max = E_max,
        ar = ar,
        eta = eta,
        d_p = d_p,
        U_loss_tank = U_loss_tank,
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

    Modelica.Blocks.Sources.RealExpression m_flow_chg(y = m_flow) 
        annotation(Placement(
            visible = true, 
            transformation(
                origin = {-80, 80},
                extent = {{-10, -10}, {10, 10}},
                rotation = 0)));

    SolarTherm.Models.Fluid.Pumps.PumpSimple pump(
        redeclare package Medium = Medium)
        annotation(Placement(
            visible = true,
            transformation(
                origin = {-50, 50},
                extent = {{-10, -10}, {10, 10}},
                rotation = 0)));

    SolarTherm.Models.Fluid.Sources.FluidSink2 Sink(
        redeclare package Medium = Medium)
        annotation(Placement(
            visible = true, 
            transformation(
                origin = {-50, -50}, 
                extent = {{10, -10}, {-10, 10}}, 
                rotation = 0)));

    //Ambient conditions
    Modelica.Blocks.Sources.RealExpression T_amb(y = 298.15) 
        annotation(Placement(
            visible = true, 
            transformation(
                origin = {-50, 8},
                extent = {{-10, -10}, {10, 10}},
                rotation = 0)));
    Modelica.Blocks.Sources.RealExpression P_amb(y = 101325) 
        annotation(Placement(
            visible = true,
            transformation(
                origin = {-50,-8},
                extent = {{-10, -10}, {10, 10}},
                rotation = 0)));

equation

    connect(T_amb.y, TES.T_amb)
        annotation(
            Line(
                points = {{-38, 8}, {-10, 8}},
                color = {0, 0, 127}));
    connect(P_amb.y, TES.p_amb)
        annotation(
            Line(
                points = {{-38, -8}, {-10, -8}},
                color = {0, 0, 127}));
    connect(TES.fluid_b, Sink.port_a) 
        annotation(
            Line(
                points = {{0, -16}, {0, -50}, {-40, -50}},
                color = {0, 127, 255}));
    connect(Source.ports[1], pump.fluid_a) 
        annotation(
            Line(
                points = {{-74, 50}, {-60, 50}},
                color = {0, 127, 255}));
    connect(m_flow_chg.y, pump.m_flow)
        annotation(
            Line(
                points = {{-68, 80}, {-50, 80}, {-50, 58}},
                color = {0, 0, 127}));
    connect(pump.fluid_b, TES.fluid_a)
        annotation(
            Line(
                points = {{-40, 50}, {0, 50}, {0, 16}},
                color = {0, 127, 255}));

annotation(
    experiment(StopTime = 4800, StartTime = 0, Tolerance = 1e-4, Interval = 60),
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
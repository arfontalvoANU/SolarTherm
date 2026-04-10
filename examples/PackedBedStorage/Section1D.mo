within PackedBedStorage;
model Section1D "Heat transfer model of thermocline tank with spherical fillers"
  extends Modelica.Icons.Example;
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  import Tables = Modelica.Blocks.Tables;

  //Initialize Material Packages
  replaceable package Medium = SolarTherm.Materials.Air_CoolProp_Table_1bar constrainedby SolarTherm.Materials.PartialMaterial "Fluid Package";

  // Filler options
  parameter SI.SpecificHeatCapacity cps = 1068.0 "Filler heat capacity (J/kg/K)";
  parameter SI.Density rhos = 2680.0 "Filler density (kg/m3)";
  parameter SI.ThermalConductivity ks = 2.5 "Filler thermal conductivity (W/m/K)";

  //Tank Design parameters
  parameter Real epsilon = 0.4 "Porosity";
  parameter Real ds = 0.03 "Diameter of sphere (particle) (m)";

  //Temperature Bounds
  parameter SI.Temperature T_min = 613 "Design cold Temperature of everything in the tank (K)";
  parameter SI.Temperature T_max = 1173 "Design hot Temperature of everything in the tank (K)";
  parameter SI.Temperature T_start = 293 "Initial (uniform) temperature of all components (K), defaults to T_min";
  final parameter SI.Temperature T_stop_charging = T_max - 0.85*(T_max - T_min);
  final parameter SI.Temperature T_stop_discharging = T_min + 0.85*(T_max - T_min);
  parameter SI.Power P_charging = 10e6;
  parameter SI.Power P_discharging = 5e6;
  parameter SI.MassFlowRate m_flow_chg = P_charging/(hf_max - hf_min);
  parameter SI.MassFlowRate m_flow_dis = P_discharging/(hf_max - hf_min);

  //Calculated Tank Design Parameters
  parameter SI.Length H_tank = 8.0;
  parameter SI.Diameter D_tank = 3.44;
  parameter Real f_area = 0.95;
  parameter SI.Area A = 0.25 * CN.pi * D_tank^2 "Cross sectional area of tank";
  parameter SI.Volume V = A*H_tank;
  parameter SI.Mass ms = (1-epsilon)*V*rhos;
  parameter SI.Density rhof_avg = 0.5*(rhof_max + rhof_min);
  parameter SI.Mass mf = epsilon*V*rhof_avg;
  parameter SI.Energy E_max = mf*(hf_max-hf_min) + ms*cps*(T_max-T_min);
  parameter SI.Energy E_init = mf*(hf_min - h_start) + ms*cps*(T_min - T_start);

  //Thermal Losses
  parameter SI.Area A_loss_tank = CN.pi*D_tank*D_tank*0.5 + CN.pi*D_tank*H_tank "Heat loss area (m2)";
  parameter SI.CoefficientOfHeatTransfer U_wall = 0.678 "Cylinder wall heat loss coeff (W/m2K)";
  parameter SI.CoefficientOfHeatTransfer U_top = 0.0 "Top circle heat loss coeff (W/m2K)";
  parameter SI.CoefficientOfHeatTransfer U_bot = 0.0 "Bottom circle heat loss coeff (W/m2K)";

  //Inititalize temperature and enthalpy profile
  parameter SI.Temperature Tf_start[Nz] = fill(T_start,Nz);
  parameter SI.Temperature Ts_start[Nz] = fill(T_start, Nz);
  parameter SI.SpecificEnthalpy h_start = Medium.h_Tf(T_start, 0.0);
  parameter SI.SpecificEnthalpy hf_start[Nz] = fill(h_start, Nz) "Defaults to uniform";
  //Property bounds
    //Fluid
  parameter SI.SpecificEnthalpy hf_min = Medium.h_Tf(T_min,0) "Starting enthalpy of the HTF";
  parameter SI.SpecificEnthalpy hf_max = Medium.h_Tf(T_max,0) "Starting enthalpy of the HTF";
  parameter SI.Density rhof_min = Medium.rho_Tf(T_min,0);
  parameter SI.Density rhof_max = Medium.rho_Tf(T_max,0);

    //Discretization
  parameter Integer Nz = 200 "Number of finite volume elements in fluid";
  parameter SI.Length dz = H_tank / Nz "discretization vertical length of fluid";
  parameter SI.Length z[Nz] = linspace(0, H_tank, Nz);

  // Control
  parameter Real h_standby = 12;
  final parameter SI.Time t_standby = 3600*h_standby;

  //Initialise Fluid Array
  SI.Temperature Tf[Nz] "(K)";
  SI.SpecificEnthalpy hf[Nz](start = hf_start) "J/kg";

  //Inlet and outlet enthalpies and temperatures
  SI.SpecificEnthalpy h_in(start = hf_max) "Inlet Enthalpy depends on mass flow direction (J/kg)";
  SI.Temperature T_top(start = T_start);
  SI.Temperature T_bot(start = T_start);

  //Mass flow rates and superficial velocity
  SI.MassFlowRate m_flow "kg/s";
  SI.Velocity uf[Nz] "Superficial fluid velocity in packed bed (m/s)";

  //Initialise Particle
  SI.Temperature Ts[Nz](start = Ts_start) "Temperature of particle elements";

  //Calculated Pumping Losses
  SI.Pressure p_drop_total "Sum of all pressure drops";

  //Convection Properties
  Real Re[Nz] "Reynolds";
  Real Pr[Nz] "Prandtl";
  Real Nu[Nz] "Nusselt";
  Real hv[Nz] "Volumetric heat transfer coeff (W/m3K)";

  //Pressure Drop
  SI.Pressure p_drop[Nz] "Pressure drop across each mesh element";

  //Fluid Properties
  SI.ThermalConductivity kf[Nz] "W/mK";
  SI.DynamicViscosity muf[Nz] "Pa.s";
  SI.SpecificHeatCapacity cpf[Nz] "J/kgK";
  SI.Density rhof[Nz] "kg/m3";

  // Importing weather data file
  parameter String weather_file = Modelica.Utilities.Files.loadResource("resources/BARRA-output-local--38.21-146.47-2024.motab");

  // Weather Input
  Modelica.Blocks.Sources.CombiTimeTable weather(
      columns = {2,3},
      fileName = weather_file,
      tableName = "weather",
      tableOnFile = true,
  smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative);

  SI.Temperature T_amb;

  // Control Variables
  Modelica.Blocks.Continuous.LimPID pid_chg(
    Ti = 60,
    k = 1.0,
    yMin = 0,
    yMax = m_flow_chg,
    y_start = m_flow_chg,
    initType = Modelica.Blocks.Types.InitPID.InitialOutput,
    limitsAtInit = true);
  Integer state;
  SI.Time t_next_event;

  // Utilisation
  SI.Energy Ei[Nz];
  SI.Energy E;

protected
  Medium.State fluid[Nz]"Fluid object array";

algorithm
    when Tf[1] > T_stop_charging then
        t_next_event := time + t_standby;
        state := 1;
    end when;
    when time > t_next_event and state < 2 then
        t_next_event := time + t_standby;
        state := 2;
    end when;
    when Tf[Nz] < T_stop_discharging and state > 0 then
        t_next_event := time + t_standby;
        state := 0;
    end when;

initial equation
  for i in 1:Nz loop
    fluid[i].h = hf_start[i];
    Ei[i] = 0;
  end for;
  m_flow = -m_flow_chg;
  t_next_event = 1e6;
  state = 0;

equation
    // Time-dependent ambient temperature
    T_amb = Modelica.SIunits.Conversions.from_degC(weather.y[1]);

    // Controlled
    pid_chg.u_m = Tf[1];
    pid_chg.u_s = T_stop_charging;

    // State Logic
    if state == 2 then
        m_flow = m_flow_dis;
        h_in = hf_min;
    else
        m_flow = -pid_chg.y;
        h_in = hf_max;
    end if;

    // Top and Bottom temperatures
    T_top = Tf[Nz];
    T_bot = Tf[1];

    if m_flow < 0 then
        //Bottom
        rhof[1] * der(hf[1]) =
        -2*kf[1]*kf[2]/(kf[1]+kf[2]) * (Tf[1]-Tf[2])/(dz^2)
        +rhof[1]*uf[1]/epsilon*(hf[1]-hf[2])/dz
        -hv[1]*(Tf[1] - Ts[1])/epsilon
        -U_bot*(Tf[1]-T_amb)/(epsilon*dz)
        -U_wall*CN.pi*D_tank*(Tf[1]-T_amb)/(epsilon*A);

        //Middle
        for i in 2:Nz - 1 loop
            rhof[i] * der(hf[i]) =
            -2*kf[i]*kf[i-1]/(kf[i]+kf[i-1]) * (Tf[i]-Tf[i-1])/(dz^2)
            -2*kf[i]*kf[i+1]/(kf[i]+kf[i+1]) * (Tf[i]-Tf[i+1])/(dz^2)
            +rhof[i]*uf[i]/epsilon*(hf[i]-hf[i+1])/dz
            -hv[i]*(Tf[i]-Ts[i])/epsilon
            -U_wall*CN.pi*D_tank*(Tf[i]-T_amb)/(epsilon*A);
        end for;

        //Top
        rhof[Nz] * der(hf[Nz]) =
        -2*kf[Nz-1]*kf[Nz]/(kf[Nz-1]+kf[Nz]) * (Tf[Nz]-Tf[Nz-1])/(dz^2)
        +rhof[Nz]*uf[Nz]/epsilon*(hf[Nz]-h_in)/dz
        -hv[Nz]*(Tf[Nz]-Ts[Nz])/epsilon
        -U_wall*CN.pi*D_tank*(Tf[Nz]-T_amb)/(epsilon*A)
        -U_top*(Tf[Nz]-T_amb)/(epsilon*dz);
    else
        //Bottom
        rhof[1] * der(hf[1]) =
        -2*kf[1]*kf[2]/(kf[1]+kf[2]) * (Tf[1]-Tf[2])/(dz^2)
        +rhof[1]*uf[1]/epsilon*(h_in-hf[1])/dz
        -hv[1]*(Tf[1] - Ts[1])/epsilon
        -U_bot*(Tf[1]-T_amb)/(epsilon*dz) 
        -U_wall*CN.pi*D_tank*(Tf[1]-T_amb)/(epsilon*A);

        //Middle
        for i in 2:Nz - 1 loop
            rhof[i] * der(hf[i]) = 
            -2*kf[i]*kf[i-1]/(kf[i]+kf[i-1]) * (Tf[i]-Tf[i-1])/(dz^2)
            -2*kf[i]*kf[i+1]/(kf[i]+kf[i+1]) * (Tf[i]-Tf[i+1])/(dz^2)
            +rhof[i]*uf[i]/epsilon*(hf[i-1]-hf[i])/dz
            -hv[i]*(Tf[i]-Ts[i])/epsilon
            -U_wall*CN.pi*D_tank*(Tf[i]-T_amb)/(epsilon*A);
        end for;

        //Top
        rhof[Nz] * der(hf[Nz]) = 
        -2*kf[Nz-1]*kf[Nz]/(kf[Nz-1]+kf[Nz]) * (Tf[Nz]-Tf[Nz-1])/(dz^2)
        +rhof[Nz]*uf[Nz]/epsilon*(hf[Nz-1]-hf[Nz])/dz
        -hv[Nz]*(Tf[Nz]-Ts[Nz])/epsilon
        -U_wall*CN.pi*D_tank*(Tf[Nz]-T_amb)/(epsilon*A)
        -U_top*(Tf[Nz]-T_amb)/(epsilon*dz);
    end if;

  //Fluid Property evaluation
  for i in 1:Nz loop
    hf[i] = fluid[i].h;
    Tf[i] = fluid[i].T;
    cpf[i] = fluid[i].cp;
    rhof[i] = fluid[i].rho;
    kf[i] = fluid[i].k;
    muf[i] = fluid[i].mu;
    uf[i] = m_flow / (rhof[i] * A * f_area);
  end for;

  //Convection Equations
  for i in 1:Nz loop
    if abs(m_flow) > 1e-12 then
      Re[i] = rhof[i] * abs(uf[i]) * ds / muf[i];
      Pr[i] = cpf[i] * muf[i] / kf[i];
      Nu[i] = 2.0 + 1.1 * (Re[i] ^ 0.6) * (Pr[i] ^ (1 / 3)); //Wakao and Kaguei
    else
      Re[i] = 0;
      Pr[i] = 0;
      Nu[i] = 2.0;
    end if;
    hv[i] = 6.0 * (1.0 - epsilon) * Nu[i] * kf[i] / (ds * ds); //Note that filler surface area correction factor is applied elsewhere.
  end for;

  //Particle energy balance
  cps * der(Ts[1]) = hv[1] * (Tf[1] - Ts[1]) / ((1 - epsilon) * rhos) + ks / rhos * (2*(Ts[2] - Ts[1])) / (dz*dz);
  for i in 2:Nz-1 loop
    cps * der(Ts[i]) = hv[i] * (Tf[i] - Ts[i]) / ((1 - epsilon) * rhos) + ks / rhos * (Ts[i+1] - 2*Ts[i] + Ts[i-1]) / (dz*dz);
  end for;
  cps * der(Ts[Nz]) = hv[Nz] * (Tf[Nz] - Ts[Nz]) / ((1 - epsilon) * rhos) + ks / rhos * (2*(Ts[Nz-1] - Ts[Nz])) / (dz*dz);

  //Calculated Pumping losses
  for i in 1:Nz loop
    p_drop[i] = dz*(((600*((1-epsilon)^2)*muf[i]*abs(m_flow))/((epsilon^3)*(ds^2)*rhof[i]*CN.pi*(D_tank^2)))+((28*(1-epsilon)*(m_flow^2))/((epsilon^3)*ds*rhof[i]*CN.pi*CN.pi*(D_tank^4))));
  end for;

  p_drop_total = sum(p_drop);

  for i in 1:Nz loop
    if Tf[i] >= T_min and Ts[i] >= T_min then
      der(Ei[i]) = rhof[i]*A*dz*epsilon*der(hf[i]) + rhos*A*dz*(1-epsilon)*cps*der(Ts[i]);
    else
      der(Ei[i]) = 0;
    end if;
  end for;
  E = sum(Ei);

annotation(
    experiment(StopTime = 172800, StartTime = 0, Tolerance = 1e-6, Interval = 60),
    Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false)),
    Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false)),
    Documentation(info =
        "<html>
        <ul>
        <li> <i>Dec 2020</i> by Z. Kee:<br> Resealed first version. </li>
        <li> <i>Feb 2026</i> by A. Fontalvo:<br> Simplification for verification purposes. </li>
        </ul>
        </html>"));
end Section1D;

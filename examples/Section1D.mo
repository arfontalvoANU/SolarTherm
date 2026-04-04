within examples;
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
  parameter SI.Power P_charging = 10e6;
  parameter SI.MassFlowRate m_flow_case = 13.7860;
  parameter SI.MassFlowRate m_flow_case_check = P_charging/(hf_max - hf_min);

  //Calculated Tank Design Parameters
  parameter SI.Length H_tank = 8.0;
  parameter SI.Diameter D_tank = 3.44;
  parameter Real f_area = 0.95;
  parameter SI.Area A = 0.25 * CN.pi * D_tank^2 "Cross sectional area of tank";

  //Thermal Losses
  parameter SI.Temperature T_amb = 298.15 "Ambient temperature (K)";
  parameter SI.Area A_loss_tank = CN.pi*D_tank*D_tank*0.5 + CN.pi*D_tank*H_tank "Heat loss area (m2)";
  parameter SI.CoefficientOfHeatTransfer U_loss_tank = 0.0 "Heat loss coeff of surfaces (W/m2K)";
  parameter SI.CoefficientOfHeatTransfer U_wall = U_loss_tank "Cylinder wall heat loss coeff (W/m2K)";
  parameter SI.CoefficientOfHeatTransfer U_top = U_loss_tank "Top circle heat loss coeff (W/m2K)";
  parameter SI.CoefficientOfHeatTransfer U_bot = U_loss_tank "Bottom circle heat loss coeff (W/m2K)";

  //Inititalize temperature and enthalpy profile
  parameter SI.Temperature Tf_start[Nz] = fill(T_start,Nz);
  parameter SI.Temperature Ts_start[Nz] = fill(T_start, Nz);
  parameter SI.SpecificEnthalpy hf_start[Nz] = fill(Medium.h_Tf(T_start, 0.0), Nz) "Defaults to uniform";
  //Property bounds
    //Fluid
  parameter SI.SpecificEnthalpy hf_min = Medium.h_Tf(T_min,0) "Starting enthalpy of the HTF";
  parameter SI.SpecificEnthalpy hf_max = Medium.h_Tf(T_max,0) "Starting enthalpy of the HTF";
  parameter SI.Density rhof_min = Medium.rho_Tf(T_min,0);
  parameter SI.Density rhof_max = Medium.rho_Tf(T_max,0);

    //Discretization
  parameter SI.Length dz = H_tank / Nz "discretization vertical length of fluid";

  parameter Integer Nz = 400 "Number of finite volume elements in fluid";

  //Initialise Fluid Array
  parameter SI.Length z[Nz] = linspace(0, H_tank, Nz);
  SI.Temperature Tf[Nz] "(K)";
  SI.SpecificEnthalpy hf[Nz](start = hf_start) "J/kg";

  //Inlet and outlet enthalpies and temperatures
  SI.SpecificEnthalpy h_in "Inlet Enthalpy depends on mass flow direction (J/kg)";

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
  Medium.State fluid[Nz]"Fluid object array";//(each h_start = hf_min) 

  SI.Energy E_f[Nz];
  SI.Energy E_s[Nz];
  SI.Power P_f[Nz];
  SI.Power P_s[Nz];
  SI.Power Pf;
  SI.Power Ps;

initial equation
  for i in 1:Nz loop
    fluid[i].h = hf_start[i];
  end for;
  m_flow = -m_flow_case;
  h_in = hf_max;

algorithm
    when Tf[1] >= T_stop_charging then //Start discharging
        m_flow := +m_flow_case;
        h_in := hf_min;
    end when;

equation

    if m_flow < 0 then
        //Bottom
//        der(rhof[1]) * hf[1] + 
        rhof[1] * der(hf[1]) =
        -2*kf[1]*kf[2]/(kf[1]+kf[2]) * (Tf[1]-Tf[2])/(dz^2)
        +rhof[1]*uf[1]/epsilon*(hf[1]-hf[2])/dz
        -hv[1]*(Tf[1] - Ts[1])/epsilon
        -U_bot*(Tf[1]-T_amb)/(epsilon*dz) 
        -U_wall*CN.pi*D_tank*(Tf[1]-T_amb)/(epsilon*A);

        //Middle
        for i in 2:Nz - 1 loop
//            der(rhof[i]) * hf[i] + 
            rhof[i] * der(hf[i]) =
            -2*kf[i]*kf[i-1]/(kf[i]+kf[i-1]) * (Tf[i]-Tf[i-1])/(dz^2)
            -2*kf[i]*kf[i+1]/(kf[i]+kf[i+1]) * (Tf[i]-Tf[i+1])/(dz^2)
            +rhof[i]*uf[i]/epsilon*(hf[i]-hf[i+1])/dz
            -hv[i]*(Tf[i]-Ts[i])/epsilon
            -U_wall*CN.pi*D_tank*(Tf[i]-T_amb)/(epsilon*A);
        end for;

        //Top
//        der(rhof[Nz]) * hf[Nz] + 
        rhof[Nz] * der(hf[Nz]) =
        -2*kf[Nz-1]*kf[Nz]/(kf[Nz-1]+kf[Nz]) * (Tf[Nz]-Tf[Nz-1])/(dz^2)
        +rhof[Nz]*uf[Nz]/epsilon*(hf[Nz]-h_in)/dz
        -hv[Nz]*(Tf[Nz]-Ts[Nz])/epsilon
        -U_wall*CN.pi*D_tank*(Tf[Nz]-T_amb)/(epsilon*A)
        -U_top*(Tf[Nz]-T_amb)/(epsilon*dz);
    else
        //Bottom
//        der(rhof[1]) * hf[1] + 
        rhof[1] * der(hf[1]) =
        -2*kf[1]*kf[2]/(kf[1]+kf[2]) * (Tf[1]-Tf[2])/(dz^2)
        +rhof[1]*uf[1]/epsilon*(h_in-hf[1])/dz
        -hv[1]*(Tf[1] - Ts[1])/epsilon
        -U_bot*(Tf[1]-T_amb)/(epsilon*dz) 
        -U_wall*CN.pi*D_tank*(Tf[1]-T_amb)/(epsilon*A);

        //Middle
        for i in 2:Nz - 1 loop
//            der(rhof[i]) * hf[i] + 
            rhof[i] * der(hf[i]) = 
            -2*kf[i]*kf[i-1]/(kf[i]+kf[i-1]) * (Tf[i]-Tf[i-1])/(dz^2)
            -2*kf[i]*kf[i+1]/(kf[i]+kf[i+1]) * (Tf[i]-Tf[i+1])/(dz^2)
            +rhof[i]*uf[i]/epsilon*(hf[i-1]-hf[i])/dz
            -hv[i]*(Tf[i]-Ts[i])/epsilon
            -U_wall*CN.pi*D_tank*(Tf[i]-T_amb)/(epsilon*A);
        end for;

        //Top
//        der(rhof[Nz]) * hf[Nz] + 
        rhof[Nz] * der(hf[Nz]) = 
        -2*kf[Nz-1]*kf[Nz]/(kf[Nz-1]+kf[Nz]) * (Tf[Nz]-Tf[Nz-1])/(dz^2)
        +rhof[Nz]*uf[Nz]/epsilon*(hf[Nz-1]-hf[Nz])/dz
        -hv[Nz]*(Tf[Nz]-Ts[Nz])/epsilon
        -U_wall*CN.pi*D_tank*(Tf[Nz]-T_amb)/(epsilon*A)
        -U_top*(Tf[Nz]-T_amb)/(epsilon*dz);
    end if;

  //Fluid Property evaluation SolarSalt
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
    P_f[i] = rhof[i] * der(hf[i]) * dz * A * epsilon;
    P_s[i] = rhos * cps * der(Ts[i]) * dz * A * (1-epsilon);
    der(E_f[i]) = P_f[i];
    der(E_s[i]) = P_s[i];
  end for;
  Pf = sum(P_f);
  Ps = sum(P_s);

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

annotation(
    experiment(StopTime = 5000, StartTime = 0, Tolerance = 1e-6, Interval = 1),
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

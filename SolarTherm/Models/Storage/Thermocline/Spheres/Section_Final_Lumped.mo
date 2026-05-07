within SolarTherm.Models.Storage.Thermocline.Spheres;
model Section_Final_Lumped "Heat transfer model of thermocline tank with spherical fillers"
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  import Tables = Modelica.Blocks.Tables;

  //Initialize Material Packages
  replaceable package Fluid_Package = SolarTherm.Materials.Air_CoolProp_Table_1bar constrainedby SolarTherm.Materials.PartialMaterial "Fluid Package";
  replaceable package Filler_Package = SolarTherm.Materials.Steatite constrainedby SolarTherm.Materials.PartialMaterial "Filler Package";

  //Fluid Material States
  Fluid_Package.State fluid_in "Model which calculates properties at inlet of the section";
  Fluid_Package.State fluid_out "Model which calculates properties at outlet of the section";

  //Interfacial heat transfer Settings
  parameter Integer Correlation = 1 "1=WakaoKaguei, 2=MelissariArgyropolus, 3=Conservative, 4=Bellan, 5=Laminar, 6=Laminar+Turbulent, 7 = Nie";

  //Height offset for plotting purposes
  parameter SI.Length z_offset = 0.0 "Amount of height offset if there is a tank below it";

  //Tank Design parameters
  parameter SI.Energy E_max = 144e9 "Design storage capacity";
  parameter Real epsilon = 0.4 "Porosity";
  parameter Real ds = 0.02 "Diameter of sphere (particle) (m)";

  //Temperature Bounds
  parameter SI.Temperature T_min = 293 "Design cold Temperature of everything in the tank (K)";
  parameter SI.Temperature T_max = 823 "Design hot Temperature of everything in the tank (K)";
  parameter SI.Temperature T_start = 293 "Initial (uniform) temperature of all components (K), defaults to T_min";

  //Calculated Tank Design Parameters
  parameter SI.Length H_tank = 1.2;
  parameter SI.Diameter D_tank = 0.148;
  parameter Real f_area = 0.95;
  parameter SI.Area A = 0.25 * CN.pi * D_tank * D_tank "Cross sectional area of tank";
  parameter SI.Volume V = A*H_tank;
  parameter SI.Mass ms = (1-epsilon)*V*rhos;
  parameter SI.Mass mf = epsilon*V*rhof_avg;
  final parameter SI.Energy Emax = mf*(hf_max-hf_min) + ms*cps*(T_max-T_min);

  //Thermal Losses
  SI.Temperature T_amb;
  parameter SI.Area A_loss_tank = CN.pi*D_tank*D_tank*0.5 + CN.pi*D_tank*H_tank "Heat loss area (m2)";
  parameter SI.CoefficientOfHeatTransfer U_wall = 0.678 "Cylinder wall heat loss coeff (W/m2K)";
  parameter SI.CoefficientOfHeatTransfer U_top = 0.0 "Top circle heat loss coeff (W/m2K)";
  parameter SI.CoefficientOfHeatTransfer U_bot = 0.0 "Bottom circle heat loss coeff (W/m2K)";

  //Inititalize temperature and enthalpy profile
  parameter SI.Temperature Tf_start[Nz] = fill(T_start,Nz);
  parameter SI.Temperature Ts_start[Nz] = fill(T_start, Nz);
  parameter SI.Temperature T_e_start = T_start;
  parameter SI.SpecificEnthalpy hf_start[Nz] = fill(Fluid_Package.h_Tf(T_start, 0.0), Nz) "Defaults to uniform";
  //parameter SI.SpecificEnthalpy h_p_start[Nz] = fill(Filler_Package.h_Tf(T_start, 0.0), Nz) "Defaults to uniform";
  //Property bounds
    //Fluid
  parameter SI.SpecificEnthalpy hf_min = Fluid_Package.h_Tf(T_min,0) "Starting enthalpy of the HTF";
  parameter SI.SpecificEnthalpy hf_max = Fluid_Package.h_Tf(T_max,0) "Starting enthalpy of the HTF";
  parameter SI.SpecificHeatCapacity cpf_min = Fluid_Package.cp_Tf(T_min,0);
  parameter SI.SpecificHeatCapacity cpf_max = Fluid_Package.cp_Tf(T_max,0);
  parameter SI.SpecificHeatCapacity cpf_avg = (cpf_min + cpf_max) / 2;
  parameter SI.Density rhof_min = Fluid_Package.rho_Tf(T_min,0);
  parameter SI.Density rhof_max = Fluid_Package.rho_Tf(T_max,0);
  parameter SI.Density rhof_avg = (rhof_min + rhof_max) / 2;
  // Filler options
  parameter SI.SpecificHeatCapacity cps = 1068.0 "Filler heat capacity (J/kg/K)";
  parameter SI.Density rhos = 2680.0 "Filler density (kg/m3)";
  parameter SI.ThermalConductivity ks = 2.5 "Filler thermal conductivity (W/m/K)";

    //Discretization
  parameter SI.Length dz = H_tank / Nz "discretization vertical length of fluid";

  parameter Integer Nz = 25 "Number of finite volume elements in fluid";

  //Initialise Fluid Array
  parameter SI.Length z_f[Nz] = Z_position(H_tank, Nz) .+ z_offset;
  SI.Temperature Tf[Nz] "(K)";
  SI.SpecificEnthalpy hf[Nz](start = hf_start) "J/kg";

  //Plotting
  parameter Real ZDH[Nz] = Relative_Tank_Axes(H_tank, Nz) "Non-dimensional tank vertical axis";

  //Operational Controls
  Integer State(start = 2) "operational state 2=standby, 3=discharge, 1=charge";

  //Inlet and outlet enthalpies and temperatures
  SI.SpecificEnthalpy h_in "Inlet Enthalpy depends on mass flow direction (J/kg)";
  SI.SpecificEnthalpy h_out "Outlet Enthalpy depends on mass flow direction (J/kg)";
  SI.Temperature T_in "Inlet Temperature depends on mass flow direction";
  SI.Temperature T_out "Outlet Temperature depends on mass flow direction";

  //Mass flow rates and superficial velocity
  SI.MassFlowRate m_flow(start=0.0) "kg/s";
  SI.Velocity uf_avg "Average fluid velocity in packed bed (m/s)";

  //Analytics
  SI.Energy Ei[Nz] "Energy stored at each section";
  SI.Energy E(start = 0.0) "Make sure the tank starts from T_min for this to be correct";
  Real Level(start = 0.0) "Tank energy charging level (0-1)";
  SI.HeatFlowRate Q_loss_total "Heat loss from the entire surface area";

  //Initialise Particle
  SI.Temperature Ts[Nz](start = Ts_start) "Temperature of particle elements";

  //Calculated Pumping Losses
  parameter Real eff_pump = 1.0 "Pump electricity to work efficiency";
  SI.Pressure p_drop_total "Sum of all pressure drops";
  SI.Power W_loss_pump "losses due to pressure drop";

  //Filler Surface Area Correction
  parameter Real f_surface = 1.0 "Don't touch this";

  //Filler mass-liquid fraction
  //Real f_p[Nz](start=fill(0.0,Nz)) "Mass liquid fraction of filler";

  //Measured outlet temperature
  Real T_outlet_degC "Outlet temperature in degrees Celcius";

//protected
  //Convection Properties
  Real Re[Nz] "Reynolds";
  Real Pr[Nz] "Prandtl";
  Real Nu[Nz] "Nusselt";
  Real hv[Nz] "Volumetric heat transfer coeff (W/m3K)";

  //Filler Properties
  //SI.SpecificEnthalpy h_p[Nz](start=h_p_start) "J/kg";
  //SI.ThermalConductivity k_p[Nz] "W/mK";

  //Filler Geometry
  parameter Real N_spheres_total = (Nz * 6 * (1-epsilon) * A * dz / (CN.pi * (ds^3))) "Total number of spheres in the tank";

  //Pressure Drop
  SI.Pressure p_drop[Nz] "Pressure drop across each mesh element";

  //Thermal Losses
  SI.HeatFlowRate Q_loss_wall[Nz] "Heat loss from the wall";
  SI.HeatFlowRate Q_loss_top "Heat loss from the top";
  SI.HeatFlowRate Q_loss_bot "Heat loss from the bottom";

  //Fluid Properties
  SI.ThermalConductivity kf[Nz] "W/mK";
  SI.DynamicViscosity muf[Nz] "Pa.s";
  SI.SpecificHeatCapacity cpf[Nz] "J/kg.K";
  Fluid_Package.State fluid[Nz]"Fluid object array";//(each h_start = hf_min) 

  //Try filler state "Remove this if using function-based calculation"
  //Filler_Package.State filler[Nz] "Filler object array";
  Real der_hf[Nz] "Rate of change of specific enthalpy of fluid";

algorithm
  //Fluid equations
  if State == 1 then
  //Charging (Mass flows top to bottom)
  //Bottom Charging Fluid Node
    der_hf[1] :=
    (-2*kf[1]*kf[2]/(kf[1]+kf[2]) * (Tf[1]-Tf[2])/(dz^2)
        +rhof_avg*uf_avg/epsilon*(hf[1]-hf[2])/dz
        -hv[1]*(Tf[1] - Ts[1])/epsilon
        -U_bot*(Tf[1]-T_amb)/(epsilon*dz)
        -U_wall*CN.pi*D_tank*(Tf[1]-T_amb)/(epsilon*A))/(rhof_avg);
    h_out := hf[1];
  //End Bottom Charging Fluid Node
  //Middle Charging Fluid Nodes
    for i in 2:Nz - 1 loop
      der_hf[i] := 
      (-2*kf[i]*kf[i-1]/(kf[i]+kf[i-1]) * (Tf[i]-Tf[i-1])/(dz^2)
            -2*kf[i]*kf[i+1]/(kf[i]+kf[i+1]) * (Tf[i]-Tf[i+1])/(dz^2)
            +rhof_avg*uf_avg/epsilon*(hf[i]-hf[i+1])/dz
            -hv[i]*(Tf[i]-Ts[i])/epsilon
            -U_wall*CN.pi*D_tank*(Tf[i]-T_amb)/(epsilon*A))/(rhof_avg);
    end for;
  //End Middle Charging Fluid Nodes
  //Top Charging Fluid Node
    der_hf[Nz] := 
    (-2*kf[Nz-1]*kf[Nz]/(kf[Nz-1]+kf[Nz]) * (Tf[Nz]-Tf[Nz-1])/(dz^2)
        +rhof_avg*uf_avg/epsilon*(hf[Nz]-h_in)/dz
        -hv[Nz]*(Tf[Nz]-Ts[Nz])/epsilon
        -U_wall*CN.pi*D_tank*(Tf[Nz]-T_amb)/(epsilon*A)
        -U_top*(Tf[Nz]-T_amb)/(epsilon*dz))/(rhof_avg);
  //End Top Charging Fluid Node
  else
  //Discharge (Mass flows bottom to top)
  //Bottom Discharge Node
    der_hf[1] :=
    (-2*kf[1]*kf[2]/(kf[1]+kf[2]) * (Tf[1]-Tf[2])/(dz^2)
        +rhof_avg*uf_avg/epsilon*(h_in-hf[1])/dz
        -hv[1]*(Tf[1] - Ts[1])/epsilon
        -U_bot*(Tf[1]-T_amb)/(epsilon*dz) 
        -U_wall*CN.pi*D_tank*(Tf[1]-T_amb)/(epsilon*A))/(rhof_avg);
  //End Bottom Discharge Node
  //Middle Discharge Nodes
    for i in 2:Nz - 1 loop
      der_hf[i] :=
      (-2*kf[i]*kf[i-1]/(kf[i]+kf[i-1]) * (Tf[i]-Tf[i-1])/(dz^2)
            -2*kf[i]*kf[i+1]/(kf[i]+kf[i+1]) * (Tf[i]-Tf[i+1])/(dz^2)
            +rhof_avg*uf_avg/epsilon*(hf[i-1]-hf[i])/dz
            -hv[i]*(Tf[i]-Ts[i])/epsilon
            -U_wall*CN.pi*D_tank*(Tf[i]-T_amb)/(epsilon*A))/(rhof_avg);
    end for;
  //End Middle Discharge Nodes
  //Top Discharge Node
    der_hf[Nz] :=
    (-2*kf[Nz-1]*kf[Nz]/(kf[Nz-1]+kf[Nz]) * (Tf[Nz]-Tf[Nz-1])/(dz^2)
        +rhof_avg*uf_avg/epsilon*(hf[Nz-1]-hf[Nz])/dz
        -hv[Nz]*(Tf[Nz]-Ts[Nz])/epsilon
        -U_wall*CN.pi*D_tank*(Tf[Nz]-T_amb)/(epsilon*A)
        -U_top*(Tf[Nz]-T_amb)/(epsilon*dz))/(rhof_avg);
    h_out := hf[Nz];
  end if;

initial equation
  for i in 1:Nz loop
    Tf[i] = T_start;
    Ts[i] = T_start;
  end for;

equation
  for i in 1:Nz loop
    der_hf[i] = der(hf[i]);
  end for;

  //Determine which operational state: In this version, standby and discharge are lumped.
  if m_flow < 0.0 then //mass is flowing downwards so charging
    State = 1;
  else //mass is flowing upwards so discharging
    State = 3;
  end if;

  //Fluid inlet and outlet properties
  fluid_in.h = h_in;
  fluid_out.h = h_out;
  fluid_in.T = T_in;
  fluid_out.T = T_out;

  //Fluid Property evaluation SolarSalt
  for i in 1:Nz loop
    hf[i] = fluid[i].h;
    Tf[i] = fluid[i].T;
    cpf[i] = fluid[i].cp;
    kf[i] = fluid[i].k;
    muf[i] = fluid[i].mu;
  end for;
  uf_avg = m_flow / (rhof_avg * A * f_area);

  //Convection Equations
  for i in 1:Nz loop
    if abs(m_flow) > 1e-12 then //There is actually mass flowing
      Re[i] = rhof_avg * abs(uf_avg) * ds / muf[i]; //Use local superficial velocity
      Pr[i] = cpf[i] * muf[i] / kf[i];
      if Correlation == 1 then 
        Nu[i] = 2 + 1.1 * (Re[i] ^ 0.6) * (Pr[i] ^ (1 / 3)); //Wakao and Kaguei
      elseif Correlation == 2 then
        Nu[i] = 2 + 0.47 * (Re[i] ^ 0.5) * (Pr[i] ^ 0.36);//Use Melissari and Argyropolus
      elseif Correlation == 3 then
        Nu[i] = 2; //Conservative
      elseif Correlation == 4 then
        Nu[i] = 2 + 0.664*(Re[i]^0.5)*(Pr[i]^(1/3)); //Only laminar
      elseif Correlation == 5 then //laminar plus turbulent
        Nu[i] = (1 + 1.5*(1 - epsilon))*(2 + ((0.664*(Re[i]^0.5)*((Pr[i])^(1/3)))^2+((0.037*(Re[i]^0.8)*Pr[i])/(1 + 2.443*(Re[i]^(-0.1))*((Pr[i]^(2/3))-1)))^2)^0.5);
      elseif Correlation == 6 then//Nie
        Nu[i] = 0.052*(((1.0-epsilon)^0.14)/epsilon)*(Re[i]^0.86)*(Pr[i]^(1/3));
      else
        Nu[i] = (2.06/epsilon)*(Re[i]^0.425)*(Pr[i]^(1/3));
      end if;
    else
      Re[i] = 0;
      Pr[i] = 0;
      Nu[i] = 2;
    end if;
    hv[i] = (f_surface)*6*(1 - epsilon) * Nu[i] * kf[i] / (ds^2); //Note that filler surface area correction factor is applied elsewhere.
  end for;
  //Particle energy balance
  cps*der(Ts[1]) = hv[1] * (Tf[1] - Ts[1]) / ((1 - epsilon) * rhos) + ks / rhos * (2*(Ts[2] - Ts[1])) / (dz*dz);
  for i in 2:Nz-1 loop
    cps*der(Ts[i]) = hv[i] * (Tf[i] - Ts[i]) / ((1 - epsilon) * rhos) + ks / rhos * (Ts[i+1] - 2*Ts[i] + Ts[i-1]) / (dz*dz);
  end for;
  cps*der(Ts[Nz]) = hv[Nz] * (Tf[Nz] - Ts[Nz]) / ((1 - epsilon) * rhos) + ks / rhos * (2*(Ts[Nz-1] - Ts[Nz])) / (dz*dz);
  
  //Heat loss calculations, different form than the equations above as they were in terms of rho*dh/dt not m*dh/dt
  Q_loss_top = U_top*CN.pi*D_tank*D_tank*0.25*(Tf[Nz]-T_amb);
  Q_loss_bot = U_bot*CN.pi*D_tank*D_tank*0.25*(Tf[1]-T_amb);
  for i in 2:Nz-1 loop
    Q_loss_wall[i] = U_wall*CN.pi*D_tank*(Tf[i]-T_amb)*dz;
  end for;
  Q_loss_wall[1] = U_wall*CN.pi*D_tank*(Tf[1]-T_amb)*dz;
  Q_loss_wall[Nz] = U_wall*CN.pi*D_tank*(Tf[Nz]-T_amb)*dz;
  Q_loss_total = Q_loss_top + sum(Q_loss_wall) + Q_loss_bot;
  //End heat loss calculations

  //Calculated Pumping losses
  for i in 1:Nz loop
    p_drop[i] = dz*(((600*((1-epsilon)^2)*muf[i]*abs(m_flow))/((epsilon^3)*(ds^2)*rhof_avg*CN.pi*(D_tank^2)))+((28*(1-epsilon)*(m_flow^2))/((epsilon^3)*ds*rhof_avg*CN.pi*CN.pi*(D_tank^4))));
  end for;

  p_drop_total = sum(p_drop);
  W_loss_pump = (abs(m_flow)/rhof_avg)*p_drop_total/eff_pump;

  //Analyics
  for i in 1:Nz loop
    der(Ei[i]) = rhof_avg*A*dz*epsilon*der_hf[i] + rhos*A*dz*(1-epsilon)*cps*der(Ts[i]);
  end for;
  E = sum(Ei);
  Level = E / Emax;

  if m_flow > 1e-3 then //Discharging, outlet is the top
    T_outlet_degC = Tf[Nz] - 273.15;
  elseif m_flow < -1e-3 then //Charging, outlet is the bottom
    T_outlet_degC = Tf[1] - 273.15;
  else //No flow, output reference temperature
    T_outlet_degC = 25;
  end if;

  annotation (Documentation(revisions ="<html>
    <p>By Zebedee Kee on 03/12/2020</p>
    </html>",info="<html>
    <p>This model contains the heat-transfer calculations of a thermocline packed bed storage tank with spherical filler geometry. 
    This model does not contain any fluid connectors, for the CSP component with connectors, see Thermocline_Spheres_SingleTank. 
    Variables fluid_top and fluid_bot provides the enthalpy-temperature relationship of the fluid material. 
    Depending on whether m_flow is positive (discharging, fluid flowing upwards) or negative (charging, fluid flowing downwards), 
    the charging/discharging equations are applied. In this iteration of the model, discharging and standby are lumped into one state.</p>
    </html>"));
end Section_Final_Lumped;

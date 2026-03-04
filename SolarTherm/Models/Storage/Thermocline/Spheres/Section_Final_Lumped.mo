within SolarTherm.Models.Storage.Thermocline.Spheres;
model Section_Final_Lumped "Heat transfer model of thermocline tank with spherical fillers"
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  import Tables = Modelica.Blocks.Tables;

  //Initialize Material Packages
  replaceable package Fluid_Package = SolarTherm.Materials.Air_CoolProp_Table_1bar constrainedby SolarTherm.Materials.PartialMaterial "Fluid Package";
  replaceable package Filler_Package = SolarTherm.Materials.Steatite constrainedby SolarTherm.Materials.PartialMaterial "Filler Package";
  replaceable package Tank_Package =  SolarTherm.Materials.SS316L constrainedby SolarTherm.Materials.PartialMaterial "Tank Package (steel shell)";
  replaceable package Encapsulation_Package = Filler_Package constrainedby SolarTherm.Materials.PartialMaterial "Encapsulation Package, default is the same as Filler package, effectively no encapsulation";

  //Fluid Material States
  Fluid_Package.State fluid_in "Model which calculates properties at inlet of the section";
  Fluid_Package.State fluid_out "Model which calculates properties at outlet of the section";

  //Interfacial heat transfer Settings
  parameter Integer Correlation = 1 "1=WakaoKaguei, 2=MelissariArgyropolus, 3=Conservative, 4=Bellan, 5=Laminar, 6=Laminar+Turbulent, 7 = Nie";

  //Height offset for plotting purposes
  parameter SI.Length z_offset = 0.0 "Amount of height offset if there is a tank below it";

  //Tank Design parameters
  parameter SI.Energy E_max = 144e9 "Design storage capacity";
  parameter Real eta = 0.4 "Porosity";
  parameter Real d_p = 0.02 "Diameter of sphere (particle) (m)";
  parameter Real t_e = d_p/(2.0*N_p) "Thickness of encapsulation, default is such that it is at a value that preserves equidistant radii discretizations (m)";

  //Temperature Bounds
  parameter SI.Temperature T_min = 293 "Design cold Temperature of everything in the tank (K)";
  parameter SI.Temperature T_max = 823 "Design hot Temperature of everything in the tank (K)";
  parameter SI.Temperature T_start = 293 "Initial (uniform) temperature of all components (K), defaults to T_min";

  //Calculated Tank Design Parameters
  parameter SI.Length H_tank = 1.2;
  parameter SI.Diameter D_tank = 0.148;
  parameter SI.Area A = 0.25 * CN.pi * D_tank * D_tank "Cross sectional area of tank";

  //Thermal Losses
  SI.Temperature T_amb;
  parameter SI.Area A_loss_tank = CN.pi*D_tank*D_tank*0.5 + CN.pi*D_tank*H_tank "Heat loss area (m2)";
  parameter SI.CoefficientOfHeatTransfer U_loss_tank = 0.1 "Heat loss coeff of surfaces (W/m2K)";
  parameter SI.CoefficientOfHeatTransfer U_wall = U_loss_tank "Cylinder wall heat loss coeff (W/m2K)";
  parameter SI.CoefficientOfHeatTransfer U_top = U_loss_tank "Top circle heat loss coeff (W/m2K)";
  parameter SI.CoefficientOfHeatTransfer U_bot = U_loss_tank "Bottom circle heat loss coeff (W/m2K)";

  //Inititalize temperature and enthalpy profile
  parameter SI.Temperature T_f_start[N_f] = fill(T_start,N_f);
  parameter SI.Temperature T_p_start[N_f] = fill(T_start, N_f);
  parameter SI.Temperature T_e_start = T_start;
  parameter SI.SpecificEnthalpy h_f_start[N_f] = fill(Fluid_Package.h_Tf(T_start, 0.0), N_f) "Defaults to uniform";
  parameter SI.SpecificEnthalpy h_p_start[N_f] = fill(Filler_Package.h_Tf(T_start, 0.0), N_f) "Defaults to uniform";
  //Property bounds
    //Fluid
  parameter SI.SpecificEnthalpy h_f_min = Fluid_Package.h_Tf(T_min,0) "Starting enthalpy of the HTF";
  parameter SI.SpecificEnthalpy h_f_max = Fluid_Package.h_Tf(T_max,0) "Starting enthalpy of the HTF";
  parameter SI.Density rho_f_min = Fluid_Package.rho_Tf(T_min,0);
  parameter SI.Density rho_f_max = Fluid_Package.rho_Tf(T_max,0);
  parameter SI.Density rho_f_avg = (rho_f_min + rho_f_max) / 2;
    //Filler
  parameter SI.SpecificEnthalpy h_p_max = Filler_Package.h_Tf(T_max, 1.0);
  parameter SI.SpecificEnthalpy h_p_min = Filler_Package.h_Tf(T_min, 0.0);
  parameter SI.Density rho_p_min = Filler_Package.rho_Tf(T_min, 0.0);
  parameter SI.Density rho_p_max = Filler_Package.rho_Tf(T_max, 1.0);
  parameter SI.Density rho_p = min(rho_p_min,rho_p_max) "kg/m3";
    //Encapsulation
  parameter SI.SpecificEnthalpy h_e_max = Encapsulation_Package.h_Tf(T_max, 1.0);
  parameter SI.SpecificEnthalpy h_e_min = Encapsulation_Package.h_Tf(T_min, 0.0);
  parameter SI.Density rho_e_min = Encapsulation_Package.rho_Tf(T_min, 0.0);
  parameter SI.Density rho_e_max = Encapsulation_Package.rho_Tf(T_max, 1.0);
  parameter SI.Density rho_e = min(rho_e_min,rho_e_max) "kg/m3";

    //Discretization
  parameter SI.Length dz = H_tank / N_f "discretization vertical length of fluid";
  parameter SI.Length dr[N_p] = cat(1,fill(0.5 * ((d_p - 2.0*t_e) / (N_p - 1)),N_p-1),{t_e}) "radial thickness of each particle discretization, with last one being the encapsulation"; 

  parameter Integer N_f = 25 "Number of finite volume elements in fluid";
  parameter Integer N_p = 10 "Number of finite volume elements in filler, including encapsulation";

  //Initialise Fluid Array
  parameter SI.Length z_f[N_f] = Z_position(H_tank, N_f) .+ z_offset;
  SI.Temperature T_f[N_f] "(K)";
  SI.SpecificEnthalpy h_f[N_f](start = h_f_start) "J/kg";

  //Plotting
  parameter Real ZDH[N_f] = Relative_Tank_Axes(H_tank, N_f) "Non-dimensional tank vertical axis";

  //Operational Controls
  Integer State(start = 2) "operational state 2=standby, 3=discharge, 1=charge";

  //Inlet and outlet enthalpies and temperatures
  SI.SpecificEnthalpy h_in "Inlet Enthalpy depends on mass flow direction (J/kg)";
  SI.SpecificEnthalpy h_out "Outlet Enthalpy depends on mass flow direction (J/kg)";
  SI.Temperature T_in "Inlet Temperature depends on mass flow direction";
  SI.Temperature T_out "Outlet Temperature depends on mass flow direction";

  //Mass flow rates and superficial velocity
  SI.MassFlowRate m_flow(start=0.0) "kg/s";
  SI.Velocity u_avg "Average fluid velocity in packed bed (m/s)";
  SI.Velocity u_f[N_f] "Fluid velocity in packed bed (m/s)";

  //Analytics
  SI.Energy E_stored(start = 0.0) "Make sure the tank starts from T_min for this to be correct";
  Real Level(start = 0.0) "Tank energy charging level (0-1)";
  SI.HeatFlowRate Q_loss_total "Heat loss from the entire surface area";

  //Initialise Particle
  SI.Temperature T_p[N_f](start = T_p_start) "Temperature of particle elements";

  //Calculated Pumping Losses
  parameter Real eff_pump = 1.0 "Pump electricity to work efficiency";
  SI.Pressure p_drop_total "Sum of all pressure drops";
  SI.Power W_loss_pump "losses due to pressure drop";

  //Filler Surface Area Correction
  parameter Real f_surface = 1.0 "Don't touch this";

  parameter SI.Length r_p[N_p] = cat(1,Particle_Radii(d_p-2*t_e,N_p-1),{(d_p/2)-(t_e/2)}) "Radii of each particle element centre";
  //Filler mass-liquid fraction
  Real f_p[N_f](start=fill(0.0,N_f)) "Mass liquid fraction of filler";

  //Measured outlet temperature
  Real T_outlet_degC "Outlet temperature in degrees Celcius";

//protected
  //Convection Properties
  Real Pe[N_f] "Peclet Number";
  Real Bi[N_f] "Biot Number";
  Real Re[N_f] "Reynolds";
  Real Pr[N_f] "Prandtl";
  Real Nu[N_f] "Nusselt";
  Real h_v[N_f] "Volumetric heat transfer coeff (W/m3K)";

  //Filler Properties
  SI.SpecificEnthalpy h_p[N_f](start=h_p_start) "J/kg";
  SI.ThermalConductivity k_p[N_f] "W/mK";

  //Filler Geometry
  parameter Real N_spheres_total = (N_f * 6 * (1-eta) * A * dz / (CN.pi * (d_p^3))) "Total number of spheres in the tank";

  //Pressure Drop
  SI.Pressure p_drop[N_f] "Pressure drop across each mesh element";

  //Thermal Losses
  SI.HeatFlowRate Q_loss_wall[N_f] "Heat loss from the wall";
  SI.HeatFlowRate Q_loss_top "Heat loss from the top";
  SI.HeatFlowRate Q_loss_bot "Heat loss from the bottom";

  //Fluid Properties
  SI.ThermalConductivity k_f[N_f] "W/mK";
  SI.DynamicViscosity mu_f[N_f] "Pa.s";
  SI.SpecificHeatCapacity c_pf[N_f] "J/kgK";
  SI.Density rho_f[N_f] "kg/m3";
  Fluid_Package.State fluid[N_f]"Fluid object array";//(each h_start = h_f_min) 

  //Try filler state "Remove this if using function-based calculation"
  Filler_Package.State filler[N_f] "Filler object array";
  Real der_h_f[N_f] "Rate of change of specific enthalpy of fluid";

  parameter Real C_ax = 0.2;
  SI.DiffusionCoefficient D_ax[N_f];
  SI.ThermalConductivity k_f_eff[N_f] "W/mK";
  SI.Velocity u_0[N_f] "Superficial fluid velocity in packed bed (m/s)";

algorithm
  //Fluid equations
  if State == 1 then
  //Charging (Mass flows top to bottom)
  //Bottom Charging Fluid Node
    der_h_f[1] :=
    ( (-2.0*k_f_eff[1]*k_f_eff[2])*(T_f[1]-T_f[2])/((k_f_eff[1]+k_f_eff[2])*dz*dz)
    + (rho_f[1]*u_f[1])*(h_f[1]-h_f[2])/dz
    - h_v[1]*(T_f[1] - T_p[1])/eta
    - U_bot*(T_f[1]-T_amb)/(eta*dz) 
    - U_wall*CN.pi*D_tank*(T_f[1]-T_amb)/(eta*A) ) / (rho_f[1]);
  
    h_out := h_f[1];
  //End Bottom Charging Fluid Node
  //Middle Charging Fluid Nodes
    for i in 2:N_f - 1 loop
      der_h_f[i] := 
      ( 2.0*k_f_eff[i - 1]*k_f_eff[i]*(T_f[i-1]-T_f[i])/((k_f_eff[i-1]+k_f_eff[i])*dz*dz)
      - 2.0*k_f_eff[i]*k_f_eff[i+1]*(T_f[i]-T_f[i + 1])/((k_f_eff[i]+k_f_eff[i+1])*dz*dz)
      + (rho_f[i]*u_f[i])*(h_f[i]-h_f[i+1])/dz
      - h_v[i]*(T_f[i]-T_p[i])/eta
      - U_wall*CN.pi*D_tank*(T_f[i]-T_amb)/(eta*A) )/ (rho_f[i]);
    end for;
  //End Middle Charging Fluid Nodes
  //Top Charging Fluid Node
    der_h_f[N_f] := 
    (2.0*k_f_eff[N_f-1]*k_f_eff[N_f]*(T_f[N_f-1]-T_f[N_f])/((k_f_eff[N_f-1]+k_f_eff[N_f])*dz*dz)
    + (rho_f[N_f]*u_f[N_f])*(h_f[N_f]-h_in)/dz
    - h_v[N_f]*(T_f[N_f]-T_p[N_f])/eta
    - U_wall*CN.pi*D_tank*(T_f[N_f]-T_amb)/(eta*A)
    - U_top*(T_f[N_f]-T_amb)/(eta*dz) ) / (rho_f[N_f]);
  //End Top Charging Fluid Node
  else
  //Discharge (Mass flows bottom to top)
  //Bottom Discharge Node
    der_h_f[1] :=
    (-2.0*k_f_eff[1]*k_f_eff[2]*(T_f[1]-T_f[2])/((k_f_eff[1]+k_f_eff[2])*dz*dz)
    + (rho_f[1]*u_f[1])*(h_in-h_f[1])/dz
    - h_v[1]*(T_f[1]-T_p[1])/eta
    - U_bot*(T_f[1]-T_amb)/(eta*dz)
    - U_wall*CN.pi*D_tank*(T_f[1]-T_amb)/(eta*A) )/ (rho_f[1]);
  //End Bottom Discharge Node
  //Middle Discharge Nodes
    for i in 2:N_f - 1 loop
      der_h_f[i] :=
      ( 2.0*k_f_eff[i-1]*k_f_eff[i]*(T_f[i-1]-T_f[i])/((k_f_eff[i-1]+k_f_eff[i])*dz*dz)
      - 2.0*k_f_eff[i]*k_f_eff[i + 1]*(T_f[i]-T_f[i+1])/((k_f_eff[i]+k_f_eff[i+1])*dz*dz)
      + (rho_f[i]*u_f[i])*(h_f[i-1]-h_f[i])/dz
      - h_v[i]*(T_f[i]-T_p[i])/eta
      - U_wall*CN.pi*D_tank*(T_f[i]-T_amb)/(eta*A) ) / (rho_f[i]);
    end for;
  //End Middle Discharge Nodes
  //Top Discharge Node
    der_h_f[N_f] :=
    ( 2.0*k_f_eff[N_f-1]*k_f_eff[N_f]*(T_f[N_f-1]-T_f[N_f])/((k_f_eff[N_f-1]+k_f_eff[N_f])*dz*dz)
    + (rho_f[N_f]*u_f[N_f])*(h_f[N_f-1]-h_f[N_f])/dz
    - h_v[N_f]*(T_f[N_f]-T_p[N_f])/eta
    - U_wall*CN.pi*D_tank*(T_f[N_f]-T_amb)/(eta*A)
    - U_top*(T_f[N_f]-T_amb)/(eta*dz) ) / (rho_f[N_f]);
  
    h_out := h_f[N_f];
  end if;

initial equation
  for i in 1:N_f loop
    fluid[i].h = h_f_start[i];
    filler[i].h = h_p_start[i];
  end for;

equation
  for i in 1:N_f loop
    der_h_f[i] = der(h_f[i]);
  end for;

  //Determine which operational state: In this version, standby and discharge are lumped.
  if m_flow < 0.0 then //mass is flowing downwards so charging
    State = 1;
  else //mass is flowing upwards so discharging
    State = 3;
  end if;

  u_avg = m_flow / (eta * rho_f_avg * A); //positive if flowing upwards (discharge)

  //Fluid inlet and outlet properties
  fluid_in.h = h_in;
  fluid_out.h = h_out;
  fluid_in.T = T_in;
  fluid_out.T = T_out;

  //Fluid Property evaluation SolarSalt
  for i in 1:N_f loop
    h_f[i] = fluid[i].h;
    T_f[i] = fluid[i].T;
    c_pf[i] = fluid[i].cp;
    rho_f[i] = fluid[i].rho;
    k_f[i] = fluid[i].k;
    mu_f[i] = fluid[i].mu;
    u_f[i] = m_flow / (eta * rho_f[i] * A);              // interstitial
    u_0[i] = m_flow / (rho_f[i] * A);                    // superficial
    D_ax[i] = C_ax * abs(u_0[i]) * d_p;                  // m2/s
    k_f_eff[i] = fluid[i].k + rho_f[i] * c_pf[i] * D_ax[i];  // W/(m.K)
  end for;

  //Particle Property evaluation
  for i in 1:N_f loop
    filler[i].h = h_p[i];
    T_p[i] = filler[i].T;
    f_p[i] = filler[i].f;
    k_p[i] = filler[i].k;
  end for;

  //Convection Equations
  for i in 1:N_f loop
    if abs(m_flow) > 1e-12 then //There is actually mass flowing
      //Re[i] = abs(m_flow) / A * d_p / mu_f[i]; //Use local superficial velocity
      Re[i] = rho_f[i] * abs(u_0[i]) * d_p / mu_f[i]; //Use local superficial velocity
      Pr[i] = c_pf[i] * mu_f[i] / k_f[i];
      if Correlation == 1 then 
        Nu[i] = 2.0 + 1.1 * (Re[i] ^ 0.6) * (Pr[i] ^ (1 / 3)); //Wakao and Kaguei
      elseif Correlation == 2 then
        Nu[i] = 2.0 + 0.47 * (Re[i] ^ 0.5) * (Pr[i] ^ 0.36);//Use Melissari and Argyropolus
      elseif Correlation == 3 then
        Nu[i] = 2.0; //Conservative
      elseif Correlation == 4 then
        Nu[i] = 2.0 + 0.664*(Re[i]^0.5)*(Pr[i]^(1/3)); //Only laminar
      elseif Correlation == 5 then //laminar plus turbulent
        Nu[i] = (1.0+1.5*(1.0-eta))*(2.0+((0.664*(Re[i]^0.5)*((Pr[i])^(1.0/3.0)))^2.0+((0.037*(Re[i]^0.8)*Pr[i])/(1.0+2.443*(Re[i]^(-0.1))*((Pr[i]^(2/3))-1)))^2)^0.5);
      elseif Correlation == 6 then//Nie
        Nu[i] = 0.052*(((1.0-eta)^0.14)/eta)*(Re[i]^0.86)*(Pr[i]^(1/3));
      else
        Nu[i] = (2.06/eta)*(Re[i]^0.425)*(Pr[i]^(1/3));
      end if;
    else
      Re[i] = 0;
      Pr[i] = 0;
      Nu[i] = 2.0;
    end if;
    Bi[i] = (Nu[i]*k_f[i])/(6.0*k_p[i]); //Use outermost shell conductivity
    Pe[i] = Re[i]*Pr[i];
    h_v[i] = (f_surface)*6.0 * (1.0 - eta) * Nu[i] * k_f[i] / (d_p * d_p); //Note that filler surface area correction factor is applied elsewhere.
  end for;
  //Particle energy balance
  der(h_p[1]) = h_v[1] * (T_f[1] - T_p[1]) / ((1 - eta) * rho_p) + k_p[1] / rho_p * (2*(T_p[2] - T_p[1])) / (dz*dz);
  for i in 2:N_f-1 loop
    der(h_p[i]) = h_v[i] * (T_f[i] - T_p[i]) / ((1 - eta) * rho_p) + k_p[i] / rho_p * (T_p[i+1] - 2*T_p[i] + T_p[i-1]) / (dz*dz);
  end for;
  der(h_p[N_f]) = h_v[N_f] * (T_f[N_f] - T_p[N_f]) / ((1 - eta) * rho_p) + k_p[N_f] / rho_p * (2*(T_p[N_f-1] - T_p[N_f])) / (dz*dz);
  
  //Heat loss calculations, different form than the equations above as they were in terms of rho*dh/dt not m*dh/dt
  Q_loss_top = U_top*CN.pi*D_tank*D_tank*0.25*(T_f[N_f]-T_amb);
  Q_loss_bot = U_bot*CN.pi*D_tank*D_tank*0.25*(T_f[1]-T_amb);
  for i in 2:N_f-1 loop
    Q_loss_wall[i] = U_wall*CN.pi*D_tank*(T_f[i]-T_amb)*dz;
  end for;
  Q_loss_wall[1] = U_wall*CN.pi*D_tank*(T_f[1]-T_amb)*dz;
  Q_loss_wall[N_f] = U_wall*CN.pi*D_tank*(T_f[N_f]-T_amb)*dz;
  Q_loss_total = Q_loss_top + sum(Q_loss_wall) + Q_loss_bot;
  //End heat loss calculations

  //Calculated Pumping losses
  for i in 1:N_f loop
    p_drop[i] = dz*(((600*((1-eta)^2)*mu_f[i]*abs(m_flow))/((eta^3)*(d_p^2)*rho_f[i]*CN.pi*(D_tank^2)))+((28*(1-eta)*(m_flow^2))/((eta^3)*d_p*rho_f[i]*CN.pi*CN.pi*(D_tank^4))));
  end for;

  p_drop_total = sum(p_drop);
  W_loss_pump = (abs(m_flow)/rho_f_avg)*p_drop_total/eff_pump;

  //Analyics
  der(E_stored) = abs(m_flow) * (h_in - h_out) - Q_loss_total;
  Level = E_stored / E_max;

  if m_flow > 1.0e-3 then //Discharging, outlet is the top
    T_outlet_degC = T_f[N_f] - 273.15;
  elseif m_flow < -1.0e-3 then //Charging, outlet is the bottom
    T_outlet_degC = T_f[1] - 273.15;
  else //No flow, output reference temperature
    T_outlet_degC = 25.0;
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

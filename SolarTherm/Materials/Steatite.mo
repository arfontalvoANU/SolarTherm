within SolarTherm.Materials;
package Steatite "MgO constant cp and density, variable thermal conductivity"
  extends SolarTherm.Materials.PartialMaterial(MM = 40.3044e-4, T_melt = 3105.0, cost = 1.0, year = 2020);
  //Properties of Steatite are defined by constant cp_1, rho_1
  constant SI.SpecificHeatCapacity cp_1 = 1068 "Specific heat capacity of solid (J/kgK)";
  constant SI.Density rho_1 = 2680 "Density of solid (kg/m3)";
  
  //Thermal Conductivity lookup tables
  constant SI.Temperature T_data[19] = {400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0, 1050.0, 1100.0, 1150.0, 1200.0, 1250.0, 1300.0} "Absolute temperature table (K)";
  constant SI.ThermalConductivity k_data[19] = {29.4,26.6,23.8,21.6,19.7,18.3,16.6,15.3,14.3,13.1,12.3,11.7,11.1,10.5,10.0,9.5,8.9,8.3,8.1} "Thermal conductivity table (K)";

  redeclare model State "A model which calculates state and properties"
    SI.SpecificEnthalpy h "Specific Enthalpy wrt 298.15K (J/kg)";
    SI.Temperature T "Temperature (K)";
    Real f "Liquid Mass Fraction";
    SI.Density rho "Density (kg/m3)";
    SI.ThermalConductivity k "Thermal conductivity (W/mK)";
  equation
    T = 298.15 + h/cp_1;
    f = 0.0;
    rho = rho_1;
    k = 2.5;
  end State;

  redeclare function h_Tf "find specific enthalpy from Temperature"
    input SI.Temperature T "Absolute temperature (K)";
    input Real f "Liquid mass fraction";
    output SI.SpecificEnthalpy h "Specific Enthalpy (J/kg)";
  algorithm
    h := cp_1*(T-298.15);
  end h_Tf;
    
  redeclare function rho_Tf "find density from temperature"
    input SI.Temperature T "Absolute temperature (K)";
    input Real f "Liquid mass fraction";
    output SI.Density rho "Density (kg/m3)";
  algorithm
    rho := rho_1;
  end rho_Tf;
  
  function k_Tf "find thermal conductivity from temperature"
    input SI.Temperature T "Absolute temperature (K)";
    input Real f "Liquid mass fraction (-) not used";
    output SI.ThermalConductivity k "Thermal conductivity (W/mK)";
  algorithm
    k := Modelica.Math.Vectors.interpolate(T_data,k_data,T);
  end k_Tf;
  annotation(
    Documentation(revisions= "<html>
    <ul>
    <li> A. Fontalvo Lascano (February 2026): <br> Resealed first version. </li>
    </ul>
    </html>"
  ));
end Steatite;

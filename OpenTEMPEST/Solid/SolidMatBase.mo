within OpenTEMPEST.Solid;
partial package SolidMatBase
  "determines properties of solid material as function of temperature"

  constant OpenTEMPEST.Solid.solidMaterialsRecord data;

  model BaseProperties "model that calls function for T-dependent k,cp,rho with coefficients in 'data' stored in 'SolidMaterialsData'"

  Modelica.SIunits.Density rho "density, kg/m³";
    Modelica.SIunits.SpecificHeatCapacity cp "heat capacity, J/kg/K";
    Modelica.SIunits.ThermalConductivity k_long "heat conductivity in longitudinal direction of layering";
    Modelica.SIunits.ThermalConductivity k_trans "heat conductivity in transversal direction of layering";
    input Modelica.SIunits.Temperature T "temperature, K";
    parameter Modelica.SIunits.ThermalConductivity kCustom_trans = 1 "heat conductivity in transversal direction of layering";
    parameter Modelica.SIunits.ThermalConductivity kCustom_long = 1 "heat conductivity in longitudinal direction of layering";
    parameter Modelica.SIunits.SpecificHeatCapacity cpCustom = 1 "heat capacity, J/kg/K";
    parameter Modelica.SIunits.Density rhoCustom = 1 "density, kg/m³";

  equation

    if not data.name == "Custom" then
      k_trans =OpenTEMPEST.Blocks.Functions.polySolid(T, data.cK_trans);
      k_long =OpenTEMPEST.Blocks.Functions.polySolid(T, data.cK_long);
      cp =OpenTEMPEST.Blocks.Functions.polySolid(T, data.cCp);
      rho =OpenTEMPEST.Blocks.Functions.polySolid(T, data.cRho);
    else
      k_trans = kCustom_trans;
      k_long = kCustom_long;
      cp = cpCustom;
      rho = rhoCustom;
    end if;

  end BaseProperties;
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end SolidMatBase;

within OpenTEMPEST.SOC.Electrochem.Interfaces;
connector ReactionHeatInfoPort
  Modelica.SIunits.EnthalpyFlowRate H_r;
  Modelica.SIunits.Voltage Vop;
  Modelica.SIunits.Voltage Vid;
  Modelica.SIunits.Current I;
  Modelica.SIunits.Temperature TPEN;
  Modelica.SIunits.EnthalpyFlowRate LHVFlow;

  annotation (
    Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})),
    Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics={  Polygon(origin = {9.74, 0}, fillColor = {230, 230, 230},
            fillPattern =                                                                                                   FillPattern.Solid, points = {{-89.7391, 98}, {-89.7391, -98}, {90.2609, 0}, {90.2609, 0}, {-89.7391, 98}})}));
end ReactionHeatInfoPort;

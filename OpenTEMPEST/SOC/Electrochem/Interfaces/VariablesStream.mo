within OpenTEMPEST.SOC.Electrochem.Interfaces;
connector VariablesStream "Connector for connecting the cell model layers"
  parameter Integer nspecies;
  Modelica.Media.Interfaces.Types.MoleFraction Y[nspecies];
  Modelica.SIunits.AbsolutePressure P;
  Modelica.SIunits.Current I, I_H, I_C;

  annotation (
    Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})),
    Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics={  Polygon(origin={
              7.74,0},                                                                                                                                                 fillColor = {230, 230, 230},
            fillPattern =                                                                                                   FillPattern.Solid, points = {{-89.7391, 98}, {-89.7391, -98}, {90.2609, 0}, {90.2609, 0}, {-89.7391, 98}}), Text(
          extent={{44,-72},{-90,58}},
          lineColor={28,108,200},
          textString="FV")}));
end VariablesStream;

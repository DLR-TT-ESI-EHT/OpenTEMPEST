within OpenTEMPEST.BOP.GasConditioning;
model HeatStorage
  extends Modelica.Icons.UnderConstruction; // Tempest model

  ThermoPower.Gas.FlangeA infl annotation (
    Placement(transformation(extent = {{-120, -20}, {-80, 20}})));
  ThermoPower.Gas.FlangeB outfl annotation (
    Placement(transformation(extent = {{80, -20}, {120, 20}})));
  ThermoPower.Gas.Flow1DFV GasFlow(redeclare model HeatTransfer =
        ThermoPower.Thermal.HeatTransferFV.DittusBoelter)
                                   annotation (
    Placement(transformation(extent = {{-14, -52}, {14, -26}})));
  TEMPEST.Heat.StorageModel storageModel
    annotation (Placement(transformation(extent={{-20,14},{20,-20}})));
equation
  connect(GasFlow.infl, infl) annotation (
    Line(points = {{-14, -39}, {-60, -39}, {-60, 8.88178e-16}, {-100, 8.88178e-16}}, color = {159, 159, 223}, smooth = Smooth.None));
  connect(GasFlow.outfl, outfl) annotation (
    Line(points = {{14, -39}, {60, -39}, {60, 8.88178e-16}, {100, 8.88178e-16}}, color = {159, 159, 223}, smooth = Smooth.None));
  connect(GasFlow.wall, storageModel.internalBoundary) annotation (
    Line(points = {{1.77636e-15, -32.5}, {0, -32.5}, {0, -8.1}, {1.33227e-15, -8.1}}, color = {255, 127, 0}, smooth = Smooth.None));
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics),
    Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent = {{-80, 22}, {80, -20}},
            lineThickness =                                                                                                   1, pattern = LinePattern.None, lineColor = {0, 0, 0},
            fillPattern =                                                                                                   FillPattern.Solid), Ellipse(extent = {{-74, 60}, {-70, -60}}, lineColor = {175, 175, 175},
            lineThickness =                                                                                                   1,
            fillPattern =                                                                                                   FillPattern.Solid, fillColor = {0, 0, 255}), Ellipse(extent = {{70, 60}, {74, -60}}, lineColor = {175, 175, 175},
            lineThickness =                                                                                                   1,
            fillPattern =                                                                                                   FillPattern.Solid, fillColor = {0, 0, 255}), Line(points = {{-72, 60}, {72, 60}}, color = {175, 175, 175}, thickness = 1, smooth = Smooth.None), Line(points = {{-72, -60}, {72, -60}}, color = {175, 175, 175}, thickness = 1, smooth = Smooth.None), Rectangle(extent = {{-70, -20}, {70, -60}}, lineColor = {175, 175, 175},
            lineThickness =                                                                                                   1, fillColor = {0, 0, 255},
            fillPattern =                                                                                                   FillPattern.Solid), Rectangle(extent = {{-70, 60}, {70, 22}}, lineColor = {175, 175, 175},
            lineThickness =                                                                                                   1, fillColor = {0, 0, 255},
            fillPattern =                                                                                                   FillPattern.Solid), Rectangle(extent = {{-70, 24}, {70, -22}},
            lineThickness =                                                                                                   1, fillColor = {0, 0, 255},
            fillPattern =                                                                                                   FillPattern.Solid, pattern = LinePattern.None), Text(extent = {{-44, 22}, {44, -8}}, lineColor = {0, 0, 0},
            lineThickness =                                                                                                   1, fillColor = {0, 0, 255},
            fillPattern =                                                                                                   FillPattern.Solid, textString = "Thermal Storage")}));
end HeatStorage;

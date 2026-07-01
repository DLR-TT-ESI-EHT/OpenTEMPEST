within OpenTEMPEST.Heat;
model RadiativeHeatFlowSplitter
  "Splits a radiative heat flow using only one reference temperature"

  import SI = Modelica.SIunits;

  parameter Integer nCombinedPorts(min=2) "number of heat ports that are combined, the firsts' temperature is used for radiation";
  parameter SI.SpectralEmissivity epsilon1 "emissivity of surface 1";
  parameter SI.SpectralEmissivity epsilonC "emissivity of the combined surface (side2)";
  parameter SI.Area A1 "Area of surface 1";
  parameter SI.Area AC[nCombinedPorts] = {1,1} "Area of combined surfaces";

  replaceable parameter Enumerations.RadiationSpecialCases specialCase=
      Enumerations.RadiationSpecialCases.case3 "Radiation Special Case";

  OpenTEMPEST.Heat.RadiativeHeatTransferGenericEnumerate
    radiativeHeatTransferGeneric(
    epsilon1=epsilon1,
    epsilon2=epsilonC,
    A1=A1,
    A2=sum(AC),
    specialCase=specialCase) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={0,0})));
  ThermoPower.Thermal.HT side1 annotation (
    Placement(transformation(extent={{70,-20},{110,20}}),   iconTransformation(extent={{-8,20},
            {8,36}})));
  ThermoPower.Thermal.HT side2[nCombinedPorts] annotation (Placement(transformation(extent={{
            -110,-20},{-70,20}}), iconTransformation(extent={{-8,-36},{8,-20}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow[
    nCombinedPorts]
    annotation (Placement(transformation(extent={{-46,-10},{-66,10}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
    prescribedTemperature
    annotation (Placement(transformation(extent={{-32,-10},{-12,10}})));
equation

  prescribedTemperature.T = side2[1].T;
  for i in 1:nCombinedPorts loop
    prescribedHeatFlow[i].Q_flow = AC[i]/sum(AC) * prescribedTemperature.port.Q_flow;
  end for;

  connect(radiativeHeatTransferGeneric.side1, side1)
    annotation (Line(points={{2.8,0},{90,0}}, color={191,0,0}));
  connect(prescribedHeatFlow.port, side2)
    annotation (Line(points={{-66,0},{-90,0}}, color={191,0,0}));
  connect(prescribedTemperature.port, radiativeHeatTransferGeneric.side2)
    annotation (Line(points={{-12,0},{-2.8,0}}, color={191,0,0}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={                               Rectangle(extent = {{-80, 20}, {80, -20}}, lineColor = {0, 0, 0}, fillColor = {135, 135, 135},
            fillPattern =                                                                                                   FillPattern.CrossDiag)}),
                                                                 Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end RadiativeHeatFlowSplitter;

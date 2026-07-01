within OpenTEMPEST.Heat;
model RadHT2DFV
  "Distributed Radiation block model - uses multiple single blocks"

  parameter Integer nX = 5 "Number of volumes in x-direction";
  parameter Integer nY = 5 "Number of volumes in y-direction";
  replaceable parameter Enumerations.RadiationSpecialCases specialCase=
      Enumerations.RadiationSpecialCases.case3;
  parameter Real epsilon1;
  parameter Real epsilon2;
  parameter Real A1;
  parameter Real A2;

  OpenTEMPEST.Heat.RadiativeHeatTransferGenericEnumerate
    radiativeHeatTransferGenericEnumerate[nX,nY](
    each epsilon1=epsilon1,
    each epsilon2=epsilon2,
    each A1=A1,
    each A2=A2,
    each specialCase=specialCase)
    annotation (Placement(transformation(extent={{-60,6},{-40,26}})));
  OpenTEMPEST.Heat.HTs_DHTVolumes2D hTs_DHTVolumes2D(i=nX, j=nY) annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-50,46})));
  OpenTEMPEST.Heat.HTs_DHTVolumes2D hTs_DHTVolumes2D1(i=nX, j=nY) annotation (
      Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={-50,-20})));
  OpenTEMPEST.Heat.DHTVolumes2D side1(i=nX, j=nY) annotation (Placement(
        transformation(extent={{-20,40},{0,48}}), iconTransformation(extent={{-20,
            40},{20,48}})));
  OpenTEMPEST.Heat.DHTVolumes2D side2(i=nX, j=nY) annotation (Placement(
        transformation(extent={{-20,-8},{0,0}}), iconTransformation(extent={{-20,
            -8},{20,0}})));
equation
  connect(hTs_DHTVolumes2D.HT_ports, radiativeHeatTransferGenericEnumerate.side1)
    annotation (Line(points={{-50,34},{-50,18.8}}, color={191,0,0}));
  connect(hTs_DHTVolumes2D.DHT_port, side1) annotation (Line(points={{-50,57},{-52,
          57},{-52,76},{-10,76},{-10,44}}, color={0,0,0}));
  connect(hTs_DHTVolumes2D1.HT_ports, radiativeHeatTransferGenericEnumerate.side2)
    annotation (Line(points={{-50,-8},{-50,13.2}}, color={191,0,0}));
  connect(hTs_DHTVolumes2D1.DHT_port, side2) annotation (Line(points={{-50,-31},
          {-50,-40},{-10,-40},{-10,-4}}, color={0,0,0}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
          extent={{-100,40},{100,0}},
          lineColor={0,0,0},
          fillColor={95,95,95},
          fillPattern=FillPattern.CrossDiag)}),                  Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end RadHT2DFV;

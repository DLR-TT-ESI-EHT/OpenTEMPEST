within OpenTEMPEST.Heat;
model RadHTFV "Distributed Radiation block model - uses multiple single blocks"

  parameter Integer N(min=1) "Number of discretisation Units";
  replaceable parameter Enumerations.RadiationSpecialCases specialCase=
      Enumerations.RadiationSpecialCases.case3;
  parameter Real epsilon1;
  parameter Real epsilon2;
  parameter Real A1;
  parameter Real A2;

  OpenTEMPEST.Heat.RadiativeHeatTransferGenericEnumerate
    radiativeHeatTransferGenericEnumerate[N](
    each epsilon1=epsilon1,
    each epsilon2=epsilon2,
    each A1=A1,
    each A2=A2,
    each specialCase=specialCase)
    annotation (Placement(transformation(extent={{-60,6},{-40,26}})));
  OpenTEMPEST.Heat.HTs_DHT hTs_DHT(N=N) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=-90,
        origin={-50,42})));
  ThermoPower.Thermal.DHTVolumes side1(N=N) annotation (Placement(
        transformation(extent={{-20,40},{0,48}}), iconTransformation(extent={
            {-20,40},{20,48}})));
  ThermoPower.Thermal.DHTVolumes side2(N=N) annotation (Placement(
        transformation(extent={{-20,-8},{0,0}}), iconTransformation(extent={{
            -20,-8},{20,0}})));
  OpenTEMPEST.Heat.HTs_DHT hTs_DHT1(N=N) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={-50,-6})));
equation
  connect(hTs_DHT.HT_ports, radiativeHeatTransferGenericEnumerate.side1)
    annotation (Line(points={{-49.8,30},{-49.8,26},{-50,26},{-50,18.8}},
                                                                       color={191,
          0,0}));
  connect(hTs_DHT.DHT_port, side1) annotation (Line(points={{-50,53},{-50,69.5},
          {-10,69.5},{-10,44}},       color={255,127,0}));
  connect(radiativeHeatTransferGenericEnumerate.side2, hTs_DHT1.HT_ports)
    annotation (Line(points={{-50,13.2},{-50,6},{-49.8,6}},        color={191,0,
          0}));
  connect(hTs_DHT1.DHT_port, side2) annotation (Line(points={{-50,-17},{-50,
          -34.5},{-10,-34.5},{-10,-4}}, color={255,127,0}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
          extent={{-100,40},{100,0}},
          lineColor={0,0,0},
          fillColor={95,95,95},
          fillPattern=FillPattern.CrossDiag)}),                  Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end RadHTFV;

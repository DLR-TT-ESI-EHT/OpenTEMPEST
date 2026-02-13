within OpenTEMPEST.Heat;
model DHT_2DHT
  ThermoPower.Thermal.DHTVolumes DHTSingle(N=n)
    annotation (Placement(transformation(extent={{-10,60},{10,80}}),
        iconTransformation(extent={{-60,40},{60,60}})));

  parameter Integer n "number of volumens at single port";
  parameter Real r "split ratio of the volumes";

  ThermoPower.Thermal.DHTVolumes DHT1(N=integer(r*n))
    annotation (Placement(transformation(extent={{-80,-80},{-60,-60}}),
        iconTransformation(extent={{-80,-60},{-20,-40}})));
  ThermoPower.Thermal.DHTVolumes DHT2(N=n - DHT1.N)
    annotation (Placement(transformation(extent={{60,-60},{80,-40}}),
        iconTransformation(extent={{20,-60},{80,-40}})));
protected
  OpenTEMPEST.Heat.HTs_DHT hTs_DHT(N=n) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={0,30})));
  OpenTEMPEST.Heat.HTs_DHT hTs_DHT1(N=DHT1.N) annotation (Placement(
        transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={-40,-30})));
  OpenTEMPEST.Heat.HTs_DHT hTs_DHT2(N=DHT2.N) annotation (Placement(
        transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={40,-30})));
equation
  connect(hTs_DHT.DHT_port, DHTSingle) annotation (Line(points={{8.88178e-16,41},
          {0,41},{0,70}}, color={255,127,0}));
  connect(hTs_DHT1.DHT_port, DHT1) annotation (Line(points={{-40,-41},{-56,-41},
          {-56,-70},{-70,-70}}, color={255,127,0}));
  connect(hTs_DHT2.DHT_port, DHT2) annotation (Line(points={{40,-41},{60,-41},{60,
          -50},{70,-50}}, color={255,127,0}));
  connect(hTs_DHT.HT_ports[1:DHT1.N], hTs_DHT1.HT_ports[:]) annotation (Line(points={{-0.2,
          18},{-20,18},{-20,-18},{-40.2,-18}}, color={191,0,0}));
  connect(hTs_DHT.HT_ports[DHT1.N+1:n], hTs_DHT2.HT_ports[:]) annotation (Line(points={{-0.2,
          18},{17.9,18},{17.9,-18},{39.8,-18}}, color={191,0,0}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Line(
          points={{-30,40},{-50,-40}},
          color={255,127,0},
          thickness=1),
        Line(
          points={{30,40},{50,-40}},
          color={255,127,0},
          thickness=1),
        Text(
          textString="",
          extent={{-38,56},{-40,54}},
          lineColor={255,127,0},
          lineThickness=1),
        Text(
          extent={{-20,34},{22,26}},
          lineColor={0,0,0},
          lineThickness=1,
          textString="n = %n"),
        Text(
          extent={{-78,-22},{-20,-40}},
          lineColor={0,0,0},
          lineThickness=1,
          textString="n = %DHT1.N")}),                           Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end DHT_2DHT;

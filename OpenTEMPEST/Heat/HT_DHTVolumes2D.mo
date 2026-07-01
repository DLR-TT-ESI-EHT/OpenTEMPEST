within OpenTEMPEST.Heat;
model HT_DHTVolumes2D "HT to DHT adaptor"
  parameter Integer i=1 "number of volumes in direction 1";
  parameter Integer j=1 "number of volumes in direction 2";
  ThermoPower.Thermal.HT HT_port annotation (Placement(transformation(extent={{-140,-20},{-100,20}}, rotation=0),
        iconTransformation(extent={{-140,-20},{-100,20}})));
  OpenTEMPEST.Heat.DHTVolumes2D DHT_port(i=i, j=j) annotation (Placement(
        transformation(extent={{100,-40},{120,40}}, rotation=0)));
equation
  for k in 1:i loop
    for n in 1:j loop
      DHT_port.T[k,n] = HT_port.T "Uniform temperature distribution on DHT side";
    end for;
  end for;
  sum(DHT_port.Q) + HT_port.Q_flow = 0 "Energy balance";
  annotation (Icon(graphics={
        Polygon(
          points={{-100,100},{-100,-100},{100,100},{-100,100}},
          lineColor={185,0,0},
          fillColor={185,0,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{100,100},{100,-100},{-100,-100},{100,100}},
          lineColor={0,0,0},
          fillColor={255,240,0},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-84,20},{2,82}},
          lineColor={255,255,255},
          lineThickness=1,
          textString="HT"),
        Text(
          extent={{-30,-128},{90,0}},
          lineColor={0,0,0},
          lineThickness=1,
          textString="2D DHT"),
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          pattern=LinePattern.None)}), Diagram(graphics));
end HT_DHTVolumes2D;

within OpenTEMPEST.Heat;
model HTs_DHTVolumes2D "HT to DHT adaptor"
  parameter Integer i=1 "number of volumes in direction 1";
  parameter Integer j=1 "number of volumes in direction 2";
  ThermoPower.Thermal.HT HT_ports[i,j] annotation (Placement(transformation(extent={{-140,-20},{-100,20}}, rotation=0),
        iconTransformation(extent={{-140,-20},{-100,20}})));
  OpenTEMPEST.Heat.DHTVolumes2D DHT_port(i=i, j=j) annotation (Placement(
        transformation(extent={{100,-40},{120,40}}, rotation=0)));
equation
  for k in 1:i loop
    for n in 1:j loop
      assert(cardinality(HT_ports[k, n]) <= 1, "
     each HT_ports[i] of boundary shall at most be connected to one component.
     If two or more connections are present, increase n to add an additional port.");
      DHT_port.T[k,n] = HT_ports[k,n].T
        "Uniform temperature distribution on DHT side";
      DHT_port.Q[k,n] + HT_ports[k,n].Q_flow = 0 "Energy balance";
    end for;
  end for;

  annotation (Icon(graphics={
        Polygon(
          points={{100,101},{100,-99},{-100,-99},{100,101}},
          lineColor={0,0,0},
          fillColor={255,240,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-100,100},{-100,-100},{100,100},{-100,100}},
          lineColor={185,0,0},
          fillColor={185,0,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          pattern=LinePattern.None),
        Text(
          extent={{-30,-128},{90,0}},
          lineColor={0,0,0},
          lineThickness=1,
          textString="2D DHT"),
        Text(
          extent={{-84,12},{4,86}},
          lineColor={255,255,255},
          lineThickness=1,
          textString="HTs")}),         Diagram(graphics));
end HTs_DHTVolumes2D;

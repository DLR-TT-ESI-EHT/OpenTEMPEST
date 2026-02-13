within OpenTEMPEST.Heat;
model HTs_DHT "Several HTs to DHT adapterFV"
  import SI = Modelica.SIunits;
  parameter Integer N "Number of HTs (= number of nodes on DHT side)";

  ThermoPower.Thermal.HT
     HT_ports[N] annotation (Placement(transformation(extent={{-140,-16},{-100,24}},
          rotation=0), iconTransformation(extent={{-140,-18},{-100,22}})));
  ThermoPower.Thermal.DHTVolumes
      DHT_port(N=N) annotation (Placement(transformation(extent={{100,-40},{
            120,40}}, rotation=0), iconTransformation(extent={{100,-40},{120,40}})));
equation
  for i in 1:N loop
     assert(cardinality(HT_ports[i]) <= 1, "
     each HT_ports[i] of boundary shall at most be connected to one component.
     If two or more connections are present, increase n to add an additional port.");
    DHT_port.T[i] = HT_ports[i].T
      "Uniform temperature distribution on DHT side";
    DHT_port.Q[i] + HT_ports[i].Q_flow = 0
      "Energy balance";
  end for;
  annotation (Icon(graphics={
        Polygon(
          points={{-100,98},{-100,-102},{100,98},{-100,98}},
          lineColor={185,0,0},
          fillColor={185,0,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{100,98},{100,-102},{-100,-102},{100,98}},
          lineColor={255,128,0},
          fillColor={255,128,0},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-74,8},{24,86}},
          lineColor={255,255,255},
          lineThickness=1,
          textString="HTs"),
        Text(
          extent={{-16,-86},{82,-8}},
          lineColor={255,255,255},
          lineThickness=1,
          textString="DHT"),
        Rectangle(
          extent={{-100,98},{100,-102}},
          lineColor={0,0,0},
          pattern=LinePattern.None)}), Diagram(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}})));
end HTs_DHT;

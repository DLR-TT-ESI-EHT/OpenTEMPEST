within OpenTEMPEST.Heat;
model HTs_HT "Several HTs to one HT adapter"
  import SI = Modelica.SIunits;
  parameter Integer N(min=1) "Number of HTs (= number of nodes on DHT side)";

  ThermoPower.Thermal.HT
     HT_ports[N] annotation (Placement(transformation(extent={{-140,-20},{-100,20}},
          rotation=0), iconTransformation(extent={{-140,-20},{-100,20}})));
  ThermoPower.Thermal.HT HT_port annotation (Placement(transformation(extent={{100,
            -20},{140,20}}, rotation=0), iconTransformation(extent={{100,-22},{140,
            18}})));
equation
  for i in 1:N loop
    HT_port.T = HT_ports[i].T
      "Uniform temperature distribution on DHT side";
  end for;
  HT_port.Q_flow + sum(HT_ports[:].Q_flow) = 0 "Energy balance";
    annotation (Icon(graphics={
        Polygon(
          points={{-100,100},{-100,-100},{100,100},{-100,100}},
          lineColor={185,0,0},
          fillColor={185,0,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{100,98},{100,-102},{-100,-102},{100,98}},
          lineColor={185,0,0},
          fillColor={185,0,0},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-74,10},{24,88}},
          lineColor={255,255,255},
          lineThickness=1,
          textString="HTs"),
        Text(
          extent={{-16,-86},{82,-8}},
          lineColor={255,255,255},
          lineThickness=1,
          textString="HT"),
        Rectangle(
          extent={{-100,98},{100,-102}},
          lineColor={0,0,0},
          pattern=LinePattern.None)}), Diagram(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}})));
end HTs_HT;

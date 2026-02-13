within OpenTEMPEST.Heat;
model DHT_DHTVolumes2D "DHT to DHTVolumes2D adaptor"
  parameter Integer i=1 "number of DHT volumes (direction 1)";
  parameter Integer j=1 "number elements in DHT Volume (direction 2)";
  ThermoPower.Thermal.DHTVolumes DHT_port[i](each N=j) annotation (Placement(transformation(extent={{-140,-20},{-100,20}}, rotation=0),iconTransformation(extent={{-140,-20},{-100,20}})));
  OpenTEMPEST.Heat.DHTVolumes2D DHT2D_port(i=i, j=j) annotation (Placement(
        transformation(extent={{100,-40},{120,40}}, rotation=0)));
equation
  for k in 1:i loop
    DHT2D_port.T[k,:] = DHT_port[k].T[:] "Uniform temperature distribution on DHT side";
    for jj in 1:j loop
      DHT2D_port.Q[k,jj] + DHT_port[k].Q[jj] = 0 "Energy balance";
    end for;
  end for;

  annotation (Icon(graphics={
        Polygon(
          points={{-100,100},{-100,-100},{100,100},{-100,100}},
          lineColor={0,0,0},
          fillColor={255,127,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{100,100},{100,-100},{-100,-100},{100,100}},
          lineColor={255,128,0},
          fillColor={231,231,0},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-74,10},{24,88}},
          lineColor={255,255,255},
          lineThickness=1,
          fillColor={255,127,36},
          fillPattern=FillPattern.Solid,
          textString="DHT"),
        Text(
          extent={{-16,-82},{82,-4}},
          lineColor={255,255,255},
          lineThickness=1,
          textString="2D DHT",
          fillColor={238,238,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          pattern=LinePattern.None)}));
end DHT_DHTVolumes2D;

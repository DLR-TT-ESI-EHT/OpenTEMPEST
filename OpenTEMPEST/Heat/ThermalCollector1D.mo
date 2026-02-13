within OpenTEMPEST.Heat;
model ThermalCollector1D "Collects m heat flows through 1D thermal ports"
  parameter Integer m(min=1)=3 "Number of collected heat flows";
  parameter Integer N "Number of control volumes";
  ThermoPower.Thermal.DHTVolumes dHT0[m](each N=N)
                                           annotation (Placement(transformation(extent={{-6,106},{6,94}}),
        iconTransformation(extent={{-6,106},{6,94}})));
  ThermoPower.Thermal.DHTVolumes dHT1(each N=N)
                                        annotation (Placement(transformation(extent={{-6,-106},{6,-94}}),
        iconTransformation(extent={{-6,-106},{6,-94}})));

equation
    dHT1.Q .+ sum(dHT0.Q) = zeros(N);
    dHT0.T = fill(dHT1.T, m);

  annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
            100,100}}), graphics={
        Text(
          extent={{-150,-30},{150,-70}},
          textString="%name",
          lineColor={0,0,255}),
        Text(
          extent={{-150,80},{150,50}},
          textString="m=%m"),
        Line(
          points={{0,90},{0,40}},
          color={244,125,35}),
        Rectangle(
          extent={{-60,40},{60,30}},
          lineColor={181,0,0},
          fillColor={244,125,35},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-60,30},{0,-30},{0,-90}},
          color={244,125,35}),
        Line(
          points={{0,-30},{-20,30}},
          color={244,125,35}),
        Line(
          points={{0,-30},{20,30}},
          color={244,125,35}),
        Line(
          points={{0,-30},{60,30}},
          color={244,125,35})}),
    Documentation(info="<html>
<p>
This is a model to collect the heat flows from <em>m</em> heatports to one single heatport.
</p>
</html>"));
end ThermalCollector1D;

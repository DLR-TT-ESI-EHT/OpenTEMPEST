within OpenTEMPEST.Blocks.Sources;
block VectorRamp

  parameter Real height[nout] = target-offset "Height of ramps";
  parameter Real target[nout] = fill(1,nout) "Target of the output signals";
  parameter Modelica.SIunits.Time duration(min=0.0, start=2)
    "Duration of ramp (= 0.0 gives a Step)";
  parameter Real offset[nout] = fill(1,nout) "Offset of output signals";
  parameter Modelica.SIunits.Time startTime=0
    "Output = offset for time < startTime";

  extends Modelica.Blocks.Interfaces.MO;
  Modelica.Blocks.Sources.Ramp ramp[nout](
  height=height,
  each duration=duration,
  offset=offset,
  each startTime=startTime)
  annotation (Placement(transformation(extent={{-20,-10},{0,10}})));
equation
connect(ramp.y, y) annotation (Line(points={{1,0},{110,0}}, color={0,0,127}));
annotation (Icon(graphics={     Rectangle(
        extent={{-100,-100},{100,100}},
        lineColor={0,0,127},
        fillColor={255,255,255},
        fillPattern=FillPattern.Solid),
        Line(points={{-80,68},{-80,-80}}, color={192,192,192}),
        Polygon(
          points={{-80,90},{-88,68},{-72,68},{-80,90}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(points={{-90,-70},{82,-70}}, color={192,192,192}),
        Polygon(
          points={{90,-70},{68,-62},{68,-78},{90,-70}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(points={{-80,-70},{-40,-70},{31,38}}),
        Line(points={{31,38},{86,38}})}), Documentation(revisions="<html>
<ul>
<li>&quot;before versioning&quot; by Marius Tomberg: <br>Created</li>
<li>14-07-2021 by Marius Tomberg:<br>Target can be set as parameters instead of ramp</li>
</ul>
</html>"));
end VectorRamp;

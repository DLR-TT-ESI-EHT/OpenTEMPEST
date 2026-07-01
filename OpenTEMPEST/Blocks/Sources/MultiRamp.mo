within OpenTEMPEST.Blocks.Sources;
model MultiRamp "Generate signal with n ramps"
  parameter Integer n=3 "number of ramps";
  parameter Real height[n]=ones(n) "Height of ramps";
  parameter Modelica.SIunits.Time duration[n](min=zeros(n), start=2*ones(n))=zeros(n)
    "Duration of ramp (= 0.0 gives a Step)";
  parameter Real offset=0 "Offset of output signal";
  parameter Modelica.SIunits.Time startTime[n]=zeros(n)
    "Output = offset for time < startTime";
  extends Modelica.Blocks.Interfaces.SO;
  Modelica.Blocks.Sources.Ramp ramp[n](
    height=height,
    duration=duration,
    each offset=0,
    startTime=startTime)
    annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
  Modelica.Blocks.Math.MultiSum multiSum(nu=n+1)
    annotation (Placement(transformation(extent={{10,-6},{22,6}})));
  Modelica.Blocks.Sources.Constant constOffset(k=offset)
    annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
equation
  connect(ramp.y, multiSum.u[2:n+1]) annotation (Line(points={{-59,0},{-24,0},{-24,
          0},{10,0}},       color={0,0,127}));
  connect(multiSum.y, y)
    annotation (Line(points={{23.02,0},{110,0}}, color={0,0,127}));
  connect(constOffset.y, multiSum.u[1]) annotation (Line(points={{-59,50},{-24.5,
          50},{-24.5,0},{10,0}}, color={0,0,127}));
  annotation (
    Icon(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}}), graphics={
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
        Line(points={{-80,-70},{-40,-70},{-20,-20}}),
        Text(
          extent={{-150,-150},{150,-110}},
          lineColor={0,0,0},
          textString="duration=%duration"),
        Line(points={{30,40},{80,40}}),
        Line(points={{-20,-20},{25,-20}}),
        Line(points={{25,-20},{30,40}})}),
    Diagram(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}})),
    Documentation(info="<html>
<p>
The Real output y is a ramp signal:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Blocks/Sources/Ramp.png\"
     alt=\"Ramp.png\">
</p>

<p>
If parameter duration is set to 0.0, the limiting case of a Step signal is achieved.
</p>
</html>", revisions="<html>
<ul>
<li>&quot;before versioning&quot; by Marius Tomberg: <br>Created</li>
</ul>
</html>"));
end MultiRamp;

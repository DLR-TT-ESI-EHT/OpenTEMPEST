within OpenTEMPEST.Blocks.Sources;
model Pulse
  "Generate pulses - faster than MSL pulse source. Uses MultiRamp in the background"

  parameter Integer nPulse=500 "Number of Pulses";
  parameter Real pulseHeight "Height of Pulse";
  parameter Modelica.SIunits.Time startPulse = 1 "";
  parameter Modelica.SIunits.Time pulseDuration = 1 "Time for one period (onTime + offtime)";
  parameter Real dutyCycle(min=0, max=1) = 0.5 "Width of Pulse (0 to 1)";
  parameter Real offset=0 "Offset of output signal";

protected
  parameter Integer b[2*nPulse]=0:2*nPulse-1;
  parameter Real factor[2*nPulse] = {(-1) .^i for i in b};

  parameter Real A[2*nPulse] = pulseHeight .*factor;

  parameter Real onTime = dutyCycle*pulseDuration;
  parameter Real offTime = (1-dutyCycle)*pulseDuration;

  parameter Real onFactor[2*nPulse] = {mod(i,2) for i in b};

  parameter Real switchTime[2*nPulse] = onTime.*onFactor .+ cat(1, {0},offTime.*(-onFactor[2:end].+1));

  parameter Real pulseTime[2*nPulse] = startPulse*ones(2*nPulse) + {sum(switchTime[1:i]) for i in 1:2*nPulse};

  parameter Integer n=2*nPulse "number of ramps";
  parameter Real height[n]=A "Height of ramps";
  parameter Modelica.SIunits.Time duration[n](min=zeros(n), start=2*ones(n))=zeros(n)
    "Duration of ramp (= 0.0 gives a Step)";

  parameter Modelica.SIunits.Time startTime[n]=pulseTime
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
        Line(points={{-80,-70},{-40,-70},{-40,44},{0,44},{0,-70},{40,-70},{40,
              44},{79,44}}),
        Text(
          extent={{-147,-152},{153,-112}},
          textString="period=%period")}),
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
end Pulse;

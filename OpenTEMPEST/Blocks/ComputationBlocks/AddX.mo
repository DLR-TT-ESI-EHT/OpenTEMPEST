within OpenTEMPEST.Blocks.ComputationBlocks;
block AddX "Output the sum of the X inputs"
  extends Modelica.Blocks.Icons.Block;

  parameter Integer n = 3 "Number of input signals";

  parameter Real k[n]=ones(n) "Gain of input signals";

  Modelica.Blocks.Interfaces.RealInput u[n] "Connectors of Real input signals"
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
  Modelica.Blocks.Interfaces.RealOutput y "Connector of Real output signal"
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));

equation
  y = sum(k*u);
  annotation (
    Documentation(info="<html>
<p>
This blocks computes output <strong>y</strong> as <em>sum</em> of the
three input signals <strong>u1</strong>, <strong>u2</strong> and <strong>u3</strong>:
</p>
<pre>
    <strong>y</strong> = k1*<strong>u1</strong> + k2*<strong>u2</strong> + k3*<strong>u3</strong>;
</pre>
<p>
Example:
</p>
<pre>
     parameter:   k1= +2, k2= -3, k3=1;

  results in the following equations:

     y = 2 * u1 - 3 * u2 + u3;
</pre>

</html>", revisions="<html>
<ul>
<li><i>14-07-2021</i> by Marius Tomberg: <br>First release</li>
</ul>
</html>"),
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
            100}}), graphics={
        Text(
          extent={{-82,-22},{3,18}},
          lineColor={0,0,0},
          textString="%k"),
        Text(
          extent={{2,36},{100,-44}},
          textString="+")}),
    Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
            100,100}}), graphics={Rectangle(
            extent={{-100,-100},{100,100}},
            lineColor={0,0,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
                             Text(
            extent={{-100,-20},{5,20}},
          textString="kX",
          lineColor={0,0,0}),Text(
            extent={{2,46},{100,-34}},
            textString="+")}));
end AddX;

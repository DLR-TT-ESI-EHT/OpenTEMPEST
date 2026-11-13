within OpenTEMPEST.Blocks.ComputationBlocks;
block AddX "Weighted sum of n input signals"
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
<p>Computes the output <code>y</code> as the weighted sum of the <code>n</code> inputs:</p>
<pre>    y = k[1]*u[1] + k[2]*u[2] + ... + k[n]*u[n]</pre>
</html>", revisions="<html>
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

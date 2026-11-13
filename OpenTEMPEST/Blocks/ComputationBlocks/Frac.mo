within OpenTEMPEST.Blocks.ComputationBlocks;
block Frac "Product of n input signals, each raised to a given exponent"
  extends Modelica.Blocks.Icons.Block;

  parameter Integer n=3 "Number of input signals";
  parameter Real k[n]={1,-1,1} "Exponent of input signals";
  Modelica.Blocks.Interfaces.RealInput u[n] "Connector of Real input signals"
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}}),
        iconTransformation(extent={{-140,-20},{-100,20}})));
  Modelica.Blocks.Interfaces.RealOutput y "Connector of Real output signal"
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));

equation
  y = product(u[i]^k[i] for i in 1:n);
  annotation (
    Documentation(info="<html>
<p>Computes the output <code>y</code> as the product of the <code>n</code> inputs, each raised to the exponent <code>k[i]</code>:</p>
<pre>    y = u[1]^k[1] * u[2]^k[2] * ... * u[n]^k[n]</pre>
<p>With the default <code>k = {1, -1, 1}</code> this gives <code>y = u[1]*u[3]/u[2]</code>.</p>
</html>"),
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
            100}})),
    Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
            100,100}}), graphics={Rectangle(
            extent={{-100,-100},{100,100}},
            lineColor={0,0,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}));
end Frac;

within OpenTEMPEST.Blocks.ComputationBlocks;
block Poly "Output a polynomical expression, y = c[1] + u*(c[2] + u*(c[3] + u*(c[4] + ...)))"
  extends Modelica.Blocks.Interfaces.SISO;

  parameter Real  c[:] "Polynomial coefficients";

  Modelica.Blocks.Interfaces.RealInput u "Connector of Real input signals"
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}}),
        iconTransformation(extent={{-140,-20},{-100,20}})));
  Modelica.Blocks.Interfaces.RealOutput y "Connector of Real output signal"
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));

equation
  y = Modelica.Math.Special.Internal.polyEval(c,u);
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

</html>"),
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
            100}}), graphics={Text(
          extent={{-96,96},{96,-96}},
          lineColor={162,162,162},
          textString="y = c[1] + u*(c[2] + u*(c[3] + ...))")}),
    Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
            100,100}}), graphics={Rectangle(
            extent={{-100,-100},{100,100}},
            lineColor={0,0,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}));
end Poly;

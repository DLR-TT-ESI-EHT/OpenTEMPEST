within OpenTEMPEST.Blocks.ComputationBlocks;
block Poly "Polynomial of the input, y = c[1] + u*(c[2] + u*(c[3] + ...))"
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
<p>Computes the output <code>y</code> as a polynomial of the input <code>u</code> with coefficients <code>c</code>:</p>
<pre>    y = c[1] + c[2]*u + c[3]*u^2 + ...</pre>
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

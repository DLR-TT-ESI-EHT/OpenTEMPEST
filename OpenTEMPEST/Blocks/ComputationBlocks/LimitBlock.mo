within OpenTEMPEST.Blocks.ComputationBlocks;
block LimitBlock "Limit a signal between a minimum and a maximum value"
  extends Modelica.Blocks.Icons.Block;

  parameter Real minVal;
  parameter Real maxVal;

  Modelica.Blocks.Interfaces.RealInput u "Connector of Real input signals"
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}}),
        iconTransformation(extent={{-140,-20},{-100,20}})));
  Modelica.Blocks.Interfaces.RealOutput y "Connector of Real output signal"
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));

equation
  y = min(maxVal,max(minVal,u));

algorithm

  annotation (
    Documentation(info="<html>
<p>Limits the input signal to the range <code>[minVal, maxVal]</code>:</p>
<pre>    y = min(maxVal, max(minVal, u))</pre>
</html>"),
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
            100}})),
    Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
            100,100}}), graphics={Rectangle(
            extent={{-100,-100},{100,100}},
            lineColor={0,0,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}));
end LimitBlock;

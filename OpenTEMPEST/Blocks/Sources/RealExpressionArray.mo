within OpenTEMPEST.Blocks.Sources;
block RealExpressionArray "Set output signal to a time varying Real expression"
  parameter Integer nAxes = 4 "number of axes of output array y";
  parameter Integer arrayOut[nAxes] = {3,4,30,5} "shape of output array y";
  Modelica.Blocks.Interfaces.RealOutput y[arrayOut[1],arrayOut[2],arrayOut[3],arrayOut[4]]= fill(fill(fill(fill(0,5), 30), 4),3) "Values of Real output array y"
    annotation (Dialog(group="Time varying output signal"), Placement(
        transformation(extent={{100,-10},{120,10}})));

  annotation (Icon(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}}), graphics={
        Rectangle(
          extent={{-100,40},{100,-40}},
          fillColor={235,235,235},
          fillPattern=FillPattern.Solid,
          borderPattern=BorderPattern.Raised),
        Text(
          extent={{-96,15},{96,-15}},
          textString="%y"),
        Text(
          extent={{-150,90},{150,50}},
          textString="%name",
          lineColor={0,0,255})}), Documentation(info="<html>
<p>
The (time varying) Real output signal of this block can be defined in its
parameter menu via variable <strong>y</strong>. The purpose is to support the
easy definition of Real expressions in a block diagram. For example,
in the y-menu the definition \"if time &lt; 1 then 0 else 1\" can be given in order
to define that the output signal is one, if time &ge; 1 and otherwise
it is zero. Note, that \"time\" is a built-in variable that is always
accessible and represents the \"model time\" and that
variable <strong>y</strong> is both a variable and a connector.
</p>
</html>"));

end RealExpressionArray;

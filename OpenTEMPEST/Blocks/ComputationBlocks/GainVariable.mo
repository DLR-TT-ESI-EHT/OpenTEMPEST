within OpenTEMPEST.Blocks.ComputationBlocks;
block GainVariable "Output the product of a gain value with the input signal"

  parameter Real k(start=1,unit="1")=0.7  "Gain value multiplied with input signal";

public
  Modelica.Blocks.Interfaces.RealInput u "Input signal connector"
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
  Modelica.Blocks.Interfaces.RealOutput y "Output signal connector" annotation (
     Placement(transformation(extent={{100,-12},{126,14}}), iconTransformation(
          extent={{100,-20},{140,20}})));
  parameter Boolean use_in_k = true "Use connector input for the composition" annotation(Dialog(group="External inputs"), choices(checkBox=true));
protected
  Real k_internal;
public
  Modelica.Blocks.Interfaces.RealInput in_k if use_in_k "Input signal connector"  annotation (Placement(transformation(
        extent={{23,-23},{-23,23}},
        rotation=90,
        origin={-37,57}), iconTransformation(
        extent={{23,-23},{-23,23}},
        rotation=90,
        origin={-18,50})));
equation

  if use_in_k ==true then
    k_internal= in_k;
  else
    k_internal=k;
  end if;

  y = k_internal*u;

  annotation (
    Documentation(info="<html>
<p>
This block computes output <em>y</em> as
<em>product</em> of gain <em>k</em> with the
input <em>u</em>:
</p>
<pre>
    y = k * u;
</pre>

</html>"),
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
            100}}), graphics={
        Polygon(
          points={{-100,-100},{-100,100},{100,0},{-100,-100}},
          lineColor={0,0,127},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-150,-140},{150,-100}},
          textString="k=%k"),
        Text(
          extent={{-150,140},{150,100}},
          textString="%name",
          lineColor={0,0,255}),
        Text(
          extent={{-40,32},{6,-20}},
          lineColor={28,108,200},
          fillColor={0,0,0},
          fillPattern=FillPattern.None,
          textString="k")}),
    Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
            100,100}}), graphics={Polygon(
            points={{-100,-100},{-100,100},{100,0},{-100,-100}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),Text(
            extent={{-76,38},{0,-34}},
            textString="k",
            lineColor={0,0,255})}));
end GainVariable;

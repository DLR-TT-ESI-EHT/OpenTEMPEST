within OpenTEMPEST.Blocks.ComputationBlocks;
model BumplessTransferController
  "Changes the setpoint of an external controller to reduce controller error peaks and toggles a switch between automatic and manual tracking for the reference variable"

//   parameter Boolean use_in_setPoint2 = false "Select to activate manual tracking for controller"  annotation(Dialog(group="External inputs"), choices(checkBox=true));
  parameter Modelica.SIunits.Time rampCompensatorDuration = 180 "Ramp Compensator duration";
  parameter Modelica.SIunits.Time rampStartTime = 64800 "Ramp Compensator start time";

  Modelica.Blocks.Interfaces.RealOutput out_setPointVariable
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));
  Modelica.Blocks.Interfaces.RealInput in_setPoint2
    annotation (Placement(transformation(
        extent={{-13,-13},{13,13}},
        rotation=270,
        origin={20,112}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=270,
        origin={40,120})));// if use_in_setPoint2
  Modelica.Blocks.Interfaces.RealInput in_setPoint1 annotation (Placement(
        transformation(
        extent={{-13,-13},{13,13}},
        rotation=0,
        origin={-110,20}), iconTransformation(extent={{-140,-20},{-100,20}})));
  Modelica.Blocks.Math.Product rampMultiplier annotation (Placement(
        transformation(
        extent={{-5,-5},{5,5}},
        rotation=0,
        origin={-20,20})));
  Modelica.Blocks.Math.Add add(k2=-1)
    annotation (Placement(transformation(extent={{-48,6},{-36,18}})));
  Modelica.Blocks.Math.Add add1
    annotation (Placement(transformation(extent={{8,-6},{20,6}})));
  Modelica.Blocks.Interfaces.BooleanInput in_setPoint2Active  "Trigger for activation of manual tracking" annotation (
      Placement(transformation(
        extent={{-14,-14},{14,14}},
        rotation=270,
        origin={-20,110}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=270,
        origin={-40,120}))); // if use_in_setPoint2
  Modelica.Blocks.Logical.TriggeredTrapezoid triggeredTrapezoid(amplitude=1,
      rising=rampCompensatorDuration)
    annotation (Placement(transformation(extent={{-64,32},{-44,52}})));
  Modelica.Blocks.MathBoolean.Not not1 annotation (Placement(transformation(extent={{-84,38},{-76,46}})));
//   Modelica.Blocks.Interfaces.RealInput in_setPoint2_internal;
//   Modelica.Blocks.Interfaces.BooleanInput in_setPoint2Active_internal;
equation

//    if not use_in_setPoint2 then
//     in_setPoint2_internal = automaticControlDummy.k;
//     in_setPoint2Active_internal = false;
//    end if;

  // Connect internal connector
//   connect(in_setPoint2_internal, in_setPoint2);
//   connect(in_setPoint2_internal, add.u2);
//   connect(in_setPoint2_internal, add1.u2);
//   connect(in_setPoint2Active_internal, in_setPoint2Active);

  connect(rampMultiplier.y, add1.u1) annotation (Line(points={{-14.5,20},{6.8,20},
          {6.8,3.6}},              color={0,0,127}));
  connect(in_setPoint1, add.u1) annotation (Line(points={{-110,20},{-60,20},{-60,
          15.6},{-49.2,15.6}}, color={0,0,127}));
  connect(in_setPoint2, add.u2) annotation (Line(points={{20,112},{20,98},{-72,98},
          {-72,8},{-60,8},{-60,8.4},{-49.2,8.4}},
                                        color={0,0,127}));
  connect(add.y, rampMultiplier.u2) annotation (Line(points={{-35.4,12},{-32,12},
          {-32,17},{-26,17}}, color={0,0,127}));
  connect(in_setPoint2, add1.u2) annotation (Line(points={{20,112},{20,98},{-72,
          98},{-72,-3.6},{6.8,-3.6}},                         color={0,0,127}));
  connect(triggeredTrapezoid.y, rampMultiplier.u1) annotation (Line(points={{-43,
          42},{-34,42},{-34,23},{-26,23}}, color={0,0,127}));
  connect(not1.y, triggeredTrapezoid.u) annotation (Line(points={{-75.2,42},{-66,42}}, color={255,0,255}));
  connect(not1.u, in_setPoint2Active)
    annotation (Line(points={{-85.6,42},{-100,42},{-100,100},{-20,100},{-20,110}},
                                                                         color={255,0,255}));
  connect(add1.y, out_setPointVariable)
    annotation (Line(points={{20.6,0},{110,0}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                                        Text(
        extent={{-150,-100},{150,-140}},
        textString="%name",
        lineColor={0,0,255}), Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={28,108,200},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid)}),                      Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
The model is created to achieve a smooth transition between two signals: it is originally tought for manual/automatic tracking for a controller, but then generalized to take any kind of signals.
It requires two real inputs (setPoint1 and setPoint2), and a boolean input (in_setPoint2Active) that when TRUE will use setPoint2, and when false setPoint1. 
The transition is regulated by a ramp that smoothen the differences between the signal, and the tuning parameter is the duration of said ramp (rampCompensatorDuration)
</html>", revisions="<html>
<ul>
<li><i>19-11-2021</i> by <a href=\"mailto:daniele.fortunati@dlr.de\">Daniele Fortunati</a>: <br> First release</li>
<li><i>19-01-2022</i> by <a href=\"mailto:marius.tomberg@dlr.de\">Marius Tomberg</a>: <br>Icon updated</li>
<li><i>19-11-2024</i> by <a href=\"mailto:rene.lorenz@dlr.de\">René Lorenz</a>: <br>Bugfix: transition was immediate not smooth (lin. ramp). And in_setPoint2Active is always required now</li>
</ul>
</html>"));
end BumplessTransferController;

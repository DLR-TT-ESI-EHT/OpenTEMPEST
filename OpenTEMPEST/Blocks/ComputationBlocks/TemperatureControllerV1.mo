within OpenTEMPEST.Blocks.ComputationBlocks;
model TemperatureControllerV1

  parameter Modelica.SIunits.Temperature TFuelMax = 900;
  parameter Modelica.SIunits.Temperature TFuelMin = 650;
  parameter Modelica.SIunits.Temperature TAirMax = 900;
  parameter Modelica.SIunits.Temperature TAirMin = 600;
  parameter Modelica.SIunits.MassFlowRate mfAirMax = 25/1000;
  parameter Modelica.SIunits.MassFlowRate mfAirMin = 4/1000;

  parameter Modelica.SIunits.Temperature TITIn = 825;
  parameter Modelica.SIunits.MassFlowRate mfITAir = 4/1000;

  Modelica.SIunits.Temperature TFuel_internal;
  Modelica.SIunits.Temperature TAir_internal;
  Modelica.SIunits.MassFlowRate mfAir_internal;
  Modelica.Blocks.Interfaces.RealInput in_TITIn
    "1 for maximum heating and 0 for maximum cooling" annotation (Placement(
        transformation(extent={{-118,24},{-78,64}}), iconTransformation(extent=
            {{-140,34},{-100,74}})));

  Modelica.Blocks.Interfaces.RealInput u
    "1 for maximum heating and 0 for maximum cooling"
                                         annotation (Placement(transformation(
          extent={{-140,-20},{-100,20}}), iconTransformation(extent={{-140,-20},
            {-100,20}})));
  Modelica.Blocks.Interfaces.RealOutput TFuel annotation (Placement(
        transformation(extent={{100,60},{120,80}}), iconTransformation(extent={{
            100,60},{120,80}})));
  Modelica.Blocks.Interfaces.RealOutput TAir annotation (Placement(
        transformation(extent={{100,-10},{120,10}}), iconTransformation(extent={
            {100,-10},{120,10}})));
  Modelica.Blocks.Interfaces.RealOutput mfAir annotation (Placement(
        transformation(extent={{100,-80},{120,-60}}), iconTransformation(extent=
           {{100,-80},{120,-60}})));

  Modelica.Blocks.Interfaces.BooleanInput in_IT
    "True isothermal operation if desired" annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={0,-98}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={0,-120})));
  BumplessTransferController bumplessTFuel(
      rampCompensatorDuration=600)
    annotation (Placement(transformation(extent={{48,60},{68,80}})));
  BumplessTransferController bumplessTAir(
      rampCompensatorDuration=600)
    annotation (Placement(transformation(extent={{48,-10},{68,10}})));
  BumplessTransferController bumplessMfAir(
      rampCompensatorDuration=600)
    annotation (Placement(transformation(extent={{48,-80},{68,-60}})));
  Modelica.Blocks.Sources.RealExpression realExpression(y=TFuel_internal)
    annotation (Placement(transformation(extent={{0,62},{20,82}})));
  Modelica.Blocks.Sources.RealExpression realExpression1(y=TAir_internal)
    annotation (Placement(transformation(extent={{0,-8},{20,12}})));
  Modelica.Blocks.Sources.RealExpression realExpression2(y=mfAir_internal)
    annotation (Placement(transformation(extent={{0,-78},{20,-58}})));
  Modelica.Blocks.Sources.Constant const(k=TITIn)
    annotation (Placement(transformation(extent={{18,22},{38,42}})));
  Modelica.Blocks.Sources.Constant const1(k=mfITAir)
    annotation (Placement(transformation(extent={{20,-50},{40,-30}})));
equation
  // input value u is 1 for maximum heating and 0 for maximum cooling
  TFuel_internal = (TFuelMax-TFuelMin)*u+TFuelMin;
  TAir_internal  = (TAirMax -TAirMin) *u+TAirMin;
  mfAir_internal = (mfAirMax-mfAirMin)*(4*u^2-4*u+1)+mfAirMin;

  connect(bumplessTFuel.out_setPointVariable, TFuel)
    annotation (Line(points={{69,70},{110,70}}, color={0,0,127}));
  connect(TAir, bumplessTAir.out_setPointVariable)
    annotation (Line(points={{110,0},{69,0}}, color={0,0,127}));
  connect(mfAir, bumplessMfAir.out_setPointVariable)
    annotation (Line(points={{110,-70},{69,-70}}, color={0,0,127}));
  connect(bumplessTFuel.in_setPoint1, realExpression.y)
    annotation (Line(points={{46,70},{34,70},{34,72},{21,72}},
                                               color={0,0,127}));
  connect(bumplessMfAir.in_setPoint1, realExpression2.y)
    annotation (Line(points={{46,-70},{34,-70},{34,-68},{21,-68}},
                                                 color={0,0,127}));
  connect(bumplessTAir.in_setPoint1, realExpression1.y)
    annotation (Line(points={{46,0},{34,0},{34,2},{21,2}},
                                             color={0,0,127}));
  connect(in_IT, bumplessMfAir.in_setPoint2Active) annotation (Line(points={{0,-98},
          {30,-98},{30,-92},{54,-92},{54,-58}},      color={255,0,255}));
  connect(bumplessMfAir.in_setPoint2Active, bumplessTAir.in_setPoint2Active)
    annotation (Line(points={{54,-58},{54,-40},{52,-40},{52,-23},{54,-23},{54,
          12}},                                 color={255,0,255}));
  connect(bumplessTAir.in_setPoint2Active, bumplessTFuel.in_setPoint2Active)
    annotation (Line(points={{54,12},{54,30},{52,30},{52,47},{54,47},{54,82}},
                                               color={255,0,255}));
  connect(in_TITIn, bumplessTFuel.in_setPoint2) annotation (Line(points={{-98,44},
          {-16,44},{-16,50},{62,50},{62,82}},       color={0,0,127}));
  connect(bumplessTAir.in_setPoint2, bumplessTFuel.in_setPoint2)
    annotation (Line(points={{62,12},{62,30},{64,30},{64,46.8},{62,46.8},{62,82}},
                                                   color={0,0,127}));
  connect(const1.y, bumplessMfAir.in_setPoint2) annotation (Line(points={{41,-40},
          {54,-40},{54,-38},{62,-38},{62,-58}},        color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                                Rectangle(
        extent={{-100,-100},{100,100}},
        lineColor={0,0,127},
        fillColor={255,255,255},
        fillPattern=FillPattern.Solid)}), Diagram(coordinateSystem(
          preserveAspectRatio=false)),
    Documentation(revisions="<html>
<ul>
<li>10.08.2022 by Marius Tomberg<br>Isothermal option added</li>
</ul>
</html>"));
end TemperatureControllerV1;

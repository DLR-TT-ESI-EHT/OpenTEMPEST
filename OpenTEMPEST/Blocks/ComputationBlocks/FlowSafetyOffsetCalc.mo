within OpenTEMPEST.Blocks.ComputationBlocks;
model FlowSafetyOffsetCalc
  "Increases/descreases mass flow depending of the trend"

  extends Modelica.Blocks.Interfaces.SISO;

  parameter Modelica.SIunits.Time tPrediction = 30 "Time the offset is calcualted for";
  parameter Modelica.SIunits.Time tPt1 = 5 "Time constant of the PT1 block for smoothing";

  parameter Real yStart=0 "Start value of the output";

  parameter Real minVal "Min value of the output";
  parameter Real maxVal "Max value of the output";

  Modelica.Blocks.Continuous.Derivative derivative
    annotation (Placement(transformation(extent={{-80,-40},{-60,-20}})));
  Modelica.Blocks.Sources.Constant const(k=tPrediction)
    annotation (Placement(transformation(extent={{-80,-70},{-60,-50}})));
  Modelica.Blocks.Math.MultiProduct
                               multiProduct(nu=4)
    annotation (Placement(transformation(extent={{-18,-46},{2,-26}})));
  Modelica.Blocks.Math.Add add
    annotation (Placement(transformation(extent={{12,-10},{32,10}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder(T=tPt1, y_start=yStart)
    annotation (Placement(transformation(extent={{76,-10},{96,10}})));
  Modelica.Blocks.Math.Sign signU
    annotation (Placement(transformation(extent={{-44,60},{-24,80}})));
  Modelica.Blocks.Math.Sign signDerU
    annotation (Placement(transformation(extent={{-44,30},{-24,50}})));
  LimitBlock limitBlock(minVal=minVal, maxVal=maxVal) annotation (Placement(
        transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={54,0})));
equation
  connect(add.u2, multiProduct.y)
    annotation (Line(points={{10,-6},{10,-36},{3.7,-36}}, color={0,0,127}));
  connect(u, derivative.u) annotation (Line(points={{-120,0},{-92,0},{-92,-30},{
          -82,-30}}, color={0,0,127}));
  connect(u, add.u1) annotation (Line(points={{-120,0},{-92,0},{-92,6},{10,6}},
        color={0,0,127}));
  connect(firstOrder.y, y)
    annotation (Line(points={{97,0},{110,0}}, color={0,0,127}));
  connect(signDerU.u, derivative.y) annotation (Line(points={{-46,40},{-52,40},{
          -52,-30},{-59,-30}}, color={0,0,127}));
  connect(signU.u, u) annotation (Line(points={{-46,70},{-92,70},{-92,0},{-120,0}},
        color={0,0,127}));
  connect(const.y, multiProduct.u[1]) annotation (Line(points={{-59,-60},{-40,-60},
          {-40,-30.75},{-18,-30.75}}, color={0,0,127}));
  connect(derivative.y, multiProduct.u[2]) annotation (Line(points={{-59,-30},{-40,
          -30},{-40,-34.25},{-18,-34.25}}, color={0,0,127}));
  connect(signDerU.y, multiProduct.u[3]) annotation (Line(points={{-23,40},{-18,
          40},{-18,-37.75}}, color={0,0,127}));
  connect(signU.y, multiProduct.u[4]) annotation (Line(points={{-23,70},{-18,70},
          {-18,-41.25}}, color={0,0,127}));
  connect(add.y, limitBlock.u) annotation (Line(points={{33,0},{37.5,0},{37.5,1.9984e-15},
          {42,1.9984e-15}}, color={0,0,127}));
  connect(limitBlock.y, firstOrder.u) annotation (Line(points={{65,-8.88178e-16},
          {69.5,-8.88178e-16},{69.5,0},{74,0}}, color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=2400,
      Tolerance=1e-05,
      __Dymola_Algorithm="Cvode"),
    __Dymola_experimentSetupOutput,
    __Dymola_experimentFlags(
      Advanced(GenerateVariableDependencies=false, OutputModelicaCode=true),
      Evaluate=false,
      OutputCPUtime=true,
      OutputFlatModelica=false));
end FlowSafetyOffsetCalc;

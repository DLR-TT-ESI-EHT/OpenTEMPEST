within OpenTEMPEST.Flow;
model ForcedFlowSplitter2
  parameter Real splitRatio = 0.5 "flow ratio going out of flange  B1";
  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
    annotation(choicesAllMatching = true);

  parameter Modelica.SIunits.Pressure pStart=101325;
  parameter Modelica.SIunits.Temperature Tstart=300;
  parameter Modelica.SIunits.MassFraction Xstart[Medium.nX]=Medium.reference_X;

  ThermoPower.Gas.FlangeA flangeA(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
  ThermoPower.Gas.FlangeB flangeB1(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{40,10},{60,30}})));
  ThermoPower.Gas.FlangeB flangeB2(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{40,-30},{60,-10}})));
  ThermoPower.Gas.SensW sensW(redeclare package Medium = Medium,
      allowFlowReversal=false)
    annotation (Placement(transformation(extent={{-40,-6},{-20,14}})));
  Modelica.Blocks.Sources.Constant const(k=splitRatio)
    annotation (Placement(transformation(extent={{-78,32},{-58,52}})));
  Modelica.Blocks.Math.Product product
    annotation (Placement(transformation(extent={{-32,26},{-12,46}})));
  Modelica.Blocks.Sources.Constant const1(k=1 - splitRatio)
    annotation (Placement(transformation(extent={{-44,-24},{-36,-16}})));
  Modelica.Blocks.Math.Product product1
    annotation (Placement(transformation(extent={{-18,-26},{-10,-18}})));
  ThermoPower.Gas.SourceMassFlow sourceMassFlow1(
    redeclare package Medium = Medium,
    p0=pStart,
    T=Tstart,
    Xnom=Xstart,
    use_in_w0=true,
    use_in_T=true,
    use_in_X=true)
    annotation (Placement(transformation(extent={{12,26},{24,14}})));
  ThermoPower.Gas.SinkPressure sinkPressure(redeclare package Medium =
        Medium,
    p0=pStart,
    T=Tstart,
    Xnom=Xstart,                 use_in_p0=true)
    annotation (Placement(transformation(extent={{-4,6},{8,-6}})));
  ThermoPower.Gas.SourceMassFlow sourceMassFlow2(
    redeclare package Medium = Medium,
    p0=pStart,
    T=Tstart,
    Xnom=Xstart,
    use_in_w0=true,
    use_in_T=true,
    use_in_X=true)
    annotation (Placement(transformation(extent={{12,-26},{24,-14}})));
  ThermoPower.Gas.SensP sensP(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{42,24},{28,10}})));
  ThermoPower.Gas.SensT sensT(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-20,4},{-6,-10}})));
equation

  sourceMassFlow1.in_X = sensW.outlet.Xi_outflow;
  sourceMassFlow2.in_X = sensW.outlet.Xi_outflow;

  connect(sensW.inlet, flangeA) annotation (Line(points={{-36,0},{-50,0}},
                    color={159,159,223}));
  connect(const.y, product.u1)
    annotation (Line(points={{-57,42},{-34,42}}, color={0,0,127}));
  connect(product.u2, sensW.w) annotation (Line(points={{-34,30},{-48,30},{-48,20},
          {-16,20},{-16,10},{-23,10}},
                     color={0,0,127}));
  connect(flangeA, flangeA) annotation (Line(points={{-50,0},{-50,2},{-50,2},
          {-50,0}}, color={159,159,223}));
  connect(sinkPressure.flange, sensT.outlet) annotation (Line(points={{-4,0},{-6.4,
          0},{-6.4,-0.2},{-8.8,-0.2}}, color={159,159,223}));
  connect(sensT.inlet, sensW.outlet) annotation (Line(points={{-17.2,-0.2},{-20.6,
          -0.2},{-20.6,0},{-24,0}}, color={159,159,223}));
  connect(sourceMassFlow1.flange, sensP.flange) annotation (Line(points={{
          24,20},{30,20},{30,19.8},{35,19.8}}, color={159,159,223}));
  connect(flangeB1, sourceMassFlow1.flange)
    annotation (Line(points={{50,20},{24,20}}, color={159,159,223}));
  connect(sourceMassFlow2.flange, flangeB2)
    annotation (Line(points={{24,-20},{50,-20}}, color={159,159,223}));
  connect(sourceMassFlow1.in_w0, product.y) annotation (Line(points={{14.4,
          17},{2.2,17},{2.2,36},{-11,36}}, color={0,0,127}));
  connect(product1.u2, sensW.w) annotation (Line(points={{-18.8,-24.4},{-18.8,-26.2},
          {-23,-26.2},{-23,10}}, color={0,0,127}));
  connect(product1.u1, const1.y) annotation (Line(points={{-18.8,-19.6},{-26.4,-19.6},
          {-26.4,-20},{-35.6,-20}}, color={0,0,127}));
  connect(product1.y, sourceMassFlow2.in_w0) annotation (Line(points={{-9.6,
          -22},{2,-22},{2,-17},{14.4,-17}}, color={0,0,127}));
  connect(sensP.p, sinkPressure.in_p0) annotation (Line(points={{30.1,12.8},{14.05,
          12.8},{14.05,-3.57},{-1.87,-3.57}}, color={0,0,127}));
  connect(sourceMassFlow2.in_T, sensT.T) annotation (Line(points={{18,-17},
          {6,-17},{6,-7.2},{-8.1,-7.2}}, color={0,0,127}));
  connect(sourceMassFlow1.in_T, sensT.T) annotation (Line(points={{18,17},{
          6,17},{6,-7.2},{-8.1,-7.2}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
          extent={{-40,32},{40,-30}},
          fillColor={170,213,255},
          fillPattern=FillPattern.Solid,
          pattern=LinePattern.None)}),                           Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end ForcedFlowSplitter2;

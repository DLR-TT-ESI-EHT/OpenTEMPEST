within OpenTEMPEST.Flow;
model ForcedFlowMultiSplitter

  parameter Integer nOutlets(min = 2)=3;
  //parameter Real flowCoeff[nOutlets-1] = {1,1} "mulitplicator for the flows of the first n-1 outlets";

  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
    annotation(choicesAllMatching = true);

  ThermoPower.Gas.FlangeA inlet(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
  ThermoPower.Gas.FlangeB outlet[nOutlets](redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{40,-10},{60,10}})));
  ThermoPower.Gas.ThroughMassFlow throughMassFlow[nOutlets-1](
    redeclare package Medium = Medium,
    each w0=0.00001,
    each use_in_w0=true)
    annotation (Placement(transformation(extent={{12,-10},{32,10}})));
  ThermoPower.Gas.SensW sensW(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-34,-6},{-14,14}})));
  Modelica.Blocks.Sources.Constant const(k=1/nOutlets)
    annotation (Placement(transformation(extent={{-100,32},{-80,52}})));
  Modelica.Blocks.Math.Product product
    annotation (Placement(transformation(extent={{-42,26},{-22,46}})));
  OpenTEMPEST.Flow.MultiSplitter multiSplitter(nOutlets=nOutlets, redeclare
      package Medium = Medium)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
  Modelica.Blocks.Math.Product productVar[nOutlets - 1]
    annotation (Placement(transformation(extent={{-10,32},{10,52}})));
  Modelica.Blocks.Sources.Constant constVar[nOutlets - 1](k=1)
    annotation (Placement(transformation(extent={{-42,56},{-22,76}})));
equation

  for i in 1:nOutlets-1 loop
    connect(multiSplitter.outlet[i],throughMassFlow[i].inlet);
    connect(throughMassFlow[i].outlet, outlet[i]);
    //connect(throughMassFlow[i].in_w0, product.y);
    connect(product.y, productVar[i].u2)
    annotation (Line(points={{-21,36},{-12,36}}, color={0,0,127}));
  end for;
  connect(multiSplitter.outlet[nOutlets],outlet[nOutlets]);

  connect(sensW.inlet, inlet) annotation (Line(points={{-30,0},{-50,0}}, color={159,159,223}));
  connect(const.y, product.u1)
    annotation (Line(points={{-79,42},{-44,42}}, color={0,0,127}));
  connect(product.u2, sensW.w) annotation (Line(points={{-44,30},{-50,30},{-50,
          20},{-10,20},{-10,10},{-17,10}},
                     color={0,0,127}));
  connect(multiSplitter.inlet, sensW.outlet) annotation (Line(points={{-5,0},{
          -18,0}},                                                                      color={159,159,223}));
  connect(productVar.y, throughMassFlow.in_w0)
    annotation (Line(points={{11,42},{16,42},{16,5}},        color={0,0,127}));

  connect(constVar.y, productVar.u1) annotation (Line(points={{-21,66},{-12,66},
          {-12,48}},          color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
          extent={{-40,32},{40,-30}},
          fillColor={170,213,255},
          fillPattern=FillPattern.Solid,
          pattern=LinePattern.None), Text(
          extent={{-40,32},{40,-30}},
          lineColor={28,108,200},
          textString="FMS")}),                                   Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end ForcedFlowMultiSplitter;

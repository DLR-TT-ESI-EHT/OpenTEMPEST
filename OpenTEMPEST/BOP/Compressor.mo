within OpenTEMPEST.BOP;
model Compressor

  //parameter Real Pi "Pressure ratio";
  parameter Real eta=0.8 "Polytropic efficiency";
  parameter Real etaM = 0.98 "Mechanical efficiency";
  parameter Real etaE = 0.98 "Electrical efficiency";

  replaceable package Medium = OpenTEMPEST.Medium.Air_Medium
                                                         constrainedby
    Modelica.Media.Interfaces.PartialMedium "Medium model"
    annotation(choicesAllMatching = true);

  Medium.BaseProperties gasIn;
  Medium.BaseProperties gasOut;

  Real kappaIn;
  Real kappaOut;

  Real nFrac;

  ThermoPower.Gas.FlangeA inlet(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
  ThermoPower.Gas.FlangeB outlet(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{80,-10},{100,10}})));

  ThermoPower.Electrical.Load load(Pnom=10000, usePowerInput=true)
    annotation (Placement(transformation(extent={{-10,40},{10,60}})));
  ThermoPower.Electrical.PowerConnection powerConnection
    annotation (Placement(transformation(extent={{-40,60},{-20,80}}),
        iconTransformation(extent={{-40,60},{-20,80}})));
  ThermoPower.Thermal.HT heatLosses
    annotation (Placement(transformation(extent={{20,60},{40,80}}),
        iconTransformation(extent={{20,60},{40,80}})));
  Modelica.Blocks.Sources.RealExpression realExpressionPEl(y=inlet.m_flow*(
        gasOut.h - gasIn.h)/etaM/etaE)
    annotation (Placement(transformation(extent={{-86,40},{-24,60}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow(T_ref=
        323.15, alpha=0.1)
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={0,-42})));
  Modelica.Blocks.Sources.RealExpression realExpressionPEl1(y=-(load.P - inlet.m_flow
        *(gasOut.h - gasIn.h)))
    annotation (Placement(transformation(extent={{-82,-42},{-20,-22}})));
  Modelica.Blocks.Interfaces.RealInput Pi "Pressure ratio" annotation (
      Placement(transformation(extent={{-10,-100},{30,-60}}),
        iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={0,-80})));
equation

  // Independent composition mass balances
  inlet.Xi_outflow = inStream(outlet.Xi_outflow);
  inStream(inlet.Xi_outflow) = outlet.Xi_outflow;

  inlet.m_flow + outlet.m_flow = 0 "Mass balance";

  // Pressure increase
  outlet.p/inlet.p = Pi "Pressure ratio";

  // Flow work
  //  y = eta * (gasOut.h-gasIn.h)

  //load.P = inlet.m_flow*(gasOut.h-gasIn.h)/etaM/etaE;
  //heatLosses.Q_flow = load.P - inlet.m_flow*(gasOut.h-gasIn.h);

  kappaIn = Medium.isentropicExponent(gasIn.state);
  kappaOut = Medium.isentropicExponent(gasOut.state);

  nFrac = eta * 0.5*(kappaIn/(kappaIn-1)+kappaOut/(kappaOut-1));

  gasOut.T = Pi^(1/nFrac) * gasIn.T;

  // Connectors

  gasIn.h = inStream(inlet.h_outflow);
  gasIn.p = inlet.p;
  gasIn.Xi = inStream(inlet.Xi_outflow);

  gasOut.h = outlet.h_outflow;
  gasOut.p = outlet.p;
  gasOut.Xi = outlet.Xi_outflow;

  inlet.h_outflow = inStream(outlet.h_outflow);

  connect(load.port, powerConnection) annotation (Line(
      points={{0,58.6},{0,70},{-30,70}},
      color={0,0,255},
      thickness=0.5));
  connect(realExpressionPEl.y, load.referencePower)
    annotation (Line(points={{-20.9,50},{-3.3,50}},
                                                  color={0,0,127}));
  connect(heatLosses, prescribedHeatFlow.port)
    annotation (Line(points={{30,70},{30,10},{0,10},{0,-52}},
                                               color={191,0,0}));
  connect(realExpressionPEl1.y, prescribedHeatFlow.Q_flow)
    annotation (Line(points={{-16.9,-32},{0,-32}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Ellipse(
          extent={{-60,60},{60,-60}},
          lineColor={0,0,0},
          lineThickness=1,
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-44,40},{56,20}},
          color={0,0,0},
          thickness=1),
        Line(
          points={{-44,-40},{56,-20}},
          color={0,0,0},
          thickness=1),
        Line(
          points={{-48,22}},
          color={0,0,0},
          thickness=1),
        Line(
          points={{60,0},{94,0}},
          color={0,0,0},
          thickness=1),
        Line(
          points={{-60,0},{-92,0}},
          color={0,0,0},
          thickness=1)}),                                        Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end Compressor;

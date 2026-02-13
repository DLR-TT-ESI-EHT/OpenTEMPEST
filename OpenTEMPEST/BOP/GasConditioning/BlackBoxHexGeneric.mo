within OpenTEMPEST.BOP.GasConditioning;
model BlackBoxHexGeneric

   import SI = Modelica.SIunits;

  replaceable package MediumA = Medium.Air_Medium         constrainedby
    Modelica.Media.Interfaces.PartialMedium
    annotation(choicesAllMatching = true);
  replaceable package MediumB = Medium.Air_Medium         constrainedby
    Modelica.Media.Interfaces.PartialMedium
    annotation(choicesAllMatching = true);

  parameter SI.HeatCapacity C = 1;
  replaceable model etaType = HexEtaTypes.ConstantEta                                        constrainedby
    HexEtaTypes.EtaBase
      annotation (choicesAllMatching=true);

  etaType eta(C_min=Cmin, C_max=Cmax, mf_A=sensAIn.mf, mf_B=sensBIn.mf, lambda_A=sensAIn.lambda, lambda_B=sensBIn.lambda, mu_A=sensAIn.eta, mu_B=sensBIn.eta,
              cp_A=sensAIn.cp, cp_B=sensBIn.cp, redeclare package MediumA = MediumA, redeclare
      package MediumB =                                                                                          MediumB);
  Real Cmin(start=1);
  Real Cmax(start=1);

  parameter Real Tstart=700;
  parameter SI.MassFraction xStartA[MediumA.nX]=MediumA.reference_X;
  parameter SI.MassFraction xStartB[MediumB.nX]=MediumB.reference_X;

  parameter Real y_start_w0 "Initial value for mass flow medium A";
  parameter Real y_start_p0 "Initial value for pressure medium B";

  SI.Temperature TAIn=sensAIn.T;
  SI.Temperature TBIn=sensBIn.T;
  SI.Temperature TAOut(start=Tstart);
  SI.Temperature TBOut(start=Tstart);

  SI.HeatFlowRate Qf;
  SI.HeatFlowRate QmaxA, QmaxB;

  ThermoPower.Gas.FlangeA inletB(redeclare package Medium = MediumB) annotation (
      Placement(transformation(rotation=0, extent={{-80,-40},{-60,-20}}),
        iconTransformation(extent={{-100,-60},{-60,-20}})));
  ThermoPower.Gas.FlangeA inletA(redeclare package Medium = MediumA)
                                                                  annotation (
      Placement(transformation(rotation=0, extent={{-80,40},{-60,60}}),
        iconTransformation(extent={{60,20},{100,60}})));
  ThermoPower.Gas.FlangeB outletA(redeclare package Medium = MediumA) annotation (
      Placement(transformation(rotation=0, extent={{80,40},{100,60}}),
        iconTransformation(extent={{-100,20},{-60,60}})));
  ThermoPower.Gas.FlangeB outletB(redeclare package Medium = MediumB) annotation (
      Placement(transformation(rotation=0, extent={{80,-40},{100,-20}}),
        iconTransformation(extent={{60,-60},{100,-20}})));
  OpenTEMPEST.Flow.SensGasProperty sensBOut(
    mfOutput=false,
    pOutput=true,
    hOutput=false,
    XOutput=false,
    YOutput=true,
    HfOutput=true,
    redeclare package Medium = MediumB,
    Tstart=Tstart,
    Xstart=xStartB,
    cpOutput=true,
    lambdaOutput=false,
    etaOutput=false,
    TOutput=true,
    rhoOutput=false)
    annotation (Placement(transformation(extent={{52,-36},{72,-16}})));
  OpenTEMPEST.Flow.SensGasProperty sensBIn(
    mfOutput=true,
    pOutput=false,
    hOutput=false,
    XOutput=true,
    YOutput=false,
    HfOutput=true,
    redeclare package Medium = MediumB,
    Tstart=Tstart,
    Xstart=xStartB,
    cpOutput=true,
    lambdaOutput=true,
    etaOutput=true,
    TOutput=true,
    rhoOutput=false)
    annotation (Placement(transformation(extent={{-52,-36},{-32,-16}})));
  OpenTEMPEST.Flow.SensGasProperty sensAIn(
    mfOutput=true,
    pOutput=false,
    hOutput=false,
    XOutput=true,
    YOutput=true,
    HfOutput=true,
    redeclare package Medium = MediumA,
    Tstart=Tstart,
    Xstart=xStartA,
    cpOutput=true,
    lambdaOutput=true,
    etaOutput=true,
    TOutput=true,
    rhoOutput=false)
    annotation (Placement(transformation(extent={{-56,44},{-36,64}})));
  OpenTEMPEST.Flow.SensGasProperty sensAOut(
    mfOutput=false,
    pOutput=true,
    hOutput=false,
    XOutput=false,
    YOutput=true,
    HfOutput=true,
    redeclare package Medium = MediumA,
    Tstart=Tstart,
    Xstart=xStartA,
    cpOutput=true,
    lambdaOutput=false,
    etaOutput=false,
    TOutput=true,
    rhoOutput=false)
    annotation (Placement(transformation(extent={{46,44},{66,64}})));
  ThermoPower.Gas.SinkPressure sinkPressureA(
    redeclare package Medium = MediumA,
    p0=100000,
    T=Tstart,
    Xnom=xStartA,
    use_in_p0=true,
    use_in_T=false)
    annotation (Placement(transformation(extent={{-26,40},{-6,60}})));
  ThermoPower.Gas.SourceMassFlow sourceMassFlowA(
    redeclare package Medium = MediumA,
    p0=100000,
    T=Tstart,
    Xnom=xStartA,
    use_in_w0=true,
    use_in_T=true,
    use_in_X=true)
    annotation (Placement(transformation(extent={{18,40},{38,60}})));
  ThermoPower.Gas.SinkPressure sinkPressureB(
    redeclare package Medium = MediumB,
    p0=100000,
    T=Tstart,
    Xnom=xStartB,
    use_in_p0=true,
    use_in_T=false)
    annotation (Placement(transformation(extent={{-16,-40},{4,-20}})));
  ThermoPower.Gas.SourceMassFlow sourceMassFlowB(
    redeclare package Medium = MediumB,
    p0=100000,
    T=Tstart,
    Xnom=xStartB,
    use_in_w0=true,
    use_in_T=true,
    use_in_X=true)
    annotation (Placement(transformation(extent={{20,-40},{40,-20}})));
  ThermoPower.Thermal.HT htAmbA annotation (Placement(transformation(extent={{-10,84},{10,104}}), iconTransformation(extent={{-6,60},{6,72}})));
  Modelica.Blocks.Sources.RealExpression rexTAOut(y=TAOut)
    annotation (Placement(transformation(extent={{4,18},{24,38}})));
  Modelica.Blocks.Sources.RealExpression rexTBOut(y=TBOut)
    annotation (Placement(transformation(extent={{4,4},{24,24}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder(T=1e-5, y_start=y_start_w0)
                                                                         annotation (Placement(transformation(extent={{-16,70},{-4,82}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder1(T=1e-5, y_start=y_start_p0)
                                                                          annotation (Placement(transformation(extent={{50,-8},{38,4}})));
  ThermoPower.Thermal.HT htAmbB annotation (Placement(transformation(extent={{-10,-104},{10,-84}}), iconTransformation(extent={{-6,-72},{6,-60}})));
equation

  htAmbA.Q_flow = 0.15*(htAmbA.T - TAOut);
  htAmbB.Q_flow = 0.15*(htAmbB.T - TBOut);

  QmaxA = sensAIn.Hf - sensAIn.mf*MediumA.h_TX(TBIn, sensAIn.x);
  QmaxB = sensBIn.Hf - sensBIn.mf*MediumB.h_TX(TAIn, sensBIn.x);

  Cmin*abs(TAIn-TBIn) = min(abs({QmaxA,QmaxB}));
  Cmax*abs(TAIn-TBIn) = max(abs({QmaxA,QmaxB}));

  Qf=eta.eta*Cmin*(TAIn-TBIn);

  // 2 thermal masses neglecting gas capacity in (Ch+1/2*Cs)dTh/dt=H_h_in-H_h_out-Q+Qamb
  0.5*C*der(TAOut) = sensAIn.Hf - sensAIn.mf*MediumA.h_TX(TAOut, sensAIn.x) - Qf + htAmbA.Q_flow;
  0.5*C*der(TBOut) = sensBIn.Hf - sensBIn.mf*MediumB.h_TX(TBOut, sensBIn.x) + Qf + htAmbB.Q_flow;

  connect(sensBOut.outlet, outletB)
    annotation (Line(points={{68,-30},{90,-30}}, color={159,159,223}));
  connect(inletB, sensBIn.inlet)
    annotation (Line(points={{-70,-30},{-48,-30}}, color={159,159,223}));
  connect(inletA, sensAIn.inlet)
    annotation (Line(points={{-70,50},{-52,50}}, color={159,159,223}));
  connect(sensAOut.outlet, outletA)
    annotation (Line(points={{62,50},{90,50}}, color={159,159,223}));
  connect(sensAIn.outlet, sinkPressureA.flange)
    annotation (Line(points={{-40,50},{-26,50}}, color={159,159,223}));
  connect(sensAOut.inlet, sourceMassFlowA.flange)
    annotation (Line(points={{50,50},{38,50}}, color={159,159,223}));
  connect(sensBIn.outlet, sinkPressureB.flange)
    annotation (Line(points={{-36,-30},{-16,-30}}, color={159,159,223}));
  connect(sensBOut.inlet, sourceMassFlowB.flange)
    annotation (Line(points={{56,-30},{40,-30}}, color={159,159,223}));

  connect(sensBIn.x, sourceMassFlowB.in_X) annotation (Line(points={{-35,-23},{-25.5,
          -23},{-25.5,-25},{36,-25}}, color={0,0,127}));
  connect(sensAIn.x, sourceMassFlowA.in_X) annotation (Line(points={{-39,57},{36,
          57},{36,55},{34,55}}, color={0,0,127}));
  connect(sensBIn.mf, sourceMassFlowB.in_w0)
    annotation (Line(points={{-35,-15},{24,-15},{24,-25}}, color={0,0,127}));
  connect(rexTAOut.y, sourceMassFlowA.in_T) annotation (Line(points={{25,28},{30,
          28},{30,55},{28,55}}, color={0,0,127}));
  connect(rexTBOut.y, sourceMassFlowB.in_T) annotation (Line(points={{25,14},{28,
          14},{28,-25},{30,-25}}, color={0,0,127}));
  connect(sensAOut.p, sinkPressureA.in_p0) annotation (Line(points={{63,59},{63,64},{-22.45,64},{-22.45,55.95}}, color={0,0,127}));
  connect(sensAIn.mf, firstOrder.u) annotation (Line(points={{-39,65},{-26,65},{-26,76},{-17.2,76}}, color={0,0,127}));
  connect(firstOrder.y, sourceMassFlowA.in_w0) annotation (Line(points={{-3.4,76},{22,76},{22,55}}, color={0,0,127}));
  connect(sensBOut.p, firstOrder1.u) annotation (Line(points={{69,-21},{69,-2},{51.2,-2}}, color={0,0,127}));
  connect(firstOrder1.y, sinkPressureB.in_p0) annotation (Line(points={{37.4,-2},{-12.45,-2},{-12.45,-24.05}}, color={0,0,127}));

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Rectangle(
          extent={{-80,60},{80,-60}},
          lineColor={28,108,200},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-80,12},{80,-6}},
          lineColor={28,108,200},
          fillColor={175,175,175},
          fillPattern=FillPattern.CrossDiag),
        Text(
          extent={{70,-10},{-72,16}},
          lineColor={28,108,200},
          textString="HEX",
          textStyle={TextStyle.Bold})}),                         Diagram(coordinateSystem(preserveAspectRatio=false)));
end BlackBoxHexGeneric;

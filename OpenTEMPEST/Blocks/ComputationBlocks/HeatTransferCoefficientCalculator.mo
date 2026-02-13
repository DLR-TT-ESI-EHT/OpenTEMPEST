within OpenTEMPEST.Blocks.ComputationBlocks;
model HeatTransferCoefficientCalculator

  import SI = Modelica.SIunits;

  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
    annotation(choicesAllMatching = true);
  parameter Medium.AbsolutePressure p=101325 "Sink pressure";
  parameter Medium.Temperature T=300 "Source temperature";
  parameter Medium.MassFraction X[Medium.nX]=Medium.reference_X
    "Gas composition";
  parameter Medium.MassFlowRate mf=0 "Mass flowrate";

  parameter SI.Distance L "Flow channel length";
  parameter SI.Length H "height of flow channel between two parallel plates";
  parameter SI.Length W "width of flow channel between two parallel plates";
  parameter SI.Length D "diameter of flow channel in a pipe";

  parameter Boolean pipe=false
    "= true if flow through a pipe, false if between two parallel plates";

  Real Nu;
  Real alpha;

  ThermoPower.Gas.Flow1DFV flow1DFV(redeclare package Medium = Medium, N=10,
    L=L,
    A=if pipe then Modelica.Constants.pi*D^2/4 else W*H,
    omega=if pipe then L*Modelica.Constants.pi*D else 2*L*(W + H),
    Dhyd=if pipe then D else H,
    wnom=mf,
    FFtype=ThermoPower.Choices.Flow1D.FFtypes.NoFriction,
    pstart=p,
    Tstartbar=T,
    Xstart=X,
    noInitialPressure=true,
    redeclare model HeatTransfer =
        ThermoPower.Thermal.HeatTransferFV.DittusBoelter)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
  ThermoPower.Gas.SourceMassFlow sourceMassFlow(redeclare package Medium =
        Medium,
    p0=p,
    T=T,
    Xnom=X,
    w0=mf)
    annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
  ThermoPower.Gas.SinkPressure sinkPressure(
    redeclare package Medium = Medium,
    p0=p,
    T=T,
    Xnom=X) annotation (Placement(transformation(extent={{80,-10},{100,10}})));
equation

  // Ebene Platte – laminare Grenzschicht, isotherme Oberfläche
  Nu = 0.664 * flow1DFV.heatTransfer.Re[4]^0.5 * (flow1DFV.heatTransfer.Pr[5])^(1/3);
  alpha = Nu*flow1DFV.heatTransfer.mu[4]/flow1DFV.Dhyd;

  connect(sourceMassFlow.flange, flow1DFV.infl)
    annotation (Line(points={{-80,0},{-10,0}}, color={159,159,223}));
  connect(sinkPressure.flange, flow1DFV.outfl)
    annotation (Line(points={{80,0},{10,0}}, color={159,159,223}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end HeatTransferCoefficientCalculator;

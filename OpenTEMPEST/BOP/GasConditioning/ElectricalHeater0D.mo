within OpenTEMPEST.BOP.GasConditioning;
model ElectricalHeater0D

  import SI = Modelica.SIunits;

  replaceable package Medium = OpenTEMPEST.Medium.Fuel_CH4
                                                        constrainedby
    Modelica.Media.Interfaces.PartialMedium annotation (choicesAllMatching=true);

  parameter SI.Length l=1 "length of flow channel";
  parameter SI.Length d=0.02 "diameter of flow channel";
  parameter Integer nHeatingElements;
  parameter SI.Power Pnom "Nominal active power consumption";

  parameter SI.Temperature TStart = 25+273.15;
  parameter SI.Pressure pStart = 101325;
  parameter SI.MassFraction xStart[Medium.nX] = Medium.reference_X;

//  parameter

  ThermoPower.Gas.FlangeA inlet(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-100,-10},{-80,10}}),
        iconTransformation(extent={{-100,-20},{-60,20}})));
  ThermoPower.Gas.FlangeB outlet(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{80,0},{100,20}}),
        iconTransformation(extent={{60,-20},{100,20}})));
  SimpleGasHeater simpleGasHeater(
    redeclare package medium = Medium,
    l=l,
    d=d/nHeatingElements,
    TStart=TStart,
    pStart=pStart,
    xStart=xStart,
    pipe=true,
    nParallel=nHeatingElements)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
  ThermoPower.Electrical.PowerConnection powerConnection
    annotation (Placement(transformation(extent={{40,60},{60,80}}),
        iconTransformation(extent={{20,40},{60,80}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={0,30})));
  ThermoPower.Electrical.Load load(Pnom=Pnom,
                                   usePowerInput=true)
    annotation (Placement(transformation(extent={{16,40},{36,60}})));
  Modelica.Blocks.Interfaces.RealInput power
    annotation (Placement(transformation(extent={{-70,36},{-30,76}}),
        iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={-40,60})));
equation
  connect(simpleGasHeater.inlet, inlet)
    annotation (Line(points={{-8,0},{-90,0}}, color={159,159,223}));
  connect(simpleGasHeater.outlet, outlet)
    annotation (Line(points={{8,0},{50,0},{50,10},{90,10}},
                                            color={159,159,223}));
  connect(prescribedHeatFlow.port, simpleGasHeater.ht[1]) annotation (Line(
        points={{-1.77636e-15,20},{0,20},{0,4}}, color={191,0,0}));
  connect(powerConnection, load.port) annotation (Line(
      points={{50,70},{26,70},{26,58.6}},
      color={0,0,255},
      thickness=0.5));
  connect(power, load.referencePower) annotation (Line(points={{-50,56},{-13,56},
          {-13,50},{22.7,50}}, color={0,0,127}));
  connect(power, prescribedHeatFlow.Q_flow) annotation (Line(points={{-50,56},{-24,
          56},{-24,40},{1.9984e-15,40}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
          extent={{-60,40},{60,-40}},
          lineColor={28,108,200},
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid)}),                      Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end ElectricalHeater0D;

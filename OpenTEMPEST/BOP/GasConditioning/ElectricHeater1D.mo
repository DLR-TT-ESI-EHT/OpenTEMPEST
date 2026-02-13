within OpenTEMPEST.BOP.GasConditioning;
model ElectricHeater1D

  import SI = Modelica.SIunits;
  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
    annotation(choicesAllMatching = true);
  parameter Integer N(min=3) = 5;
  parameter Integer Ntubes = 1 "number of individual tubes";
  parameter Boolean LUDS = true "Set LU discretisation";

  parameter Modelica.SIunits.Length L "length of the heater tube";
  parameter SI.Length Dhyd "hydraulic diameter of the heater tube";
  parameter SI.Area A "cross sectional area of the heater ";

  //Transport related properties//
  parameter Real  dpdx "nominal pressure drop";

  // intialisation options
  parameter Modelica.SIunits.Temperature TStartIn
    "Inlet temperature to the heater starting value";
  parameter Modelica.SIunits.Temperature TStartOut
    "starting outlet temperature";
  parameter Modelica.SIunits.Pressure pStart
    "start pressure for initialization";
  parameter Modelica.SIunits.MassFraction xStart[Medium.nX]=Medium.reference_X
    "starting composition for initilization";

  // performance parameters
  parameter Real Efficiency(max=1, min=0)
    "thermal efficiency of the heater";

  ThermoPower.Gas.FlangeA flangeA(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
  ThermoPower.Gas.FlangeB flangeB(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));
  OpenTEMPEST.Flow.Flow1D2w flow1D2w(
    redeclare package Medium = Medium,
    N=N,
    Ntubes=Ntubes,
    TStartIn=TStartIn,
    TStartOut=TStartOut,
    pStart=pStart,
    xStart=xStart,
    Rh=Dhyd/2,
    omega=2*A/Dhyd,
    L=L,
    dpdxNom=dpdx,
    LUDS=LUDS,
    internalHT=false)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
  ThermoPower.Electrical.Load load(Pnom(displayUnit="kW") = 5000, usePowerInput=
       true) annotation (Placement(transformation(extent={{56,58},{76,38}})));
public
  ThermoPower.Thermal.HeatSource1DFV heatSource1DFV(final Nw=N)
    annotation (Placement(transformation(extent={{-10,8},{10,28}})));
  Modelica.Blocks.Math.Gain efficiency(k=Efficiency) annotation (Placement(
        transformation(
        extent={{-4,-4},{4,4}},
        rotation=-90,
        origin={0,32})));
  Modelica.Blocks.Interfaces.RealInput PowerIn annotation (Placement(
        transformation(
        extent={{-18,18},{18,-18}},
        rotation=-90,
        origin={0,58}), iconTransformation(
        extent={{-10,10},{10,-10}},
        rotation=-90,
        origin={0,38})));
  ThermoPower.Electrical.PowerConnection port annotation (Placement(
        transformation(rotation=0, extent={{40,-50},{60,-30}})));
  ThermoPower.Thermal.DHTVolumes heatloss(final N=N) annotation (Placement(
        transformation(extent={{-4,-40},{6,-30}}), iconTransformation(
          extent={{-8,-44},{6,-30}})));
equation
  connect(PowerIn,efficiency. u) annotation (Line(points={{0,58},{0,36.8},{
          8.88178e-16,36.8}}, color={0,0,127}));
  connect(heatSource1DFV.power,efficiency. y) annotation (Line(points={{0,
          22},{0,27.6},{-8.88178e-16,27.6}}, color={0,0,127}));
  connect(heatSource1DFV.wall, flow1D2w.wall)
    annotation (Line(points={{0,15},{0,4.2}}, color={255,127,0}));
  connect(efficiency.u, load.referencePower) annotation (Line(points={{0,36.8},{
          32,36.8},{32,48},{62.7,48}}, color={0,0,127}));
  connect(flow1D2w.wall2, heatloss) annotation (Line(points={{0,-4.1},{0,-20},{0,
          -35},{1,-35}}, color={255,127,0}));
  connect(flow1D2w.outfl, flangeB)
    annotation (Line(points={{9.4,0},{100,0}}, color={159,159,223}));
  connect(flow1D2w.infl, flangeA) annotation (Line(points={{-9.4,-0.2},{-51.7,-0.2},
          {-51.7,0},{-100,0}}, color={159,159,223}));
  connect(load.port, port) annotation (Line(
      points={{66,39.4},{68,39.4},{68,-40},{50,-40}},
      color={0,0,255},
      thickness=0.5));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                                                                Rectangle(
          extent={{-92,28},{90,-30}},
          lineColor={28,108,200},
          radius=10,
          fillColor={238,46,47},
          fillPattern=FillPattern.Solid)}),                      Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end ElectricHeater1D;

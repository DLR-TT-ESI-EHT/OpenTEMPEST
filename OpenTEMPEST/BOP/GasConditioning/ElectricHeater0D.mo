within OpenTEMPEST.BOP.GasConditioning;
model ElectricHeater0D
  extends ThermoPower.Icons.Gas.Tube;
  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
    "Medium model"
    annotation(choicesAllMatching = true);

  parameter Modelica.SIunits.Power PMax = 5000 "Maximal power";

  parameter Real PI_k(min=0, unit="1") = 20 "Gain of controller";
  parameter Modelica.SIunits.Time PI_Ti(min=Modelica.Constants.small)=300
    "Time constant of controller's Integrator block" annotation (Dialog(enable=
          controllerType == .Modelica.Blocks.Types.SimpleController.PI or
          controllerType == .Modelica.Blocks.Types.SimpleController.PID));
  parameter Modelica.SIunits.HeatCapacity C=7800*500*0.05*0.01 "Heat capacity of element (= cp*m)";
  parameter Modelica.SIunits.Temperature TStart(displayUnit="degC")
    "Temperature of element";
  SimpleGasHeater heater(
    redeclare package medium = Medium,
    l=0.5,
    TStart=1103.15) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={30,40})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
    annotation (Placement(transformation(
        extent={{-6,-6},{6,6}},
        rotation=-90,
        origin={30,64})));
  Modelica.Blocks.Continuous.LimPID PID(
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=PI_k,
    Ti=PI_Ti,
    yMax=PMax,
    yMin=0) annotation (Placement(transformation(extent={{76,76},{66,86}})));
protected
  OpenTEMPEST.Flow.SensGasProperty sensOut(
    mfOutput=false,
    pOutput=false,
    hOutput=false,
    XOutput=false,
    YOutput=false,
    HfOutput=false,
    redeclare package Medium = Medium,
    cpOutput=false,
    lambdaOutput=false,
    etaOutput=false,
    TOutput=true,
    rhoOutput=false)
    annotation (Placement(transformation(extent={{62,34},{82,54}})));
  ThermoPower.Electrical.Load load(Pnom=5000, usePowerInput=true)
    annotation (Placement(transformation(extent={{34,84},{14,104}})));
public
  Modelica.Blocks.Interfaces.RealInput in_TSet(unit="K")
                                               annotation (Placement(
        transformation(rotation=0, extent={{100,94},{120,114}}),
        iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={-50,50})));
  ThermoPower.Gas.FlangeA fl_inlet(redeclare package Medium = Medium)
    annotation (Placement(transformation(rotation=0, extent={{-100,0},{-80,20}}),
        iconTransformation(extent={{-120,-20},{-80,20}})));
  ThermoPower.Electrical.PowerConnection powerConnection annotation (Placement(
        transformation(rotation=0, extent={{0,60},{20,80}}), iconTransformation(
          extent={{-20,40},{20,80}})));
  ThermoPower.Gas.FlangeB fl_outlet(redeclare package Medium =
        Medium) annotation (Placement(transformation(rotation=0,
          extent={{80,-20},{100,0}}), iconTransformation(extent={{80,-20},{120,20}})));
  Modelica.Blocks.Interfaces.RealOutput out_TAct(unit="K") annotation (Placement(
        transformation(rotation=0, extent={{106,50},{126,70}}),
        iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={50,50})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor(C=C, T(
        start=TStart))
    annotation (Placement(transformation(extent={{40,52},{60,72}})));
  ThermoPower.Thermal.HT htAmb annotation (Placement(transformation(extent={{40,2},{
            60,22}}),       iconTransformation(extent={{60,-60},{80,-40}})));
equation
  connect(heater.ht[1], prescribedHeatFlow.port)
    annotation (Line(points={{30,44},{30,58}}, color={191,0,0}));
  connect(heater.outlet, sensOut.inlet)
    annotation (Line(points={{38,40},{66,40}}, color={159,159,223}));
  connect(sensOut.T, PID.u_m) annotation (Line(points={{81,50},{86,50},{86,74},{
          71,74},{71,75}}, color={0,0,127}));
  connect(PID.y, load.referencePower)
    annotation (Line(points={{65.5,81},{27.3,81},{27.3,94}}, color={0,0,127}));
  connect(in_TSet, PID.u_s)
    annotation (Line(points={{110,104},{110,81},{77,81}}, color={0,0,127}));
  connect(fl_inlet, heater.inlet) annotation (Line(points={{-90,10},{-34,10},{-34,
          40},{22,40}}, color={159,159,223}));
  connect(powerConnection, load.port)
    annotation (Line(points={{10,70},{10,102.6},{24,102.6}}, color={0,0,255}));
  connect(fl_outlet, sensOut.outlet) annotation (Line(points={{90,-10},{84,-10},
          {84,40},{78,40}}, color={159,159,223}));
  connect(prescribedHeatFlow.Q_flow, load.referencePower) annotation (Line(
        points={{30,70},{30,95},{27.3,95},{27.3,94}}, color={0,0,127}));
  connect(out_TAct, PID.u_m) annotation (Line(points={{116,60},{102,60},{102,66},
          {86,66},{86,74},{71,74},{71,75}}, color={0,0,127}));
  connect(prescribedHeatFlow.port, heatCapacitor.port) annotation (Line(points={
          {30,58},{40,58},{40,52},{50,52}}, color={191,0,0}));
  connect(heatCapacitor.port, htAmb)
    annotation (Line(points={{50,52},{50,52},{50,12}}, color={191,0,0}));
  annotation (Diagram(coordinateSystem(extent={{-100,-100},{100,100}})),
                                                                       Icon(
        coordinateSystem(extent={{-100,-100},{100,100}}), graphics={
                   Text(extent={{-100,-40},{100,-80}}, textString="%name")}),
    Documentation(revisions="<html>
<ul>
<li>12.08.22 by Marius Tomberg:<br>Heat capacity and heat loss port added</li>
</ul>
</html>"));
end ElectricHeater0D;

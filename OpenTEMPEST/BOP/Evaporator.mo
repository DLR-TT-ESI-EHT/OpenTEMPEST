within OpenTEMPEST.BOP;
model Evaporator
 extends Modelica.Icons.UnderConstruction; // check media

  replaceable package gasPhaseMedium =
      Modelica.Media.Interfaces.PartialMedium                          annotation(choicesAllMatching = true);

  parameter Modelica.SIunits.Volume Vd "Drum internal volume";
  parameter Modelica.SIunits.Mass Mm "Drum metal mass";
  parameter Modelica.SIunits.SpecificHeatCapacity cm
    "Specific heat capacity of the drum's metal";

  parameter ThermoPower.Units.HydraulicConductance Kv "Nominal hydraulic conductance of the control valve";

  parameter Modelica.SIunits.Pressure pStart "Pressure start value"
    annotation (Dialog(tab="Initialisation"));
  parameter Modelica.SIunits.Volume VlStart "Start value of drum water volume"
    annotation (Dialog(tab="Initialisation"));

  ThermoPower.Water.DrumEquilibrium drumEquilibrium(
    Vd=Vd,
    Mm=Mm,
    cm=cm,
    pstart=pStart,
    Vlstart=VlStart)
    annotation (Placement(transformation(extent={{-76,-10},{-56,10}})));
  ThermoPower.Water.FlangeA inlet
    annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
  ThermoPower.Water.SinkPressure sinkPressure(use_in_p0=true)
    annotation (Placement(transformation(extent={{10,-10},{30,10}})));
  ThermoPower.Gas.SourceMassFlow sourceMassFlow(
    redeclare package Medium = gasPhaseMedium,
    Xnom={1e-5,1e-5,1e-5,1e-5,1,1e-5},          use_in_w0=true, use_in_T=true)
    annotation (Placement(transformation(extent={{38,10},{58,-10}})));
  ThermoPower.Gas.FlangeB outlet(redeclare package Medium = gasPhaseMedium)
    annotation (Placement(transformation(extent={{80,-10},{100,10}})));
  ThermoPower.Gas.SensP sensP(redeclare package Medium = gasPhaseMedium)
    annotation (Placement(transformation(extent={{80,-6},{60,14}})));
  ThermoPower.Water.ValveLin valveLin(Kv=Kv)
    annotation (Placement(transformation(extent={{-54,-10},{-34,10}})));
  ThermoPower.Water.SensT sensT
    annotation (Placement(transformation(extent={{-12,6},{8,-14}})));
  ThermoPower.Water.SensW sensW
    annotation (Placement(transformation(extent={{-32,6},{-12,-14}})));
  ThermoPower.Thermal.HT hT
    annotation (Placement(transformation(extent={{-12,80},{8,100}}),
        iconTransformation(extent={{-10,40},{10,60}})));
  Modelica.Blocks.Interfaces.RealInput controlValvePosition
    annotation (Placement(transformation(extent={{20,20},{60,60}}),
        iconTransformation(extent={{40,40},{60,60}})));
equation
  connect(drumEquilibrium.feed, inlet) annotation (Line(points={{-75,-4.4},{-75.5,
          -4.4},{-75.5,0},{-90,0}}, color={0,0,255}));
  connect(sourceMassFlow.flange, outlet)
    annotation (Line(points={{58,0},{90,0}}, color={159,159,223}));
  connect(sourceMassFlow.flange, sensP.flange)
    annotation (Line(points={{58,0},{70,0}}, color={159,159,223}));
  connect(sensP.p, sinkPressure.in_p0)
    annotation (Line(points={{63,10},{16,10},{16,8.4}}, color={0,0,127}));
  connect(valveLin.inlet, drumEquilibrium.steam)
    annotation (Line(points={{-54,0},{-54,7.2},{-59.2,7.2}}, color={0,0,255}));
  connect(sensT.outlet, sinkPressure.flange)
    annotation (Line(points={{4,0},{10,0}}, color={0,0,255}));
  connect(sensT.T, sourceMassFlow.in_T)
    annotation (Line(points={{6,-10},{48,-10},{48,-5}}, color={0,0,127}));
  connect(valveLin.outlet, sensW.inlet)
    annotation (Line(points={{-34,0},{-28,0}}, color={0,0,255}));
  connect(sensW.outlet, sensT.inlet)
    annotation (Line(points={{-16,0},{-8,0}}, color={0,0,255}));
  connect(sensW.w, sourceMassFlow.in_w0) annotation (Line(points={{-14,-10},{-10,
          -10},{-10,-16},{42,-16},{42,-5}}, color={0,0,127}));
  connect(hT, drumEquilibrium.wall)
    annotation (Line(points={{-2,90},{-2,22},{-66,22},{-66,-9}},
                                                       color={191,0,0}));
  connect(valveLin.cmd, controlValvePosition) annotation (Line(points={{-44,8},{
          -46,8},{-46,40},{40,40}},  color={0,0,127}));
  annotation (                                                             Icon(
        graphics={Rectangle(
          extent={{-80,40},{80,-40}},
          lineColor={28,108,200},
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid)}));
end Evaporator;

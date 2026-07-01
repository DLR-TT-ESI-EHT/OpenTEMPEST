within OpenTEMPEST.BOP.GasConditioning;
model CoolingCycle "Heat exchanger with integrated water cycle."
  import SI = Modelica.SIunits;
  import Medi = Modelica.Media.Interfaces.PartialMedium;

  replaceable package MediumHeating = Modelica.Media.Interfaces.PartialMedium
    annotation(choicesAllMatching = true);

  replaceable package MediumCooling = Modelica.Media.Interfaces.PartialMedium
   annotation(choicesAllMatching = true);

  parameter SI.Volume VExpansionTank = 1 "Inner volume";
  parameter SI.ThermalConductance UACooling = 63800/10 "Termal conductence";
  parameter SI.Distance LHeat = 1 "Heat tube lenght";
  parameter SI.Distance DHeat = 0.05 "Heat tube diameter";
  parameter SI.Length omegaHeat = 0.157
    "Perimeter of heat transfer surface (single heat tube)";
  parameter Medi.MassFlowRate wnomHeat = 0.1 "Heat Nominal mass flowrate (total)";
  parameter Medi.MassFlowRate wnomCooling = 1 "Cooling Nominal mass flowrate (total)";

  replaceable function flowCharacteristicPump =
      ThermoPower.Functions.PumpCharacteristics.baseFlow
    "Head vs. q_flow characteristic at nominal speed and density"
    annotation (choicesAllMatching=true);

  ThermoPower.Water.Pump pump(
    redeclare package Medium = Modelica.Media.Water.WaterIF97_pT,
    redeclare function flowCharacteristic = flowCharacteristicPump,
    usePowerCharacteristic=false,
    redeclare function efficiencyCharacteristic =
        ThermoPower.Functions.PumpCharacteristics.constantEfficiency (eta_nom=0.8),
    Np0=1,
    use_in_Np=false,
    rho0=1000,
    n0=500,
    CheckValve=false,
    wstart=1,
    w0=1,
    dp0(displayUnit="bar") = 100000,
    use_in_n=true) annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=90,
        origin={50,-4})));

  ThermoPower.Water.Flow1DFV flow1DFV(
    redeclare package Medium = Modelica.Media.Water.WaterIF97_pT,
    N=5,
    L=1,
    H=0,
    A=0.002,
    omega=0.157,
    wnom=1,
    dpnom=100000,
    FluidPhaseStart=ThermoPower.Choices.FluidPhase.FluidPhases.Liquid,
    pstart=2000000,
    initOpt=ThermoPower.Choices.Init.Options.noInit,
    noInitialPressure=true)
    annotation (Placement(transformation(extent={{14,18},{-6,38}})));
  ThermoPower.Water.Flow1DFV flow1DFV1(
    redeclare package Medium = Modelica.Media.Water.WaterIF97_pT,
    N=5,
    L=1,
    H=0,
    A=0.002,
    omega=0.157,
    wnom=1,
    dpnom=100000,
    HydraulicCapacitance=ThermoPower.Choices.Flow1D.HCtypes.Middle,
    FluidPhaseStart=ThermoPower.Choices.FluidPhase.FluidPhases.Liquid,
    pstart=2000000,
    noInitialPressure=true)
    annotation (Placement(transformation(extent={{-4,-24},{16,-44}})));
  ThermoPower.Water.SensT SensT_A_in(redeclare package Medium =
        Modelica.Media.Water.WaterIF97_pT) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={-48,12})));
  ThermoPower.Thermal.MetalWallFV metalWallFV(
    Nw=4,
    M=2,
    cm=900,
    WallRes=true,
    UA_ext=UACooling)
    annotation (Placement(transformation(extent={{-4,-66},{16,-46}})));
  ThermoPower.Thermal.MetalWallFV metalWallFV1(
    Nw=4,
    M=2,
    cm=900,
    WallRes=true,
    UA_ext=15000/10)
    annotation (Placement(transformation(extent={{-6,56},{14,36}})));
  ThermoPower.Water.Flow1DFV flow1DFV2(
    redeclare package Medium = Modelica.Media.Water.WaterIF97_pT,
    N=5,
    L=1,
    H=0,
    A=0.196,
    omega=1.57,
    wnom=wnomCooling,
    FluidPhaseStart=ThermoPower.Choices.FluidPhase.FluidPhases.Liquid,
    pstart=100000,
    noInitialPressure=true,
    redeclare model HeatTransfer =
        ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer)
    annotation (Placement(transformation(extent={{18,-90},{-2,-70}})));
  OpenTEMPEST.BOP.GasConditioning.ExpansionTank expansionTank(
    redeclare package Medium = Modelica.Media.Water.WaterIF97_pT,
    V=VExpansionTank,
    S=1,
    FluidPhaseStart=ThermoPower.Choices.FluidPhase.FluidPhases.Liquid,
    pstart=100000,
    gasShareStart=0.9)
    annotation (Placement(transformation(extent={{28,-44},{48,-24}})));
  ThermoPower.Gas.Flow1DFV flow1DFV3(
    redeclare package Medium = MediumHeating,
    N=5,
    L=LHeat,
    A=AHeat,
    omega=0.157,
    Dhyd=0.05,
    wnom=wnomHeat,
    noInitialPressure=true,
    redeclare model HeatTransfer =
        ThermoPower.Thermal.HeatTransferFV.ConstantThermalConductance (UA=100))
    annotation (Placement(transformation(extent={{-6,82},{14,62}})));
  Modelica.Blocks.Continuous.LimPID PID(
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=-11.59,
    Ti=8,
    yMax=180,
    yMin=40)
    annotation (Placement(transformation(extent={{-70,-40},{-50,-60}})));
  Modelica.Blocks.Interfaces.RealInput u_s annotation (Placement(transformation(
          rotation=0, extent={{-120,-50},{-100,-30}})));
  ThermoPower.Water.FlangeA waterInlet(redeclare package Medium = MediumCooling)
                                           annotation (Placement(transformation(
          rotation=0, extent={{90,-108},{110,-88}}),  iconTransformation(extent={{50,-148},
            {110,-88}})));
  ThermoPower.Water.FlangeB waterOutlet(redeclare package Medium =
        MediumCooling)                     annotation (Placement(transformation(
          rotation=0, extent={{-110,-112},{-90,-92}}), iconTransformation(
          extent={{-144,-146},{-90,-92}})));
  ThermoPower.Gas.FlangeA gasInlet(redeclare package Medium = MediumHeating)
    annotation (Placement(transformation(rotation=0, extent={{-110,90},{-90,110}}),
        iconTransformation(extent={{-144,92},{-90,146}})));
  ThermoPower.Gas.FlangeB gasOutlet(redeclare package Medium = MediumHeating)
                           annotation (Placement(transformation(rotation=0,
          extent={{88,88},{108,108}}), iconTransformation(extent={{94,90},{150,146}})));
  Modelica.Blocks.Interfaces.RealOutput TR
    annotation (Placement(transformation(extent={{-110,6},{-90,26}})));

protected
  parameter SI.Area AHeat = 3.14159*DHeat^2/4;

equation

  // AHeat = 3.14159*DHeat^2/4;
  connect(pump.outfl,flow1DFV. infl) annotation (Line(points={{57,2},{54,2},{54,
          28},{14,28}},     color={0,0,255}));
  connect(flow1DFV.outfl,SensT_A_in. inlet)
    annotation (Line(points={{-6,28},{-44,28},{-44,18}}, color={0,0,255}));
  connect(metalWallFV.int,flow1DFV1. wall) annotation (Line(points={{6,-53},{6,-39}},
                                     color={255,127,0}));
  connect(metalWallFV1.int,flow1DFV. wall)
    annotation (Line(points={{4,43},{4,33},{4,33}}, color={255,127,0}));
  connect(flow1DFV2.wall,metalWallFV. ext)
    annotation (Line(points={{8,-75},{8,-59.1},{6,-59.1}},
                                                     color={255,127,0}));
  connect(SensT_A_in.outlet,flow1DFV1. infl)
    annotation (Line(points={{-44,6},{-44,-34},{-4,-34}}, color={0,0,255}));
  connect(flow1DFV1.outfl,expansionTank. inlet)
    annotation (Line(points={{16,-34},{27.9,-34}}, color={0,0,255}));
  connect(expansionTank.outlet,pump. infl)
    annotation (Line(points={{48,-34},{52,-34},{52,-12}}, color={0,0,255}));
  connect(flow1DFV3.wall,metalWallFV1. ext)
    annotation (Line(points={{4,67},{4,49.1}}, color={255,127,0}));
  connect(SensT_A_in.T,PID. u_m) annotation (Line(points={{-54,4},{-58,4},{-58,-38},
          {-60,-38}},      color={0,0,127}));
  connect(u_s, PID.u_s) annotation (Line(points={{-110,-40},{-86,-40},{-86,-50},
          {-72,-50}}, color={0,0,127}));
  connect(waterInlet, flow1DFV2.infl) annotation (Line(points={{100,-98},{32,
          -98},{32,-80},{18,-80}},
                              color={0,0,255}));
  connect(waterOutlet, flow1DFV2.outfl) annotation (Line(points={{-100,-102},{-38,
          -102},{-38,-80},{-2,-80}}, color={0,0,255}));
  connect(gasInlet, flow1DFV3.infl) annotation (Line(points={{-100,100},{-54,
          100},{-54,72},{-6,72}},
                            color={159,159,223}));
  connect(gasOutlet, flow1DFV3.outfl) annotation (Line(points={{98,98},{18,98},{
          18,72},{14,72}}, color={159,159,223}));
  connect(SensT_A_in.T, TR) annotation (Line(points={{-54,4},{-79,4},{-79,16},{-100,
          16}}, color={0,0,127}));
  connect(gasOutlet, gasOutlet)
    annotation (Line(points={{98,98},{98,98}}, color={159,159,223}));
  connect(waterInlet, waterInlet) annotation (Line(points={{100,-98},{100,-98}},
                               color={0,0,255}));
  connect(PID.y, pump.in_n) annotation (Line(points={{-49,-50},{70,-50},{70,-6.6},
          {58,-6.6}}, color={0,0,127}));
  connect(gasInlet, gasInlet)
    annotation (Line(points={{-100,100},{-100,100}}, color={159,159,223}));
  annotation (Icon(graphics={Rectangle(extent={{-100,100},{100,-100}},
                                                                   lineColor={28,
              108,200}), Text(
          extent={{-8,12},{6,8}},
          lineColor={28,108,200},
          textString="fdh"),
        Rectangle(
          extent={{-94,80},{92,72}},
          lineColor={255,0,0},
          fillColor={255,0,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-92,-84},{94,-92}},
          lineColor={0,0,255},
          fillColor={0,0,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-86,56},{90,-72}},
          lineColor={0,255,255},
          fillColor={85,255,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-70,38},{72,-52}},
          lineColor={0,255,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-42,14},{48,-30}},
          lineColor={0,0,255},
          fillColor={0,0,255},
          fillPattern=FillPattern.None,
          textStyle={TextStyle.Bold,TextStyle.Italic,TextStyle.UnderLine},
          textString="Cooling cycle")}));
end CoolingCycle;

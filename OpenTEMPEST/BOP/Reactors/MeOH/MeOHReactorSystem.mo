within OpenTEMPEST.BOP.Reactors.MeOH;
model MeOHReactorSystem
  "System representing a MeOH reactor process with separation, compression and heat transfer for heat balancing purposes"
  extends Modelica.Icons.UnderConstruction;

  import SI = Modelica.SIunits;
   parameter Real conversion = 0.9 "ratio C/H>=0.25; const conversion rate nMeOH/(nCO+nCO2)in"
   annotation (Dialog(tab="MeOH Reactor"));
  parameter SI.Temperature TReactor=260 + 273.15
    "[K] temperature of Purge gas (unequal reactor temperature)"
    annotation (Dialog(tab="MeOH Reactor"));
  parameter SI.Temperature TPurge=TReactor
    "[K] temperature of Purge gas (unequal reactor temperature)"
    annotation (Dialog(tab="MeOH Reactor"));
  parameter SI.Temperature TProduct=40 + 273.15
    "[K] temperature of condensated product flow"
    annotation (Dialog(tab="MeOH Reactor"));
  parameter SI.Pressure pReactor=5000000
    "pressure gas Purge after flush/sparator (unequal reator pressure)"
    annotation (Dialog(tab="MeOH Reactor"));
  // parameter Real fac = 1 "factor for accoutning for additional heat sinks"
  // annotation (Dialog(tab="MeOH Reactor"));
  parameter SI.SpecificEnthalpy dhvMeOH=1165000
  annotation (Dialog(tab="MeOH Reactor"));

  MeOH.MeOHEnergyBalance meOHEnergyBalance_noSep(
    conversion=conversion,
    TReactor=TReactor,
    TPurge=TPurge,
    TProduct=TProduct,
    pReactor=pReactor,
    dhvMeOH=dhvMeOH)
    annotation (Placement(transformation(extent={{62,12},{82,32}})));
  Compressor compressor2(redeclare package Medium =
        OpenTEMPEST.Medium.Fuel_MethanolReactor)
    annotation (Placement(transformation(extent={{16,32},{36,12}})));
  ThermoPower.Electrical.Grid grid(Pgrid=999999999) annotation (Placement(
        transformation(
        extent={{7,-7},{-7,7}},
        rotation=90,
        origin={19,-1})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T=system.T_amb)
                      annotation (Placement(transformation(
        extent={{-6,6},{6,-6}},
        rotation=90,
        origin={34,-2})));
  ThermoPower.Gas.SensT sensTInter(redeclare package Medium =
        OpenTEMPEST.Medium.Fuel_MethanolReactor)
    annotation (Placement(transformation(extent={{-16,-26},{-2,-12}})));
  ThermoPower.Gas.SensP sensPInter(redeclare package Medium =
        OpenTEMPEST.Medium.Fuel_MethanolReactor)
    annotation (Placement(transformation(extent={{16,-26},{2,-12}})));
  OpenTEMPEST.Flow.SensGasProperty sensGasHighP(
    mfOutput=false,
    pOutput=false,
    HfOutput=false,
    redeclare package Medium = OpenTEMPEST.Medium.Fuel_MethanolReactor,
    cpOutput=false,
    etaOutput=false,
    rhoOutput=false)
    annotation (Placement(transformation(extent={{36,18},{50,32}})));
  ThermoPower.Gas.SensP sensPHigh(redeclare package Medium =
        OpenTEMPEST.Medium.Fuel_MethanolReactor)
    annotation (Placement(transformation(extent={{66,30},{52,44}})));
  Modelica.Blocks.Math.Gain gain(k=1/25e5)
    "write desired p_in in denominator of fraction"
    annotation (Placement(transformation(extent={{36,36},{28,44}})));
  Compressor compressor1(eta=0.9, redeclare package Medium =
        OpenTEMPEST.Medium.Fuel_MethanolReactor)
    annotation (Placement(transformation(extent={{-38,-12},{-18,-32}})));
  ThermoPower.Electrical.Grid grid1(Pgrid=999999999)
                                                    annotation (Placement(
        transformation(
        extent={{7,-7},{-7,7}},
        rotation=90,
        origin={-37,-43})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature1(T=system.T_amb)
                      annotation (Placement(transformation(
        extent={{-6,6},{6,-6}},
        rotation=90,
        origin={-20,-42})));
  Modelica.Blocks.Math.Gain gain1(k=1/1e5)
    "write desired p_in in denominator of fraction"
    annotation (Placement(transformation(extent={{-10,-8},{-18,0}})));
  OpenTEMPEST.BOP.GasConditioning.SimpleGasHeater simpleGasInterCooler(
    redeclare package medium = OpenTEMPEST.Medium.Fuel_MethanolReactor,
    l=1,
    h=0.2,
    w=0.5,
    xStart=OpenTEMPEST.Medium.Fuel_MethanolReactor.reference_X,
    nParallel=60,
    useAlphaIn=false) annotation (Placement(transformation(
        extent={{8,8},{-8,-8}},
        rotation=180,
        origin={4,14})));
  OpenTEMPEST.BOP.GasConditioning.SimpleGasHeater simpleGasCooler(
    redeclare package medium = OpenTEMPEST.Medium.Fuel_MethanolReactor,
    l=3,
    h=0.2,
    w=1,
    xStart=OpenTEMPEST.Medium.Fuel_MethanolReactor.reference_X,
    nParallel=70,
    useAlphaIn=false) annotation (Placement(transformation(
        extent={{-8,8},{8,-8}},
        rotation=-90,
        origin={-76,14})));
  LiquidSeparator liquidSeparator1 annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-62,-36})));
  inner ThermoPower.System system
    annotation (Placement(transformation(extent={{-100,80},{-80,100}})));
  Modelica.Blocks.Interfaces.RealOutput mfMeOH annotation (Placement(
        transformation(extent={{94,-40},{114,-20}}), iconTransformation(extent={{94,-62},
            {120,-36}})));
  ThermoPower.Gas.SensT sensTSeperator(redeclare package Medium =
        OpenTEMPEST.Medium.Fuel_MethanolReactor) annotation (Placement(
        transformation(
        extent={{-7,7},{7,-7}},
        rotation=-90,
        origin={-79,-11})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
    prescribedTemperature annotation (Placement(transformation(
        extent={{-6,-6},{6,6}},
        rotation=0,
        origin={-108,14})));
  Modelica.Blocks.Continuous.LimPID PID(
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=15.5,
    Ti=3.3*0.035,
    yMax=100 + 273.15,
    yMin=10 + 273.15)
            annotation (Placement(transformation(extent={{6,-6},{-6,6}},
        rotation=0,
        origin={-114,-14})));
  Modelica.Blocks.Sources.Constant const1(k=35 + 273.15)
    annotation (Placement(transformation(extent={{-6,-6},{6,6}},
        rotation=180,
        origin={-98,-14})));
  OpenTEMPEST.Heat.InsulationVarArea metalWall(
    lambda=50,
    c=500,
    rho(displayUnit="kg/m3") = 7850,
    area1=3,
    area2=3,
    depth=1.5/1000,
    TStart=298.15) annotation (Placement(transformation(
        extent={{-6.5,-5.5},{6.5,5.5}},
        rotation=90,
        origin={-86.5,13.5})));
  OpenTEMPEST.Heat.InsulationVarArea metalWall1(
    lambda=50,
    c=500,
    rho(displayUnit="kg/m3") = 7850,
    area1=3,
    area2=3,
    depth=2/1000,
    TStart=298.15) annotation (Placement(transformation(
        extent={{-6.5,-5.5},{6.5,5.5}},
        rotation=90,
        origin={-94.5,13.5})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
    prescribedTemperature1
                          annotation (Placement(transformation(
        extent={{-6,-6},{6,6}},
        rotation=-90,
        origin={4,50})));
  OpenTEMPEST.Heat.InsulationVarArea metalWall2(
    lambda=50,
    c=500,
    rho(displayUnit="kg/m3") = 7850,
    area1=3,
    area2=3,
    depth=2/1000,
    TStart=298.15) annotation (Placement(transformation(
        extent={{-6.5,-5.5},{6.5,5.5}},
        rotation=0,
        origin={3.5,33.5})));
  OpenTEMPEST.Heat.InsulationVarArea metalWall3(
    lambda=50,
    c=500,
    rho(displayUnit="kg/m3") = 7850,
    area1=3,
    area2=3,
    depth=1.5/1000,
    TStart=298.15) annotation (Placement(transformation(
        extent={{-6.5,-5.5},{6.5,5.5}},
        rotation=0,
        origin={3.5,25.5})));
  Modelica.Blocks.Continuous.LimPID PID1(
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=13.5,
    Ti=3.3*0.04,
    yMax=120 + 273.15,
    yMin=10 + 273.15)
            annotation (Placement(transformation(extent={{6,-6},{-6,6}},
        rotation=0,
        origin={28,62})));
  Modelica.Blocks.Sources.Constant const2(k=225 + 273.15)
    annotation (Placement(transformation(extent={{56,56},{44,68}})));
  ThermoPower.Thermal.HT hT annotation (Placement(transformation(extent={{-12,
            94},{12,118}}), iconTransformation(extent={{-10,38},{10,58}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
    annotation (Placement(transformation(extent={{-22,76},{-2,96}})));
  Modelica.Blocks.Sources.RealExpression realExpression(y=-simpleGasCooler.ht[1].Q_flow
         - simpleGasInterCooler.ht[1].Q_flow - meOHEnergyBalance_noSep.hT.Q_flow)
    annotation (Placement(transformation(extent={{-48,78},{-30,94}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1(alpha=0.1)
    annotation (Placement(transformation(extent={{86,40},{74,52}})));
  Modelica.Blocks.Sources.RealExpression realExpression1(y=-
        meOHEnergyBalance_noSep.hT.Q_flow)
    annotation (Placement(transformation(extent={{104,38},{92,54}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow2(alpha=0.1)
    annotation (Placement(transformation(extent={{6,6},{-6,-6}},
        rotation=90,
        origin={-60,-14})));
  Modelica.Blocks.Sources.RealExpression realExpression2(y=-liquidSeparator1.hT.Q_flow)
    annotation (Placement(transformation(extent={{6,8},{-6,-8}},
        rotation=90,
        origin={-60,6})));
  ThermoPower.Gas.FlangeA infl(redeclare package Medium =
        OpenTEMPEST.Medium.Fuel_MethanolReactor) annotation (Placement(
        transformation(extent={{-114,38},{-94,58}}), iconTransformation(extent=
            {{-122,-6},{-94,22}})));
  ThermoPower.Water.FlangeB outflWater(redeclare package Medium =
        ThermoPower.Water.StandardWater) annotation (Placement(transformation(
          extent={{-66,-102},{-46,-82}}),iconTransformation(extent={{-74,-110},
            {-46,-82}})));
  ThermoPower.Gas.FlangeB outflPurge(redeclare package Medium =
        OpenTEMPEST.Medium.Fuel_MethanolReactor) annotation (Placement(
        transformation(extent={{102,2},{122,22}}), iconTransformation(extent={{
            94,-6},{122,22}})));
equation
  connect(grid.port,compressor2. powerConnection) annotation (Line(
      points={{19,5.02},{19,15},{23,15}},
      color={0,0,255},
      thickness=0.5));
  connect(fixedTemperature.port,compressor2. heatLosses)
    annotation (Line(points={{34,4},{34,15},{29,15}},  color={191,0,0}));
  connect(sensGasHighP.outlet, meOHEnergyBalance_noSep.infl) annotation (Line(
        points={{47.2,22.2},{56.9,22.2},{56.9,20.2222},{60.6,20.2222}},
                                                              color={159,159,
          223}));
  connect(sensGasHighP.outlet, sensPHigh.flange) annotation (Line(points={{47.2,
          22.2},{56,22.2},{56,34.2},{59,34.2}}, color={159,159,223}));
  connect(sensPHigh.p,gain. u) annotation (Line(points={{54.1,41.2},{54.1,40},{36.8,
          40}},     color={0,0,127}));
  connect(gain.y,compressor2. Pi)
    annotation (Line(points={{27.6,40},{26,40},{26,30}},
                                                      color={0,0,127}));
  connect(grid1.port,compressor1. powerConnection) annotation (Line(
      points={{-37,-36.98},{-37,-29},{-31,-29}},
      color={0,0,255},
      thickness=0.5));
  connect(fixedTemperature1.port,compressor1. heatLosses) annotation (Line(
        points={{-20,-36},{-20,-29},{-25,-29}},
                                              color={191,0,0}));
  connect(gain1.y,compressor1. Pi) annotation (Line(points={{-18.4,-4},{-28,-4},
          {-28,-14}},   color={0,0,127}));
  connect(sensTInter.outlet,sensPInter. flange)
    annotation (Line(points={{-4.8,-21.8},{9,-21.8}},color={159,159,223}));
  connect(sensPInter.p,gain1. u) annotation (Line(points={{4.1,-14.8},{4.1,-4},{
          -9.2,-4}},                    color={0,0,127}));
  connect(compressor1.outlet,sensTInter. inlet) annotation (Line(points={{-19,-22},
          {-16,-22},{-16,-21.8},{-13.2,-21.8}},  color={159,159,223}));
  connect(sensPInter.flange,simpleGasInterCooler. inlet) annotation (Line(
        points={{9,-21.8},{-2.4,-21.8},{-2.4,14}},   color={159,159,223}));
  connect(compressor2.outlet, sensGasHighP.inlet) annotation (Line(points={{35,
          22},{38,22},{38,22.2},{38.8,22.2}}, color={159,159,223}));
  connect(simpleGasInterCooler.outlet,compressor2. inlet) annotation (Line(
        points={{10.4,14},{14,14},{14,22},{17,22}},
                                            color={159,159,223}));
  connect(liquidSeparator1.outflGas,compressor1. inlet) annotation (Line(
        points={{-54,-33},{-42,-33},{-42,-22},{-37,-22}},           color={
          159,159,223}));
  connect(meOHEnergyBalance_noSep.outMfProduct, mfMeOH) annotation (Line(points={{83.4,
          20.2222},{90,20.2222},{90,-30},{104,-30}},
                                                 color={0,0,127}));
  connect(simpleGasCooler.outlet, sensTSeperator.inlet) annotation (Line(points=
         {{-76,7.6},{-76,2.4},{-76.2,2.4},{-76.2,-6.8}}, color={159,159,223}));
  connect(sensTSeperator.outlet, liquidSeparator1.infl) annotation (Line(points=
         {{-76.2,-15.2},{-76.2,-36},{-70,-36}}, color={159,159,223}));
  connect(const1.y, PID.u_s)
    annotation (Line(points={{-104.6,-14},{-106.8,-14}}, color={0,0,127}));
  connect(metalWall.ht2, simpleGasCooler.ht[1]) annotation (Line(points={{-83.75,
          13.5},{-83.75,14},{-79.2,14}},                       color={191,0,0}));
  connect(prescribedTemperature.port, metalWall1.ht1) annotation (Line(points={{-102,14},
          {-102,13.37},{-97.25,13.37}},        color={191,0,0}));
  connect(metalWall1.ht2, metalWall.ht1) annotation (Line(points={{-91.75,13.5},
          {-90,13.5},{-90,13.37},{-89.25,13.37}},      color={191,0,0}));
  connect(PID.y, prescribedTemperature.T) annotation (Line(points={{-120.6,-14},
          {-120.6,14},{-115.2,14}},        color={0,0,127}));
  connect(prescribedTemperature1.port, metalWall2.ht1) annotation (Line(points={{4,44},{
          4,36.25},{3.37,36.25}},                           color={191,0,0}));
  connect(PID1.y, prescribedTemperature1.T) annotation (Line(points={{21.4,62},{
          4,62},{4,57.2}},       color={0,0,127}));
  connect(const2.y, PID1.u_s) annotation (Line(points={{43.4,62},{35.2,62}},
                           color={0,0,127}));
  connect(sensGasHighP.T, PID1.u_m) annotation (Line(points={{49.3,29.2},{49.3,
          30},{48,30},{48,48},{28,48},{28,54.8}}, color={0,0,127}));
  connect(metalWall2.ht2, metalWall3.ht1) annotation (Line(points={{3.5,30.75},{
          2.75,30.75},{2.75,28.25},{3.37,28.25}},  color={191,0,0}));
  connect(metalWall3.ht2, simpleGasInterCooler.ht[1]) annotation (Line(points={{3.5,
          22.75},{3.5,18},{4,18},{4,17.2}},      color={191,0,0}));
  connect(prescribedHeatFlow.port, hT)
    annotation (Line(points={{-2,86},{0,86},{0,106}},color={191,0,0}));
  connect(realExpression.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-29.1,
          86},{-22,86}},                     color={0,0,127}));
  connect(prescribedHeatFlow1.Q_flow, realExpression1.y)
    annotation (Line(points={{86,46},{91.4,46}}, color={0,0,127}));
  connect(prescribedHeatFlow1.port, meOHEnergyBalance_noSep.hT)
    annotation (Line(points={{74,46},{72,46},{72,28.2222}},
                                                         color={191,0,0}));
  connect(prescribedHeatFlow2.Q_flow, realExpression2.y) annotation (Line(
        points={{-60,-8},{-59.4,-8},{-59.4,-0.6},{-60,-0.6}}, color={0,0,127}));
  connect(prescribedHeatFlow2.port, liquidSeparator1.hT) annotation (Line(
        points={{-60,-20},{-62,-20},{-62,-29.8}}, color={191,0,0}));
  connect(sensTSeperator.T, PID.u_m) annotation (Line(points={{-83.2,-15.9},{-83.2,
          -26},{-114,-26},{-114,-21.2}}, color={0,0,127}));
  connect(outflPurge, outflPurge) annotation (Line(points={{112,12},{112,12},
          {112,12}},
                color={159,159,223}));
  connect(meOHEnergyBalance_noSep.outflPurge, outflPurge) annotation (Line(
        points={{83.4,24.6667},{94,24.6667},{94,12},{112,12}},
                                                     color={159,159,223}));
  connect(simpleGasCooler.inlet, infl) annotation (Line(points={{-76,20.4},{-78,
          20.4},{-78,48},{-104,48}}, color={159,159,223}));
  connect(liquidSeparator1.outflLiquid, outflWater)
    annotation (Line(points={{-54,-39},{-56,-39},{-56,-92}}, color={0,0,255}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,80}}),                                   graphics={
          Rectangle(
          extent={{-94,38},{94,-82}},
          lineColor={28,108,200},
          fillColor={79,238,238},
          fillPattern=FillPattern.HorizontalCylinder), Text(
          extent={{-68,0},{60,-40}},
          lineColor={28,108,200},
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={79,238,238},
          textString="MeOH Synthesis 
Process")}),                                                     Diagram(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,80}})),
    experiment(StopTime=10000, __Dymola_Algorithm="Cvode"),
    __Dymola_experimentSetupOutput(equidistant=false),
    __Dymola_experimentFlags(
      Advanced(GenerateVariableDependencies=false, OutputModelicaCode=true),
      Evaluate=false,
      OutputCPUtime=true,
      OutputFlatModelica=false));
end MeOHReactorSystem;

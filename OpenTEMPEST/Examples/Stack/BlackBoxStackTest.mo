within OpenTEMPEST.Examples.Stack;
model BlackBoxStackTest "Simple example model using the model blackBoxStack, adiabatic"
  extends Modelica.Icons.Example;

  import SI = Modelica.SIunits;

  // Initial conditions
  parameter SI.Pressure pStart=100000;
  parameter SI.Temperature TStart(displayUnit="K")=973.15;

  // Stack
  parameter Integer nCells=24575;
  parameter Integer nParallel = 1;
  parameter SI.Area ACell = 127.8/100^2;
  parameter SI.HeatCapacity C = 1 "Capacity of each stack if nParallel>=1";
  parameter Units.AreaSpecificResistance alphaASR=728;
  parameter SI.LinearTemperatureCoefficient betaASR = -8.29e-3;
  parameter Integer n(min=1) = 10 "Number of axially discretized units in heat ports";
  parameter Integer intermediatePlateDistance = 10 "Number of cells between two intermediate plates, set to >nCells to have none";
  parameter Boolean calcPressureDrop=false   "Define in pressure drop calculation is needed";
  parameter Boolean useDhtInletOutlet = true;
  parameter SI.PerUnit QLoss = 0.05 "Stack heat loss relative to stack DC power output";


  // Components
  SOC.Stack.BlackBoxStack stack(
    nCells=nCells,
    nParallel=nParallel,
    ACell=ACell,
    C=C,
    alphaASR=alphaASR,
    betaASR=betaASR,
    n=n,
    intermediatePlateDistance=intermediatePlateDistance,
    calcPressureDrop=calcPressureDrop,
    useDhtInletOutlet=useDhtInletOutlet,
    redeclare model Electrochem = OpenTEMPEST.SOC.Electrochem.Components.ASR_Steam (redeclare model ASRobj =
            OpenTEMPEST.SOC.Electrochem.ASR.ASR_Exponential (A=alphaASR, B=betaASR)),
    XStartGas=fuelSource.Xnom,
    XStartAir=airSource.Xnom,
    TASR(start=681.85))
    annotation (Placement(transformation(extent={{-28,-26},{28,28}})));
  BOP.PowerElectronics powerElectronics(nChannels=1, currentFactor=1)
                                        annotation (Placement(transformation(extent={{-110,-12},{-84,12}})));
  Flow.SourceMassFlowEC fuelSource(
    p0=pStart,
    T(displayUnit="K") = TStart,
    compositionNom={1,0.0001,0.0001,0.0001,0.0001,0.0001},
    RC0=0.75,
    nCells=nCells,
    minCurrentForDosing=12.78,
    use_in_I=true) annotation (Placement(transformation(extent={{-112,32},{-92,52}})));
  ThermoPower.Gas.SourceMassFlow airSource(
    redeclare package Medium = OpenTEMPEST.Medium.Air_Medium,
    w0=2.667,
    p0=pStart,
    T(displayUnit="K") = TStart)
           annotation (Placement(transformation(extent={{-112,-64},{-92,-44}},
          rotation=0)));

  Modelica.Blocks.Sources.Constant currentDensitySOFC(k=0.4) annotation (Placement(transformation(extent={{-216,32},{-196,52}})));

  ThermoPower.Gas.SinkPressure sinkAir(
    redeclare package Medium = OpenTEMPEST.Medium.Air_Medium,
    p0=pStart,
    use_in_p0=false) annotation (Placement(transformation(extent={{88,-58},{108,-38}},rotation=0)));
  ThermoPower.Gas.SinkPressure sinkFuel(redeclare package Medium = OpenTEMPEST.Medium.Fuel_CH4, p0=pStart)
    annotation (Placement(transformation(extent={{88,32},{108,52}}, rotation=0)));
  ThermoPower.Electrical.Grid grid(Pgrid=55555555555) annotation (Placement(transformation(extent={{-82,-36},{-62,-16}})));
  Modelica.Blocks.Math.Gain gain(k=ACell*1e4)
                                 annotation (Placement(transformation(extent={{-160,36},{-148,48}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow(T_ref(displayUnit="K") = 240.6)
    annotation (Placement(transformation(extent={{-72,62},{-52,82}})));
  Modelica.Blocks.Sources.RealExpression realExpression(y=-abs(powerElectronics.PDC)*QLoss)
                                                        annotation (Placement(transformation(extent={{-132,64},{-112,84}})));
equation
  // Flanges connections
  connect(airSource.flange, stack.airInlet)
    annotation (Line(points={{-92,-54},{-32,-54},{-32,-9.8},{-22.4,-9.8}}, color={159,159,223}));
  connect(fuelSource.flange, stack.fuelInlet)
    annotation (Line(points={{-92,42},{-32,42},{-32,11.8},{-22.4,11.8}}, color={159,159,223}));
  connect(stack.fuelOutlet, sinkFuel.flange)
    annotation (Line(points={{22.4,11.8},{78,11.8},{78,42},{88,42}}, color={159,159,223}));
  connect(stack.airOutlet, sinkAir.flange)
    annotation (Line(points={{22.4,-9.8},{78,-9.8},{78,-48},{88,-48}}, color={159,159,223}));

  // Electrical connections
  connect(powerElectronics.pinN[1], stack.pinP)
    annotation (Line(points={{-93.1,-10.8},{-93.1,-16},{-34,-16},{-34,-3.05},{-21,-3.05}}, color={0,0,255}));
  connect(powerElectronics.pinP[1], stack.pinN)
    annotation (Line(points={{-100.9,-10.8},{-100.9,-16},{-114,-16},{-114,16},{-34,16},{-34,5.05},{-21,5.05}}, color={0,0,255}));
  connect(powerElectronics.gridConnectionAC, grid.port)
    annotation (Line(
      points={{-85.3,-3.6},{-80,-3.6},{-80,-14},{-84,-14},{-84,-18},{-86,-18},{-86,-26},{-80.6,-26}},
      color={0,0,255},
      thickness=0.5));

  connect(currentDensitySOFC.y, gain.u) annotation (Line(points={{-195,42},{-161.2,42}}, color={0,0,127}));
  connect(gain.y, fuelSource.in_I) annotation (Line(points={{-147.4,42},{-116,42},{-116,47},{-108,47}}, color={0,0,127}));
  connect(gain.y, powerElectronics.currentInput[1])
    annotation (Line(points={{-147.4,42},{-116,42},{-116,18},{-100.9,18},{-100.9,10.8}}, color={0,0,127}));
  connect(prescribedHeatFlow.port, stack.htSideLeft)
    annotation (Line(points={{-52,72},{-11.2,72},{-11.2,18.82}}, color={191,0,0}));
  connect(realExpression.y, prescribedHeatFlow.Q_flow)
    annotation (Line(points={{-111,74},{-80,74},{-80,72},{-72,72}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end BlackBoxStackTest;

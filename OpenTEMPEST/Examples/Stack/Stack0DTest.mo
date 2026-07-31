within OpenTEMPEST.Examples.Stack;
model Stack0DTest "Simple example model using the model Stack0D, adiabatic"
  extends Modelica.Icons.Example;

  import SI = Modelica.SIunits;

  // Initial conditions
  parameter SI.Pressure pStart=100000;
  parameter SI.Temperature TStart(displayUnit="K")=973.15;

  // Stack
  parameter Integer nCells=10;
  parameter Integer nParallel = 1;
  parameter SI.Area ACell = 127.8/100^2;
  parameter SI.HeatCapacity C = 1 "Capacity of each stack if nParallel>=1";
  parameter Integer n(min=1) = 10 "Number of axially discretized units in heat ports";
  parameter Integer intermediatePlateDistance = 10 "Number of cells between two intermediate plates, set to >nCells to have none";
  parameter Boolean calcPressureDrop=false   "Define in pressure drop calculation is needed";
  parameter Boolean useDhtInletOutlet = true;

  // Components
  SOC.Stack.Stack0D stack(
    nCells=nCells,
    nParallel=nParallel,
    ACell=ACell,
    C=C,
    n=n,
    intermediatePlateDistance=intermediatePlateDistance,
    calcPressureDrop=calcPressureDrop,
    useDhtInletOutlet=useDhtInletOutlet,
    redeclare model Electrochem = OpenTEMPEST.SOC.Electrochem.Components.ASR_Steam (redeclare model ASRobj =
            OpenTEMPEST.SOC.Electrochem.ASR.ASR_log (A=stack.alphaASR, B=stack.betaASR)),
    XStartGas=fuelSource.Xnom,
    XStartAir=airSource.Xnom,
    TASR(start=681.85)) annotation (Placement(transformation(extent={{-28,-26},{28,28}})));
  BOP.PowerElectronics powerElectronics(nChannels=1, currentFactor=1)
                                        annotation (Placement(transformation(extent={{-110,-12},{-84,12}})));
  ThermoPower.Gas.SourceMassFlow airSource(
    redeclare package Medium = OpenTEMPEST.Medium.Air_Medium,
    w0=nCells*4.29048E-04,
    p0=1.1*sinkAir.p0,
    T(displayUnit="K") = TStart)
           annotation (Placement(transformation(extent={{-112,-64},{-92,-44}},
          rotation=0)));

  ThermoPower.Gas.SinkPressure sinkAir(
    redeclare package Medium = OpenTEMPEST.Medium.Air_Medium,
    p0=pStart,
    use_in_p0=false) annotation (Placement(transformation(extent={{88,-58},{108,-38}},rotation=0)));
  ThermoPower.Gas.SinkPressure sinkFuel(redeclare package Medium =
        OpenTEMPEST.Medium.Fuel_CH4,                                                            p0=pStart)
    annotation (Placement(transformation(extent={{88,32},{108,52}}, rotation=0)));
  ThermoPower.Electrical.Grid grid(Pgrid=55555555555) annotation (Placement(transformation(extent={{-82,-36},{-62,-16}})));
  Blocks.Sources.MultiRamp currentRamp(
    n=3,
    height=fill(-15, currentRamp.n),
    duration=fill(2000, currentRamp.n),
    startTime={20000,40000,60000}) annotation (Placement(transformation(extent={{-212,-4},{-192,16}})));
  Flow.SourceMassFlowEC fuelSource(
    p0=1.1*sinkFuel.p0,
    T(displayUnit="K") = TStart,
    compositionNom={0.01228,0.00001,0.00001,0.00001,0.98772,0.00001},
    RC0=0.75,
    nCells=nCells,
    minCurrentForDosing=12.78,
    use_in_I=true) annotation (Placement(transformation(extent={{-140,42},{-120,62}})));
equation
  // Flanges connections
  connect(airSource.flange, stack.airInlet)
    annotation (Line(points={{-92,-54},{-32,-54},{-32,-9.8},{-22.4,-9.8}}, color={159,159,223}));
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

  connect(fuelSource.flange, stack.fuelInlet)
    annotation (Line(points={{-120,52},{-80,52},{-80,14},{-36,14},{-36,11.8},{-22.4,11.8}}, color={159,159,223}));
  connect(currentRamp.y, fuelSource.in_I) annotation (Line(points={{-191,6},{-144,6},{-144,57},{-136,57}}, color={0,0,127}));
  connect(currentRamp.y, powerElectronics.currentInput[1])
    annotation (Line(points={{-191,6},{-144,6},{-144,20},{-100.9,20},{-100.9,10.8}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end Stack0DTest;

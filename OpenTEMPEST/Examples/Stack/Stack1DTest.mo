within OpenTEMPEST.Examples.Stack;
model Stack1DTest "Simple example model using the model Stack1D, adiabatic"
  extends Modelica.Icons.Example;

  import SI = Modelica.SIunits;

  // Initial conditions
  parameter SI.Pressure pStart=100000;
  parameter SI.Temperature TStart(displayUnit="K")=973.15;

  // Stack
  parameter Integer Ncell=5;
  parameter Integer N(min=1)=5    "Number of axially discretized units in heat ports";

  // Components
  SOC.Stack.Stack1D stack(
    N=N,
    Ncell=Ncell,
    TStart=TStart,
    pStartOutFC=pStart,
    redeclare model Electrochem = OpenTEMPEST.SOC.Electrochem.Components.ASR_Steam)
                        annotation (Placement(transformation(extent={{-28,-26},{28,28}})));
  BOP.PowerElectronics powerElectronics(nChannels=1, currentFactor=1)
                                        annotation (Placement(transformation(extent={{-110,-12},{-84,12}})));
  Flow.SourceMassFlowEC fuelSource(
    p0=1.1*sinkFuel.p0,
    T(displayUnit="K") = TStart,
    compositionNom={0.01228,0.00001,0.00001,0.00001,0.98772,0.00001},
    RC0=0.75,
    nCells=Ncell,
    minCurrentForDosing=12.78,
    use_in_I=true) annotation (Placement(transformation(extent={{-112,32},{-92,52}})));
  ThermoPower.Gas.SourceMassFlow airSource(
    redeclare package Medium = OpenTEMPEST.Medium.Air_Medium,
    w0=Ncell*4.29048E-04,
    p0=1.1*sinkAir.p0,
    T(displayUnit="K") = TStart)
           annotation (Placement(transformation(extent={{-112,-64},{-92,-44}},
          rotation=0)));

  ThermoPower.Gas.SinkPressure sinkAir(
    redeclare package Medium = OpenTEMPEST.Medium.Air_Medium,
    p0=pStart,
    use_in_p0=false) annotation (Placement(transformation(extent={{88,-58},{108,-38}},rotation=0)));
  ThermoPower.Gas.SinkPressure sinkFuel(redeclare package Medium = OpenTEMPEST.Medium.Fuel_CH4, p0=pStart)
    annotation (Placement(transformation(extent={{88,32},{108,52}}, rotation=0)));
  ThermoPower.Electrical.Grid grid(Pgrid=55555555555) annotation (Placement(transformation(extent={{-82,-36},{-62,-16}})));
  Blocks.Sources.MultiRamp currentRamp(
    n=3,
    height=fill(-15, currentRamp.n),
    duration=fill(2000, currentRamp.n),
    startTime={20000,40000,60000}) annotation (Placement(transformation(extent={{-258,-8},{-238,12}})));
equation
  // Flanges connections

  // Electrical connections
  connect(powerElectronics.gridConnectionAC, grid.port)
    annotation (Line(
      points={{-85.3,-3.6},{-80,-3.6},{-80,-14},{-84,-14},{-84,-18},{-86,-18},{-86,-26},{-80.6,-26}},
      color={0,0,255},
      thickness=0.5));

  connect(fuelSource.flange, stack.fuelIn)
    annotation (Line(points={{-92,42},{-42,42},{-42,17.2},{-31.92,17.2}}, color={159,159,223}));
  connect(airSource.flange, stack.airIn) annotation (Line(points={{-92,-54},{-31.92,-54},{-31.92,-21.14}}, color={159,159,223}));
  connect(stack.airOut, sinkAir.flange)
    annotation (Line(points={{33.6,-15.74},{78,-15.74},{78,-48},{88,-48}}, color={159,159,223}));
  connect(stack.fuelOut, sinkFuel.flange) annotation (Line(points={{30.8,17.2},{78,17.2},{78,42},{88,42}}, color={159,159,223}));
  connect(powerElectronics.pinN[1], stack.pin_p)
    annotation (Line(points={{-93.1,-10.8},{-93.1,-40},{-44,-40},{-44,-11.96},{-34.16,-11.96}}, color={0,0,255}));
  connect(powerElectronics.pinP[1], stack.pin_n) annotation (Line(points={{-100.9,-10.8},{-100.9,-16},{-84,-16},{-84,-14},{-78,-14},
          {-78,-4.94},{-34.16,-4.94}}, color={0,0,255}));
  connect(currentRamp.y, powerElectronics.currentInput[1])
    annotation (Line(points={{-237,2},{-116,2},{-116,18},{-100.9,18},{-100.9,10.8}}, color={0,0,127}));
  connect(currentRamp.y, fuelSource.in_I)
    annotation (Line(points={{-237,2},{-116,2},{-116,26},{-118,26},{-118,47},{-108,47}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(StopTime=100000, __Dymola_Algorithm="Cvode"));
end Stack1DTest;

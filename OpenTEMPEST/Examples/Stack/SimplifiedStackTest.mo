within OpenTEMPEST.Examples.Stack;
model SimplifiedStackTest
  "Simple example model using the model blackBoxStack, adiabatic"
  extends Modelica.Icons.Example;

  import SI = Modelica.SIunits;

  // Initial conditions
  parameter SI.Pressure pStart=100000;
  parameter SI.Temperature TStart(displayUnit="K")=1073.15;

  // Stack
  parameter Integer nCells=30;
  parameter Integer n(min=3)=5    "Number of axially discretized units in heat ports";

  // Components
  SOC.Stack.SimplifiedStack1D
                          stack(
    N=n,                        Ncell=nCells,
    verticalBlockSize={10,3,1,1,4,10},
    NdetailedCell=1,
    TStart=TStart,                            redeclare model Electrochem =
        OpenTEMPEST.SOC.Electrochem.Components.BV_Steam)
    annotation (Placement(transformation(extent={{-28,-26},{28,28}})));
  BOP.PowerElectronics powerElectronics(nChannels=1, currentFactor=1)
                                        annotation (Placement(transformation(extent={{-110,-12},{-84,12}})));
  ThermoPower.Gas.SourceMassFlow
                        fuelSource(
    redeclare package Medium = OpenTEMPEST.Medium.Fuel_CH4,
    p0=pStart,
    T(displayUnit="K") = TStart,
    Xnom={0.012281,1e-5,1e-5,1e-5,0.987703,1e-5},
    w0=0.0002025395)
                   annotation (Placement(transformation(extent={{-112,32},{-92,52}})));
  ThermoPower.Gas.SourceMassFlow airSource(
    redeclare package Medium = OpenTEMPEST.Medium.Air_Medium,
    w0=1e-4,
    p0=pStart,
    T(displayUnit="K") = TStart)
           annotation (Placement(transformation(extent={{-112,-64},{-92,-44}},
          rotation=0)));

  Modelica.Blocks.Sources.Ramp     currentDensitySOFC(
    height=-50,
    duration=200,
    startTime=1e3)                                           annotation (Placement(transformation(extent={{-216,32},{-196,52}})));

  ThermoPower.Gas.SinkPressure sinkAir(
    redeclare package Medium = OpenTEMPEST.Medium.Air_Medium,
    p0=pStart,
    use_in_p0=false) annotation (Placement(transformation(extent={{88,-58},{108,-38}},rotation=0)));
  ThermoPower.Gas.SinkPressure sinkFuel(redeclare package Medium =
        OpenTEMPEST.Medium.Fuel_CH4,                                                            p0=pStart)
    annotation (Placement(transformation(extent={{88,32},{108,52}}, rotation=0)));
  ThermoPower.Electrical.Grid grid(Pgrid=55555555555) annotation (Placement(transformation(extent={{-82,-36},{-62,-16}})));
equation
  // Flanges connections

  // Electrical connections
  connect(powerElectronics.gridConnectionAC, grid.port)
    annotation (Line(
      points={{-85.3,-3.6},{-80,-3.6},{-80,-14},{-84,-14},{-84,-18},{-86,-18},{-86,-26},{-80.6,-26}},
      color={0,0,255},
      thickness=0.5));

  connect(fuelSource.flange, stack.fuelIn) annotation (Line(points={{-92,42},{-46,
          42},{-46,18.82},{-35.84,18.82}}, color={159,159,223}));
  connect(airSource.flange, stack.airIn) annotation (Line(points={{-92,-54},{-46,
          -54},{-46,-17.9},{-36.4,-17.9}}, color={159,159,223}));
  connect(stack.fuelOut, sinkFuel.flange) annotation (Line(points={{33.6,18.82},
          {78,18.82},{78,42},{88,42}}, color={159,159,223}));
  connect(stack.airOut, sinkAir.flange) annotation (Line(points={{33.6,-15.74},{
          78,-15.74},{78,-48},{88,-48}}, color={159,159,223}));
  connect(currentDensitySOFC.y, powerElectronics.currentInput[1]) annotation (
      Line(points={{-195,42},{-118,42},{-118,18},{-100.9,18},{-100.9,10.8}},
        color={0,0,127}));
  connect(powerElectronics.pinN[1], stack.pin_p) annotation (Line(points={{
          -93.1,-10.8},{-93.1,-40},{-48,-40},{-48,-12.5},{-35.28,-12.5}}, color
        ={0,0,255}));
  connect(powerElectronics.pinP[1], stack.pin_n) annotation (Line(points={{
          -100.9,-10.8},{-100.9,-16},{-84,-16},{-84,-14},{-78,-14},{-78,-1.16},
          {-35.28,-1.16}}, color={0,0,255}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(StopTime=10000, __Dymola_Algorithm="Cvode"));
end SimplifiedStackTest;

within OpenTEMPEST.BOP;
model PowerElectronics
  "Model of a power electronics unit including AC/DC conversion"

  parameter Integer nChannels=3 "Number of parallel electrical channels";
  parameter Real currentFactor=1 "Liniar factor, which is multiplied with the input current";

  parameter Real aLoss = 0.04821 "Parameter for AC/DC loss and aux power calculation of power; form: a * x^2 + b*x + c";
  parameter Real bLoss = 0.009846;
  parameter Real cLoss = 0.01151;
  parameter Real desginPower = 7200 "Desgin power for AC/DC loss and aux power calculation per channel";

  parameter Boolean useVoltageInput = false "Use a voltage set value instead a current set value";

  parameter Integer nCells = 720 "Number of total cells attached (only for UAvg calc)";

protected
      Modelica.Electrical.Analog.Basic.Ground ground annotation (
        Placement(transformation(extent={{-100,-100},{-80,-80}})));
public
      Modelica.Electrical.Analog.Sources.SignalCurrent signalCurrent[nChannels] if not
    useVoltageInput                                                  annotation (
        Placement(transformation(extent={{-10,-10},{10,10}},  rotation=0,    origin={0,0})));
      Modelica.Blocks.Interfaces.RealInput currentInput[nChannels](unit="A")
                                                                   if not
    useVoltageInput
    "The set current"                                              annotation (Placement(
            transformation(
            extent={{10,-10},{-10,10}},
            rotation=90,
            origin={4,62}), iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={-30,90})));
      Modelica.Electrical.Analog.Interfaces.NegativePin pinN[nChannels]
    "Negative pins"                                                     annotation (Placement(
            transformation(extent={{20,-100},{60,-60}}), iconTransformation(extent={{20,-100},
            {40,-80}})));
      Modelica.Electrical.Analog.Interfaces.PositivePin pinP[nChannels]
    "Positives pins"                                                    annotation (Placement(
            transformation(extent={{-60,-100},{-20,-60}}), iconTransformation(
              extent={{-40,-100},{-20,-80}})));
public
      ThermoPower.PowerPlants.Buses.Sensors sensors "Sensor bus"
        annotation (Placement(transformation(extent={{80,20},{100,40}}),
        iconTransformation(extent={{80,20},{100,40}})));
protected
  Modelica.Blocks.Sources.Constant constCurrentFactor[nChannels](each k=
        currentFactor)
    annotation (Placement(transformation(extent={{-6,-6},{6,6}},
        rotation=-90,
        origin={-12,56})));
public
  ThermoPower.Electrical.PowerConnection gridConnectionAC "AC grid connection"
    annotation (Placement(transformation(rotation=0, extent={{80,-40},{100,-20}}),
        iconTransformation(extent={{80,-40},{100,-20}})));
      Modelica.Blocks.Interfaces.RealInput voltageInput[nChannels]
 if useVoltageInput
    "The set voltage" annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={-34,62}), iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={30,90})));
  Modelica.Electrical.Analog.Sources.SignalVoltage signalVoltage[nChannels]
 if useVoltageInput
    annotation (Placement(transformation(extent={{10,-34},{-10,-14}})));
  Modelica.Electrical.Analog.Sensors.MultiSensor multiSensor[nChannels]
    annotation (Placement(transformation(extent={{-82,-44},{-62,-24}})));
protected
  Blocks.ComputationBlocks.AddX addPower(n=nChannels)
    annotation (Placement(transformation(extent={{40,20},{60,40}})));
protected
  ThermoPower.Electrical.Load load(Pnom(displayUnit="kW") = 50000,usePowerInput=
       true) annotation (Placement(transformation(extent={{80,16},{100,-4}})));
protected
  Modelica.Blocks.Math.Product productCurrent[nChannels]     if not
    useVoltageInput annotation (Placement(transformation(
        extent={{-6,-6},{6,6}},
        rotation=-90,
        origin={8.88178e-16,28})));
public
  Modelica.Blocks.Interfaces.RealOutput PDC(unit="W") "DC power" annotation (Placement(
        transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={46,66}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-90,50})));

public
  Modelica.Blocks.Interfaces.RealOutput PAC(unit="W") "AC power" annotation (Placement(
        transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={74,66}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-90,10})));
  Blocks.ComputationBlocks.Ploss pLoss[nChannels](
    each Pdesign=desginPower,
    each a=aLoss,
    each b=bLoss,
    each c=cLoss)
    annotation (Placement(transformation(extent={{40,8},{48,16}})));
public
  Modelica.Blocks.Interfaces.RealOutput UAvg(unit="V")
                                             "Average cell voltage"
                                                       annotation (Placement(
        transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={28,66}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-90,-40})));
protected
  Blocks.ComputationBlocks.AddX addPowerLoss(n=nChannels)
    annotation (Placement(transformation(extent={{52,8},{60,16}})));
  Modelica.Blocks.Math.Add add(k2=-1)
    annotation (Placement(transformation(extent={{64,18},{72,26}})));
  Modelica.Blocks.Math.Abs absP[nChannels]
    annotation (Placement(transformation(extent={{28,8},{36,16}})));
protected
  Modelica.Blocks.Math.Gain gain(k=-1)
    annotation (Placement(transformation(extent={{76,4},{80,8}})));
equation

  UAvg = abs(PDC) / nCells / (abs(sum(multiSensor.i)+Modelica.Constants.small)/nChannels);

  for i in 1:nChannels loop
    connect(pinP[i], ground.p)
      annotation (Line(points={{-40,-80},{-90,-80}}, color={0,0,255}));
  end for;

      // Sensor conncetions
      for i in 1:nChannels loop
      end for;

  connect(currentInput, productCurrent.u1)
    annotation (Line(points={{4,62},{4,35.2},{3.6,35.2}}, color={0,0,127}));
  connect(constCurrentFactor.y, productCurrent.u2) annotation (Line(points={{-12,
          49.4},{-4,49.4},{-4,35.2},{-3.6,35.2}},
                                                color={0,0,127}));
  connect(productCurrent.y, signalCurrent.i) annotation (Line(points={{-3.33067e-16,
          21.4},{-3.33067e-16,13.75},{0,13.75},{0,12}},color={0,0,127}));
  connect(load.port, gridConnectionAC) annotation (Line(
      points={{90,-2.6},{90,-30}},
      color={0,0,255},
      thickness=0.5));
  connect(signalCurrent.n, pinN) annotation (Line(points={{10,0},{40,0},{40,-80}},
                     color={0,0,255}));
//   connect(signalCurrent.p, voltageSensor.n) annotation (Line(points={{-10,0},{-10,
//           -12},{-10,-24},{-8,-24}}, color={0,0,255}));
//   connect(signalCurrent.n, voltageSensor.p)
//     annotation (Line(points={{10,0},{10,-24},{8,-24}}, color={0,0,255}));
//   connect(currentInput, productPowerCalc.u1) annotation (Line(points={{4,62},{10,
//           62},{10,12},{17,12}}, color={0,0,127}));
//   connect(voltageSensor.v, productPowerCalc.u2) annotation (Line(points={{1.77636e-15,
//           -15.2},{17,-15.2},{17,6}}, color={0,0,127}));

  connect(voltageInput, signalVoltage.v) annotation (Line(points={{-34,62},{-34,
          40},{-62,40},{-62,-12},{0,-12}}, color={0,0,127}));
  connect(signalVoltage.p, pinN)
    annotation (Line(points={{10,-24},{40,-24},{40,-80}}, color={0,0,255}));
  connect(multiSensor.pv, pinN)
    annotation (Line(points={{-72,-24},{40,-24},{40,-80}}, color={0,0,255}));
  connect(multiSensor.nv, pinP)
    annotation (Line(points={{-72,-44},{-72,-80},{-40,-80}}, color={0,0,255}));
  connect(multiSensor.pc, pinP)
    annotation (Line(points={{-82,-34},{-82,-66},{-40,-66},{-40,-80}},
                                                   color={0,0,255}));
  connect(multiSensor.nc, signalVoltage.n)
    annotation (Line(points={{-62,-34},{-62,-24},{-10,-24}}, color={0,0,255}));
  connect(multiSensor.nc, signalCurrent.p)
    annotation (Line(points={{-62,-34},{-62,0},{-10,0}}, color={0,0,255}));

  connect(multiSensor.power, addPower.u) annotation (Line(points={{-83,-40},{
          -40,-40},{-40,30},{38,30}},
                                color={0,0,127}));
  connect(addPower.y, PDC) annotation (Line(points={{61,30},{64,30},{64,66},{46,
          66}}, color={0,0,127}));
  connect(pLoss.y, addPowerLoss.u)
    annotation (Line(points={{48.4,12},{51.2,12}}, color={0,0,127}));
  connect(add.u2, addPowerLoss.y) annotation (Line(points={{63.2,19.6},{63.2,
          11.8},{60.4,11.8},{60.4,12}}, color={0,0,127}));
  connect(add.u1, addPower.y) annotation (Line(points={{63.2,24.4},{63.2,27.2},
          {61,27.2},{61,30}}, color={0,0,127}));
  connect(add.y, PAC)
    annotation (Line(points={{72.4,22},{74,22},{74,66}}, color={0,0,127}));
  connect(absP.y, pLoss.u)
    annotation (Line(points={{36.4,12},{39.2,12}}, color={0,0,127}));
  connect(absP.u, addPower.u) annotation (Line(points={{27.2,12},{22,12},{22,30},
          {38,30}}, color={0,0,127}));
  connect(add.y, gain.u) annotation (Line(points={{72.4,22},{74,22},{74,6},{
          75.6,6}}, color={0,0,127}));
  connect(gain.y, load.referencePower)
    annotation (Line(points={{80.2,6},{86.7,6}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
              extent={{-80,80},{80,-80}},
              fillColor={175,175,175},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None), Text(
              extent={{-80,80},{80,-80}},
              pattern=LinePattern.None,
              fillColor={175,175,175},
              fillPattern=FillPattern.Solid,
              lineColor={0,0,0},
              textString="Power 
Electronics")}),                                                     Diagram(
            coordinateSystem(preserveAspectRatio=false)));
end PowerElectronics;

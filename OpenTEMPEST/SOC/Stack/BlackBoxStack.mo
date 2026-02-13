within OpenTEMPEST.SOC.Stack;
model BlackBoxStack
  "0D stack based on lumped ASR expression assuming full internal steam reforming (CH4 outlet is 0) with no pressure loss."

  import SI = Modelica.SIunits;

  parameter Integer nCells = 30;
  parameter Integer nParallel = 1 "parallel stacks with equal voltage";
  parameter SI.Area ACell = 127.8/100^2;
  parameter SI.HeatCapacity C = 1 "capacity of each stack if nParallel>=1";
  parameter Units.AreaSpecificResistance alphaASR=728;
  parameter SI.LinearTemperatureCoefficient betaASR=-8.29e-3;

  parameter Integer n(min = 1)=10 "number of axially discretisized units in heat ports";
  parameter Integer intermediatePlateDistance = 10 "number of cells between two intermediate plates, set to >nCellPerStack to have none" annotation (Dialog(tab="Intermediate and end plates properties"));
  parameter Boolean calcPressureDrop = false "define if pressure drop calculation needed";
  parameter Boolean useDhtInletOutlet = true "define if dht at inlet and outlet are needed. For FMU, scaling nCells in DHT generates error 'the start values for the following variables could not be set: Ncell' when changing Ncell param";

  replaceable model ASRobj =
      Electrochem.ASR.ASR_Exponential                                   constrainedby
    Electrochem.ASR.ASR_Base                                                                                                              annotation(choicesAllmatching=true);

  ASRobj asr(Tpen= TASR, A=alphaASR, B=betaASR);

  parameter SI.MassFraction [Fuel.nXi] XStartGas=Fuel.X_default;
  parameter SI.MassFraction [Air.nXi] XStartAir=Air.X_default;

  replaceable package Fuel =
      Medium.Fuel_CH4         annotation(choicesAllMatching = true);
  replaceable package Air =
      Medium.Air_Medium         annotation(choicesAllMatching = true);

  SI.Voltage UOp;
  SI.Voltage UId;
  SI.Temperature TASR(start=825); // = sensAirOut.T;
  SI.MassFlowRate mfO2cr;
  Real XOutGas[6];
  Real XOutAir[Air.nXi];
  Real YOutGuess[6];
  SI.Current I;
  Real RC;
  Real z;
  Real ASR "Area specific resistance in ohm cm^2";

  SI.SpecificGibbsFreeEnergy go_T_H;
  SI.SpecificGibbsFreeEnergy delGr_H "Specific gibbs free energy change of reaction H2/H2O electrochem reaction";
  SI.AbsolutePressure Po = 1e5;
  parameter Real eps=Modelica.Constants.eps "compare vs. 0 in a conform way, needed for FMU import since other environments implement it differently with different eps. eps=1e-6 needed for simulink import as model exchange";

  ThermoPower.Gas.FlangeA airInlet(redeclare package Medium = Air) annotation (
      Placement(transformation(rotation=0, extent={{-80,-40},{-60,-20}}),
        iconTransformation(extent={{-100,-60},{-60,-20}})));
  ThermoPower.Gas.FlangeA fuelInlet(redeclare package Medium = Fuel)
    annotation (Placement(transformation(rotation=0, extent={{-80,40},{-60,60}}),
        iconTransformation(extent={{-100,20},{-60,60}})));
  ThermoPower.Gas.FlangeB fuelOutlet(redeclare package Medium = Fuel)
    annotation (Placement(transformation(rotation=0, extent={{80,40},{100,60}}),
        iconTransformation(extent={{60,20},{100,60}})));
  ThermoPower.Gas.FlangeB airOutlet(redeclare package Medium = Air) annotation (
     Placement(transformation(rotation=0, extent={{80,-40},{100,-20}}),
        iconTransformation(extent={{60,-60},{100,-20}})));
  Modelica.Electrical.Analog.Interfaces.PositivePin pinP annotation (Placement(
        transformation(extent={{-106,-34},{-86,-14}}), iconTransformation(
          extent={{-80,-20},{-70,-10}})));
  Modelica.Electrical.Analog.Interfaces.NegativePin pinN annotation (Placement(
        transformation(extent={{-106,14},{-86,34}}),iconTransformation(extent={{-80,10},
            {-70,20}})));
  OpenTEMPEST.Flow.SensGasProperty sensAirOut(
    mfOutput=false,
    pOutput=true,
    hOutput=false,
    XOutput=false,
    YOutput=true,
    HfOutput=true,
    redeclare package Medium = Medium.Air_Medium,
    cpOutput=false,
    lambdaOutput=false,
    etaOutput=true,
    TOutput=true,
    rhoOutput=true)
    annotation (Placement(transformation(extent={{52,-36},{72,-16}})));
  OpenTEMPEST.Flow.SensGasProperty sensAirIn(
    mfOutput=true,
    pOutput=true,
    hOutput=false,
    XOutput=true,
    YOutput=true,
    HfOutput=true,
    redeclare package Medium = Medium.Air_Medium,
    pstart=140000,
    cpOutput=false,
    lambdaOutput=false,
    etaOutput=true,
    TOutput=true,
    rhoOutput=true)
    annotation (Placement(transformation(extent={{-52,-36},{-32,-16}})));
  OpenTEMPEST.Flow.SensGasProperty sensFuelIn(
    mfOutput=true,
    pOutput=false,
    hOutput=false,
    XOutput=true,
    YOutput=true,
    HfOutput=true,
    redeclare package Medium = Fuel,
    pstart=140000,
    Xstart=XStartGas,
    cpOutput=false,
    lambdaOutput=false,
    etaOutput=true,
    TOutput=true,
    rhoOutput=true)
    annotation (Placement(transformation(extent={{-56,44},{-36,64}})));
  OpenTEMPEST.Flow.SensGasProperty sensFuelOut(
    mfOutput=false,
    pOutput=true,
    hOutput=false,
    XOutput=false,
    YOutput=true,
    HfOutput=true,
    redeclare package Medium = Fuel,
    Xstart=XStartGas,
    cpOutput=false,
    lambdaOutput=false,
    etaOutput=true,
    TOutput=false,
    rhoOutput=true)
    annotation (Placement(transformation(extent={{46,44},{66,64}})));
  ThermoPower.Gas.SinkPressure sinkPressure(
    redeclare package Medium = Fuel,
    Xnom=XStartGas,
    use_in_p0=true,
    use_in_T=true)
    annotation (Placement(transformation(extent={{-26,40},{-6,60}})));
  ThermoPower.Gas.SourceMassFlow sourceMassFlowFuel(
    redeclare package Medium = Fuel,
    Xnom=XStartGas,
    w0=1e-4,
    use_in_w0=true,
    use_in_T=true,
    use_in_X=true)
    annotation (Placement(transformation(extent={{18,40},{38,60}})));
  Modelica.Blocks.Sources.RealExpression rexTASR(y=TASR)
    annotation (Placement(transformation(extent={{-60,0},{-40,20}})));
  Modelica.Blocks.Sources.RealExpression rexMfO2cr(y=mfO2cr)
    annotation (Placement(transformation(extent={{-60,16},{-40,36}})));
  Modelica.Blocks.Math.Add add
    annotation (Placement(transformation(extent={{0,74},{12,86}})));
  Modelica.Blocks.Sources.RealExpression rexXOut[6](y=XOutGas)
    annotation (Placement(transformation(extent={{24,74},{44,94}})));
  ThermoPower.Gas.SinkPressure sinkPressure1(
    redeclare package Medium = Air,
    use_in_p0=true,
    use_in_T=true)
    annotation (Placement(transformation(extent={{-16,-40},{4,-20}})));
  ThermoPower.Gas.SourceMassFlow sourceMassFlowAir(
    redeclare package Medium = Medium.Air_Medium,
    w0=1e-4,
    use_in_w0=true,
    use_in_T=true,
    use_in_X=true)
    annotation (Placement(transformation(extent={{20,-40},{40,-20}})));
  Modelica.Blocks.Math.Add add1(k1=-1)
    annotation (Placement(transformation(extent={{-2,-16},{4,-10}})));
  Modelica.Blocks.Sources.RealExpression rexXOutAir[Air.nXi](y=XOutAir)
    annotation (Placement(transformation(extent={{12,-80},{32,-60}})));
  ThermoPower.Thermal.HT htSideLeft
    annotation (Placement(transformation(extent={{-70,76},{-50,96}}),
        iconTransformation(extent={{-46,60},{-34,72}})));
  ThermoPower.Thermal.HT htSideRight
    annotation (Placement(transformation(extent={{56,-98},{76,-78}}),
        iconTransformation(extent={{34,-72},{46,-60}})));
  ThermoPower.Thermal.DHTVolumes dhtBottom(N=n) annotation (Placement(
        transformation(extent={{-10,-100},{10,-80}}), iconTransformation(extent=
           {{-10,-80},{10,-60}})));
  ThermoPower.Thermal.DHTVolumes dhtTop(N=n) annotation (Placement(
        transformation(extent={{-10,80},{10,100}}), iconTransformation(extent={{
            -10,60},{10,80}})));
  ThermoPower.Thermal.DHTVolumes dhtInlet(N=if useDhtInletOutlet then nCells else 1)        annotation (
      Placement(transformation(extent={{-100,-10},{-80,10}}),
        iconTransformation(extent={{-80,-10},{-60,10}})));
  ThermoPower.Thermal.HT htPlatesInlet[if useDhtInletOutlet then 2 + integer((floor((nCells - 1)/intermediatePlateDistance))) else 1]
    annotation (Placement(transformation(extent={{-100,-14},{-80,6}}),
        iconTransformation(extent={{-92,-6},{-80,6}})));
  ThermoPower.Thermal.DHTVolumes dhtOutlet(N=if useDhtInletOutlet then nCells else 1)        annotation (
      Placement(transformation(extent={{80,-10},{100,10}}), iconTransformation(
          extent={{60,-10},{80,10}})));
  ThermoPower.Thermal.HT htPlatesOutlet[if useDhtInletOutlet then 2 + integer((floor((nCells - 1)/intermediatePlateDistance))) else 1]
    annotation (Placement(transformation(extent={{86,-12},{106,8}}),
        iconTransformation(extent={{80,-6},{92,6}})));
  Modelica.Blocks.Interfaces.RealOutput TMEACenterOut(unit="K",displayUnit="degC") annotation (Placement(
        transformation(extent={{40,-10},{60,10}}),iconTransformation(extent={{40,-10},
            {60,10}})));

  Modelica.Blocks.Sources.RealExpression pfuel_in(y=sensFuelOut.p +
        Blocks.Functions.pressureDropDarcy(
        sensFuelIn.mf/nCells/nParallel,
        (sensFuelIn.eta + sensFuelOut.eta)/2,
        (sensFuelIn.rho + sensFuelOut.rho)/2,
        0.112*2,
        0.0012*0.07143,
        (0.112*2/(0.0012*0.07143))*(3.76e-5/(1178917.95007*1e-3*1e5*0.185))))
    if calcPressureDrop
    annotation (Placement(transformation(extent={{-26,20},{-6,40}})));
  Modelica.Blocks.Sources.RealExpression pair_in(y=sensAirOut.p +
        Blocks.Functions.pressureDropDarcy(
        sensAirIn.mf/nCells/nParallel,
        (sensAirIn.eta + sensAirOut.eta)/2,
        (sensAirIn.rho + sensAirOut.rho)/2,
        0.112*2,
        0.0017*0.07143,
        (0.112*2/(0.0017*0.07143))*(4.2e-5/(63876e-3*1e5*0.367))))
    if calcPressureDrop
    annotation (Placement(transformation(extent={{-16,-60},{4,-40}})));
protected
  parameter Real zFC[Fuel.nX] = {2,8,0,2,0,0};
  parameter Real zEC[Fuel.nX] = {0,0,2,0,2,0};

equation

  // Electric connections
  nCells*UOp = pinP.v - pinN.v;
  0 = pinP.i + pinN.i;
  I = pinP.i;

  // Energy balance
  C*nParallel * der(TASR) = sensFuelIn.Hf + sensAirIn.Hf - sensFuelOut.Hf - sensAirOut.Hf - I * nCells*UOp*nParallel +
                    htSideLeft.Q_flow + htSideRight.Q_flow + sum(dhtBottom.Q)+ sum(dhtTop.Q) +
                   sum(dhtInlet.Q) + sum(htPlatesInlet.Q_flow) + sum(dhtOutlet.Q) + sum(htPlatesOutlet.Q_flow);

  // Electrochemical reaction gas side
  z = if I >= 0 then sum(zFC * sensFuelIn.y) else sum(zEC * sensFuelIn.y);
  RC = if z>eps then abs(I)*nCells*nParallel/z/Modelica.Constants.F / (sensFuelIn.mf/sensFuelIn.M + 1e-19) else 0;
  YOutGuess[1] = if I>=0 then sensFuelIn.y[1]*(1-RC) else sensFuelIn.y[1]+sensFuelIn.y[5]*RC;
  YOutGuess[2] = if I>=0 then sensFuelIn.y[2]*(1-RC) else 0;
  YOutGuess[3] = if I>=0 then sensFuelIn.y[3]+sensFuelIn.y[2]*RC+sensFuelIn.y[4]*RC else sensFuelIn.y[3]*(1-RC);
  YOutGuess[4] = if I>=0 then sensFuelIn.y[4]*(1-RC) else sensFuelIn.y[4]+sensFuelIn.y[3]*RC;
  YOutGuess[5] = if I>=0 then sensFuelIn.y[5]+sensFuelIn.y[1]*RC+2*sensFuelIn.y[2]*RC else sensFuelIn.y[5]*(1-RC);
  YOutGuess[6] = sensFuelIn.y[6];
  XOutGas =Blocks.Functions.simpleWGSEqn(TASR, Fuel.moleToMassFractions(
    YOutGuess, Fuel.MMX));

  // Electrochemical reaction air side
  mfO2cr = Air.MMX[1] * I * nCells*nParallel / 4 / Modelica.Constants.F;
  XOutAir[1] = (sensAirIn.x[1] * sensAirIn.mf - mfO2cr) / (sensAirIn.mf - mfO2cr + 1e-12);
  XOutAir[2] = 1 - XOutAir[1];

  // Gibbs energy calculation
  go_T_H = (0.05354*(TASR) -245.9767)*1e3; // linear relationship from NASA polynomials
  delGr_H = if z>eps then go_T_H + Modelica.Constants.R*TASR*(0.5*Modelica.Math.log( ((sensFuelIn.y[5]+sensFuelOut.y[5])/2)^2)
    - 0.5*Modelica.Math.log(((sensFuelIn.y[1]+sensFuelOut.y[1])/2)^2)  - 0.5*Modelica.Math.log((sensAirIn.y[1]+sensAirOut.y[1])/2)
    - 0.5*Modelica.Math.log((sensAirIn.p+sensAirOut.p)/2/Po)) else 0;

  // Voltage calculation
  UId = (-delGr_H/(2*Modelica.Constants.F));
  ASR = asr.ASR;
  //alphaASR*exp(betaASR*(TASR-273.15)); // Riedel, 2020, https://doi.org/10.1016/j.jpowsour.2020.228682
  UOp =  (UId - I/(ACell*100^2) * ASR);

  // External heat transfer
  htSideLeft.Q_flow = 0.15 * (htSideLeft.T-TASR);
  htSideRight.Q_flow = 0.15 * (htSideRight.T-TASR);
  dhtBottom.T[:] = fill(TASR,n);
  dhtTop.T[:] = fill(TASR,n);
  dhtInlet.T[:] = fill(TASR,if useDhtInletOutlet then nCells else 1);
  htPlatesInlet[:].T = fill(TASR,if useDhtInletOutlet then 2 + integer((floor((nCells - 1)/intermediatePlateDistance))) else 1);
  dhtOutlet.T[:] = fill(TASR,if useDhtInletOutlet then nCells else 1);
  htPlatesOutlet[:].T = fill(TASR,if useDhtInletOutlet then 2 + integer((floor((nCells - 1)/intermediatePlateDistance))) else 1);

  // Sensor output
  TMEACenterOut = TASR;

  connect(sensAirOut.outlet, airOutlet)
    annotation (Line(points={{68,-30},{90,-30}}, color={159,159,223}));
  connect(airInlet, sensAirIn.inlet)
    annotation (Line(points={{-70,-30},{-48,-30}}, color={159,159,223}));
  connect(fuelInlet, sensFuelIn.inlet)
    annotation (Line(points={{-70,50},{-52,50}}, color={159,159,223}));
  connect(sensFuelOut.outlet, fuelOutlet)
    annotation (Line(points={{62,50},{90,50}}, color={159,159,223}));
  connect(sensFuelIn.outlet, sinkPressure.flange)
    annotation (Line(points={{-40,50},{-26,50}}, color={159,159,223}));
  connect(sensFuelOut.inlet, sourceMassFlowFuel.flange)
    annotation (Line(points={{50,50},{38,50}}, color={159,159,223}));
  if calcPressureDrop then
    connect(pfuel_in.y, sinkPressure.in_p0);
  else
  connect(sensFuelOut.p, sinkPressure.in_p0) annotation (Line(points={{63,59},{63,
          58},{72,58},{72,68},{-22.45,68},{-22.45,55.95}},
                                           color={0,0,127}));
  end if;
  connect(sensFuelIn.mf, add.u1) annotation (Line(points={{-39,65},{-39,83.6},{-1.2,
          83.6}}, color={0,0,127}));
  connect(rexMfO2cr.y, add.u2) annotation (Line(points={{-39,26},{-32,26},{-32,76.4},
          {-1.2,76.4}}, color={0,0,127}));
  connect(sourceMassFlowFuel.in_w0, add.y)
    annotation (Line(points={{22,55},{22,80},{12.6,80}}, color={0,0,127}));
  connect(sourceMassFlowFuel.in_X, rexXOut.y) annotation (Line(points={{34,55},{
          34,76},{45,76},{45,84}}, color={0,0,127}));
  connect(sensAirIn.outlet, sinkPressure1.flange)
    annotation (Line(points={{-36,-30},{-16,-30}}, color={159,159,223}));
  connect(sensAirOut.inlet, sourceMassFlowAir.flange)
    annotation (Line(points={{56,-30},{40,-30}}, color={159,159,223}));
  connect(sensAirIn.mf, add1.u2) annotation (Line(points={{-35,-15},{-35,-14.8},
          {-2.6,-14.8}}, color={0,0,127}));
  connect(rexMfO2cr.y, add1.u1) annotation (Line(points={{-39,26},{-22,26},{-22,
          -11.2},{-2.6,-11.2}}, color={0,0,127}));
  connect(add1.y, sourceMassFlowAir.in_w0)
    annotation (Line(points={{4.3,-13},{24,-13},{24,-25}}, color={0,0,127}));
  connect(rexXOutAir.y, sourceMassFlowAir.in_X) annotation (Line(points={{33,-70},
          {50,-70},{50,-25},{36,-25}}, color={0,0,127}));
  if calcPressureDrop then
    connect(pair_in.y, sinkPressure1.in_p0);
  else
  connect(sinkPressure1.in_p0, sensAirOut.p) annotation (Line(points={{-12.45,-24.05},
          {-12.45,-6},{69,-6},{69,-21}}, color={0,0,127}));
  end if;
  connect(rexTASR.y, sourceMassFlowFuel.in_T) annotation (
    Line(points={{-39,10},{12,10},{12,64},{28,64},{28,55}},            color = {0, 0, 127}));
  connect(rexTASR.y, sourceMassFlowAir.in_T) annotation (
    Line(points={{-39,10},{30,10},{30,-25}},        color = {0, 0, 127}));
  connect(sensAirIn.T, sinkPressure1.in_T) annotation (
    Line(points={{-33,-20},{-20,-20},{-20,-21},{-6,-21}},
                                           color = {0, 0, 127}));
  connect(sensFuelIn.T, sinkPressure.in_T) annotation (
    Line(points={{-37,60},{-26,60},{-26,59},{-16,59}},
                                          color = {0, 0, 127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(pattern = LinePattern.None,
            fillPattern =                                                                                                     FillPattern.Solid, extent = {{-60, 60}, {60, -60}}), Text(lineColor = {255, 255, 255}, pattern = LinePattern.None,
            fillPattern =                                                                                                                                                                                                        FillPattern.Solid, extent = {{40, 40}, {-40, -40}}, textString = "0D-Stack")}),                              Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    Documentation(revisions="<html>
<ul>
<li>03-02-2023 by Santiago Salas Ventura<br>for FMU import into simulink, compare vs. eps instead of 0, since simulink has different machine precision and &quot;I&gt;0 or I&lt;0&quot; is not triggered if 0 is not replaced by ~1e-10 (though dymola is more forgiving, and that works).</li>
<li>31-01-2022 by Santiago Salas Ventura<br>add ASR class for different ASR expression options, avoid cooling when I=0 with small amount of CH4 by adding conditional.</li>
<li>23-09-2021 by Marius Tomberg<br>New model.</li>
</ul>
</html>"));
end BlackBoxStack;

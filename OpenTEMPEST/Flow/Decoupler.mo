within OpenTEMPEST.Flow;
model Decoupler "Decouples two parts of a model"
  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
    annotation (choicesAllMatching=true);
  parameter Real massFlowFactor=1
    "Scaling of the mass flow (e. g. 2 --> output = input * 2)";
  parameter Real fixedOutputTemperature=-1
    "positive: T in K, negative: inlet T";
  parameter Real fixedInletPressure=-1 "positive: p in Pa, negative: outlet p";
  parameter Modelica.SIunits.MassFraction xStart[Medium.nX]=Medium.reference_X
    "Start mass fraction: Linear fading to simulated composition in startTime";
  parameter Modelica.SIunits.Pressure pStartSourceMassFlow=102600;
  parameter Modelica.SIunits.Pressure pStartSinkPressure=101325;
  parameter Modelica.SIunits.Temperature Tstart=1073.15;
  parameter Modelica.SIunits.MassFlowRate mfOutletStart=0
    "Start mass flow: Linear fading to simulated mass flow in startTime";
  parameter Modelica.SIunits.Time startTime=1 "Time for start value fading, set to 0 if no fading necessary";

  ThermoPower.Gas.FlangeA flangeA(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-100,-10},{-80,10}}),
        iconTransformation(extent={{-80,-20},{-40,20}})));
  ThermoPower.Gas.FlangeB flangeB(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{80,-10},{100,10}}),
        iconTransformation(extent={{40,-20},{80,20}})));
  ThermoPower.Gas.SensW sensW(redeclare package Medium = Medium,
      allowFlowReversal=false)
    annotation (Placement(transformation(extent={{-70,-6},{-50,14}})));
  ThermoPower.Gas.SourceMassFlow sourceMassFlow(
    redeclare package Medium = Medium,
    p0=pStartSourceMassFlow,
    T=Tstart,
    Xnom=xStart,
    use_in_w0=true,
    use_in_T=true,
    use_in_X=true)
    annotation (Placement(transformation(extent={{20,10},{40,-10}})));
  ThermoPower.Gas.SinkPressure sinkPressure(redeclare package Medium = Medium,
  p0=pStartSinkPressure,
      use_in_p0=true,
    use_in_T=true)
    annotation (Placement(transformation(extent={{-18,10},{2,-10}})));
  ThermoPower.Gas.SensP sensP(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{60,-6},{80,14}})));
  ThermoPower.Gas.SensT sensT(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-50,-6},{-30,14}})));
  ThermoPower.Gas.SensT sensT1(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{42,-6},{62,14}})));
equation

  //for i in 1:Medium.nX loop
  //sourceMassFlow.in_X[i] = firstOrderX[i].u; //sensW.outlet.Xi_outflow;
  //firstOrderX[i].y = sensW.outlet.Xi_outflow[i];
  //end for;
  //   if time < 1 then
  //     sourceMassFlow.in_X = xStart;
  //     sourceMassFlow.in_w0 = mfOutletStart;
  //   else
  //     sourceMassFlow.in_X = sensW.outlet.Xi_outflow;
  //     sourceMassFlow.in_w0 = sensW.w*massFlowFactor;
  //   end if;
  sourceMassFlow.in_X = noEvent(if time < startTime then (startTime - time)/startTime*
    xStart + time/startTime*sensW.outlet.Xi_outflow else sensW.outlet.Xi_outflow);
  sourceMassFlow.in_w0 = noEvent(if time < startTime then (startTime - time)/startTime*mfOutletStart + time/startTime*sensW.w else
    sensW.w)*massFlowFactor;
  //sourceMassFlow.in_X = smooth(0, if time<1 then xStart else sensW.outlet.Xi_outflow);
  //sourceMassFlow.in_w0 = smooth(0, if time<1 then mfOutletStart else sensW.w*massFlowFactor);

  if fixedOutputTemperature < 0 then
    sourceMassFlow.in_T = sensT.T;
  else
    sourceMassFlow.in_T = fixedOutputTemperature;
  end if;
  if fixedInletPressure < 0 then
    sinkPressure.in_p0 = sensP.p;
  else
    sinkPressure.in_p0 = fixedInletPressure;
  end if;

  connect(sensW.inlet, flangeA)
    annotation (Line(points={{-66,0},{-90,0}}, color={159,159,223}));
  connect(sinkPressure.flange, sensT.outlet)
    annotation (Line(points={{-18,0},{-34,0}}, color={159,159,223}));
  connect(sensT.inlet, sensW.outlet)
    annotation (Line(points={{-46,0},{-54,0}}, color={159,159,223}));
  connect(flangeB, sensP.flange)
    annotation (Line(points={{90,0},{70,0}}, color={159,159,223}));
  connect(flangeB, sensT1.outlet)
    annotation (Line(points={{90,0},{58,0}}, color={159,159,223}));
  connect(sensT1.inlet, sourceMassFlow.flange)
    annotation (Line(points={{46,0},{40,0}}, color={159,159,223}));
  connect(sensT1.T, sinkPressure.in_T) annotation (Line(points={{59,10},{58,10},
          {58,-18},{-8,-18},{-8,-9}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
          extent={{-40,32},{40,-30}},
          fillColor={170,213,255},
          fillPattern=FillPattern.Solid,
          pattern=LinePattern.None)}), Diagram(coordinateSystem(
          preserveAspectRatio=false)));
end Decoupler;

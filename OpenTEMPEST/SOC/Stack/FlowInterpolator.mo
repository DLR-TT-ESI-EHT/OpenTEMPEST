within OpenTEMPEST.SOC.Stack;
model FlowInterpolator

  import SI = Modelica.SIunits;

  parameter Integer nCell(min = 3) "total number of single cells";
  parameter Integer nSimplified(min = 1) "number of simplified cells among nCell";
  parameter Integer nNonUnitOrSimpCells = 2 "Number of Top and Bottom Cell models";
  parameter Boolean isRedu[nCell-nNonUnitOrSimpCells] = {true,false,true};

  parameter SI.Temperature TStart;
  parameter SI.AbsolutePressure pStart;
  parameter SI.AbsolutePressure pStartAir=pStart "Starting Pressure on the air side" annotation (
    Dialog(tab = "Initialization"));
  parameter SI.AbsolutePressure pStartFuel=pStart "Starting Pressure on the fuel side" annotation (
    Dialog(tab = "Initialization"));
  parameter SI.MassFraction XFuelStart[Fuel.nX]=Fuel.reference_X;

  parameter Integer posDetailed=integer(floor(nCell/2))
                                  "Position of the detailed cell";

  replaceable package Fuel =
      OpenTEMPEST.Medium.Fuel_CH4;
  replaceable package Air =
      OpenTEMPEST.Medium.Air_Medium;

  Real TempDistriFactor(start=0);
  SI.Temperature TDetailedAvg(start=TStart);

  Integer detailedIdx[nCell-nSimplified-nNonUnitOrSimpCells] = Modelica.Math.BooleanVectors.index(not isRedu);

  ThermoPower.Gas.FlangeB fuelExit(redeclare package Medium = Fuel) annotation (
     Placement(transformation(extent={{90,40},{110,60}}), iconTransformation(
          extent={{60,20},{100,60}})));
  ThermoPower.Gas.FlangeB airExit(redeclare package Medium = Air) annotation (
      Placement(transformation(extent={{90,-60},{110,-40}}), iconTransformation(
          extent={{60,-60},{100,-20}})));
  ThermoPower.Gas.FlangeA fuelInlet(redeclare package Medium = Fuel)
    annotation (Placement(transformation(extent={{-100,40},{-80,60}}), iconTransformation(extent={{-120,20},{-80,60}})));
  ThermoPower.Gas.FlangeA airInlet(redeclare package Medium = Air)
    annotation (Placement(transformation(extent={{-110,-60},{-90,-40}}), iconTransformation(extent={{-120,-60},{-80,-20}})));
  ThermoPower.Gas.SinkPressure sinkPressure(redeclare package Medium = Fuel,
    p0=pStartFuel,
    T=TStart,
    Xnom=XFuelStart,
      use_in_p0=true)
    annotation (Placement(transformation(extent={{-40,40},{-20,60}})));
  ThermoPower.Gas.SinkPressure sinkPressure1(redeclare package Medium = Air,
    p0=pStartAir,
    T=TStart,
      use_in_p0=true)
    annotation (Placement(transformation(extent={{-40,-60},{-20,-40}})));
  ThermoPower.Gas.SourceMassFlow sourceMassFlowFuel(
    redeclare package Medium = Fuel,
    p0=pStartFuel,
    Xnom=XFuelStart,
    w0=0.001,
    use_in_w0=true,
    use_in_T=true,
    use_in_X=true)
    annotation (Placement(transformation(extent={{20,40},{40,60}})));
  ThermoPower.Gas.SourceMassFlow sourceMassFlowAir(
    redeclare package Medium = Air,
    p0=pStartAir,
    w0=0.001,
    use_in_w0=true,
    use_in_T=true,
    use_in_X=true)
    annotation (Placement(transformation(extent={{28,-60},{48,-40}})));
protected
  ThermoPower.Gas.SensP sensP(redeclare package Medium = Fuel)
    annotation (Placement(transformation(extent={{50,44},{70,64}})));
  ThermoPower.Gas.SensP sensP1(redeclare package Medium = Air)
    annotation (Placement(transformation(extent={{50,-56},{70,-36}})));

  ThermoPower.Gas.SensW sensW1(redeclare package Medium = Air)
    annotation (Placement(transformation(extent={{-64,-56},{-44,-36}})));
  Modelica.Blocks.Math.Product product
    annotation (Placement(transformation(extent={{-10,52},{0,62}})));
  Modelica.Blocks.Math.Product product1
    annotation (Placement(transformation(extent={{-10,-42},{0,-32}})));
  Modelica.Blocks.Sources.Constant const(k=nCell/(nCell - nSimplified))
    annotation (Placement(transformation(extent={{-70,-6},{-58,6}})));
public
  ThermoPower.Thermal.HT hT[nCell]
    annotation (Placement(transformation(extent={{-100,-10},{-80,10}}), iconTransformation(extent={{-100,-10},{-80,10}})));
protected
  OpenTEMPEST.Flow.SensGasProperty sensTFuel(
    mfOutput=true,
    pOutput=false,
    hOutput=false,
    XOutput=true,
    YOutput=false,
    HfOutput=false,
    redeclare package Medium = Fuel,
    pstart=pStartFuel,
    Tstart=TStart,
    Xstart=XFuelStart,
    cpOutput=false,
    lambdaOutput=false,
    etaOutput=false,
    TOutput=true,
    rhoOutput=false)
    annotation (Placement(transformation(extent={{-80,44},{-60,64}})));
  ThermoPower.Gas.SensT sensTAir(redeclare package Medium = Air,
      allowFlowReversal=false)
    annotation (Placement(transformation(extent={{-84,-56},{-64,-36}})));
public
  Modelica.Blocks.Sources.RealExpression RealExp2(y=(TOutAirSimp*nSimplified +
        sensTAir.T*(nCell - nSimplified))/nCell)
    annotation (Placement(transformation(extent={{-144,-36},{-52,-16}})));
  Modelica.Blocks.Interfaces.RealInput TOutAirSimp
    annotation (Placement(transformation(extent={{-162,-16},{-122,24}})));
public
  Modelica.Blocks.Sources.RealExpression RealExp3(y=(TOutFuelSimp*nSimplified
         + sensTFuel.T*(nCell - nSimplified))/nCell)
    annotation (Placement(transformation(extent={{-80,72},{-60,92}})));
  Modelica.Blocks.Interfaces.RealInput TOutFuelSimp
    annotation (Placement(transformation(extent={{-156,36},{-116,76}})));
  Modelica.Blocks.Math.UnitConversions.To_degC to_degC
    annotation (Placement(transformation(extent={{-120,78},{-100,98}})));
  Modelica.Blocks.Math.UnitConversions.To_degC to_degC1
    annotation (Placement(transformation(extent={{-168,-56},{-148,-36}})));
equation

  sourceMassFlowAir.in_X = sensW1.outlet.Xi_outflow;

  //assert(sourceMassFlowAir.in_T <= 1, "TAir in interpolator too small", AssertionLevel.error);
  //assert(sourceMassFlowFuel.in_T <= 1, "TFuel in interpolator too small", AssertionLevel.error);

  for j in 1:nCell loop

    hT[j].Q_flow = 0;

  end for;

  if nNonUnitOrSimpCells==2 then
    TDetailedAvg = (sum(hT[detailedIdx].T) + hT[1].T + hT[nCell].T)/(nCell-nSimplified);
  else
    TDetailedAvg = (sum(hT[detailedIdx].T))/(nCell-nSimplified);
  end if;

  TempDistriFactor = (TDetailedAvg-sum(hT[:].T)/nCell)/(sum(hT[:].T)/nCell);

  connect(airExit, airExit)
    annotation (Line(points={{100,-50},{100,-50}}, color={159,159,223}));
  connect(sourceMassFlowAir.flange, airExit) annotation (Line(points={{48,-50},
          {100,-50}},                  color={159,159,223}));
  connect(sourceMassFlowFuel.flange, fuelExit) annotation (Line(points={{40,50},
          {100,50}},                 color={159,159,223}));
  connect(sourceMassFlowFuel.flange, sensP.flange)
    annotation (Line(points={{40,50},{60,50}}, color={159,159,223}));
  connect(sourceMassFlowAir.flange, sensP1.flange)
    annotation (Line(points={{48,-50},{60,-50}}, color={159,159,223}));
  connect(sensP1.p, sinkPressure1.in_p0) annotation (Line(points={{67,-40},{70,
          -40},{70,-10},{-40,-10},{-40,-44.05},{-36.45,-44.05}},
                                              color={0,0,127}));
  connect(sensP.p, sinkPressure.in_p0) annotation (Line(points={{67,60},{80,60},
          {80,86},{-40,86},{-40,55.95},{-36.45,55.95}},
                                       color={0,0,127}));
  connect(sinkPressure1.flange, sensW1.outlet)
    annotation (Line(points={{-40,-50},{-48,-50}}, color={159,159,223}));
  connect(product.y, sourceMassFlowFuel.in_w0) annotation (Line(points={{0.5,57},
          {23.25,57},{23.25,55},{24,55}}, color={0,0,127}));
  connect(sensW1.w, product1.u2)
    annotation (Line(points={{-47,-40},{-11,-40}}, color={0,0,127}));
  connect(product1.y, sourceMassFlowAir.in_w0) annotation (Line(points={{0.5,-37},
          {32.25,-37},{32.25,-45},{32,-45}}, color={0,0,127}));
  connect(const.y, product.u2) annotation (Line(points={{-57.4,0},{-20,0},{-20,54},
          {-11,54}}, color={0,0,127}));
  connect(product1.u1, product.u2) annotation (Line(points={{-11,-34},{-20,-34},
          {-20,54},{-11,54}}, color={0,0,127}));
  connect(sensTFuel.inlet, fuelInlet)
    annotation (Line(points={{-76,50},{-90,50}},  color={159,159,223}));
  connect(sensW1.inlet, sensTAir.outlet)
    annotation (Line(points={{-60,-50},{-68,-50}}, color={159,159,223}));
  connect(sensTAir.inlet, airInlet)
    annotation (Line(points={{-80,-50},{-100,-50}}, color={159,159,223}));
  connect(RealExp2.y, sourceMassFlowAir.in_T)
    annotation (Line(points={{-47.4,-26},{38,-26},{38,-45}}, color={0,0,127}));
  connect(RealExp3.y, sourceMassFlowFuel.in_T)
    annotation (Line(points={{-59,82},{30,82},{30,55}}, color={0,0,127}));
  connect(sensTFuel.outlet, sinkPressure.flange) annotation (Line(points={{-64,50},
          {-54,50},{-54,50},{-40,50}}, color={159,159,223}));
  connect(sensTFuel.mf, product.u1)
    annotation (Line(points={{-63,65},{-11,65},{-11,60}}, color={0,0,127}));
  connect(sensTFuel.x, sourceMassFlowFuel.in_X) annotation (Line(points={{-63,
          57},{-48,57},{-48,74},{36,74},{36,55}}, color={0,0,127}));
  connect(TOutFuelSimp, to_degC.u) annotation (Line(points={{-136,56},{-130,56},{
          -130,88},{-122,88}}, color={0,0,127}));
  connect(TOutAirSimp, to_degC1.u) annotation (Line(points={{-142,4},{-156,4},{-156,
          -46},{-170,-46}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{60,80},{-80,-80}},
          lineColor={28,108,200},
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid)}),                      Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    Documentation(revisions="<html>
<li><i>04 May 2022</i>
by Santiago Salas Ventura:<br>
Define parameters pStartFuel and pStartAir which default to pStart; merge request !240.</li>
</html>"));
end FlowInterpolator;

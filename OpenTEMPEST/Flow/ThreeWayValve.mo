within OpenTEMPEST.Flow;
model ThreeWayValve
  "time controlled three-way valve for switching between two flow paths 'purge' and 'use'"

  import SI = Modelica.SIunits;
  parameter Boolean useTimeInput = true "Use the input connector for opening direction";
  parameter SI.Time fixTime=1e4 "Fixed start time when the input connector not used" annotation (Dialog(group="Valve Opening", enable=not useTimeInput));
  parameter SI.MassFraction X7 = 0 "mass fraction value of methanol added in reactor path";
  SI.MassFraction XUse[sourceMFUse.Medium.nX] "mass fractions for gas with added species methanol";
  Real sumX "auxiliary internal variable";

protected
  Modelica.Blocks.Interfaces.RealInput Tstart "Use connector input for the start time";

public
  Modelica.Blocks.Interfaces.RealInput startTime if useTimeInput annotation (Placement(
        transformation(
        origin={0,72},
        extent={{-10,-10},{10,10}},
        rotation=270), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={0,60})));
  ThermoPower.Gas.FlangeA infl(redeclare package Medium = Medium.Fuel_CH4)
    annotation (Placement(transformation(extent={{-128,-14},{-100,14}}),
        iconTransformation(extent={{-120,-20},{-80,20}})));
  ThermoPower.Gas.FlangeB outflPurge(redeclare package Medium = Medium.Fuel_CH4)
    annotation (Placement(transformation(extent={{88,-14},{116,14}}),
        iconTransformation(extent={{80,-20},{120,20}})));
  ThermoPower.Gas.SinkPressure sinkPressure(redeclare package Medium =
        Medium.Fuel_CH4, use_in_p0=true)
    annotation (Placement(transformation(extent={{-44,-10},{-24,10}})));
  ThermoPower.Gas.FlangeB outflUse(redeclare package Medium =
        Medium.Fuel_MethanolReactor) annotation (Placement(transformation(
          extent={{-14,-90},{14,-62}}), iconTransformation(extent={{-20,-120},{
            20,-80}})));
  OpenTEMPEST.Flow.SensGasProperty sensXTM(
    pOutput=false,
    hOutput=false,
    YOutput=false,
    HfOutput=false,
    redeclare package Medium = Medium.Fuel_CH4,
    cpOutput=false,
    lambdaOutput=false,
    etaOutput=false,
    rhoOutput=false)
    annotation (Placement(transformation(extent={{-96,-4},{-82,10}})));
  ThermoPower.Gas.SensP sensPPurge(redeclare package Medium = Medium.Fuel_CH4)
    annotation (Placement(transformation(extent={{84,-4},{68,12}})));
  ThermoPower.Gas.SensP sensPUse(redeclare package Medium =
        Medium.Fuel_MethanolReactor) annotation (Placement(transformation(
        extent={{8,8},{-8,-8}},
        rotation=-90,
        origin={-14,-50})));
  ThermoPower.Gas.SourceMassFlow sourceMFPurge(
    redeclare package Medium = Medium.Fuel_CH4,
    use_in_w0=true,
    use_in_T=true,
    use_in_X=true)
    annotation (Placement(transformation(extent={{40,-10},{60,10}})));
  ThermoPower.Gas.SourceMassFlow sourceMFUse(
    redeclare package Medium = Medium.Fuel_MethanolReactor,
    Xnom={0.12616,1e-05,0.2969,0.5769,1e-05,1e-05,1e-05},
    use_in_w0=true,
    use_in_T=true,
    use_in_X=true) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={0,-36})));
  inner ThermoPower.System system
    annotation (Placement(transformation(extent={{-100,80},{-80,100}})));

equation

    sumX = sum(sensXTM.x[:])+X7;
    for i in 1:6 loop
      XUse[i]= sensXTM.x[i]/sumX;
    end for;
    XUse[7] = X7;

  if time < Tstart then
    sourceMFPurge.in_w0 = sensXTM.mf;
//     sourceMFUse.in_w0 = sensXTM.mf;
    sourceMFUse.in_T = 800; //
    sourceMFUse.in_X = sourceMFUse.Xnom;
//     sinkPressure.in_p0 = sensPPurge.p;
  else
    sourceMFPurge.in_w0 = 0;
//     sourceMFUse.in_w0 = sensXTM.mf; //+XUse[7]*sensXTM.mf;
    sourceMFUse.in_T = sensXTM.T;
     sourceMFUse.in_X = XUse;
//     sinkPressure.in_p0 = sensPUse.p;
  end if;

 // Valve opening
  connect(startTime,Tstart);
  if not useTimeInput then
    Tstart = fixTime;
  end if;

  connect(infl, sensXTM.inlet) annotation (Line(points={{-114,0},{-104,0},{-104,
          0.2},{-93.2,0.2}}, color={159,159,223}));

  connect(sourceMFPurge.flange, outflPurge) annotation (Line(points={{60,0},{82,0},{82,1.77636e-15},{
          102,1.77636e-15}}, color={159,159,223}));
  connect(sensXTM.x, sourceMFPurge.in_X) annotation (Line(points={{-84.1,5.1},{-80,
          5.1},{-80,20},{56,20},{56,5}}, color={0,0,127}));

  connect(sourceMFUse.flange, outflUse) annotation (Line(points={{0,-46},{0,-76}}, color={159,159,223}));

  connect(sourceMFPurge.flange, sensPPurge.flange) annotation (Line(points={{60,0},{76,0},{76,0.8}}, color={159,159,223}));
  connect(sourceMFUse.flange, sensPUse.flange)
    annotation (Line(points={{0,-46},{0,-51},{-10.8,-51},{-10.8,-50}},
                                                                   color={159,159,223}));
  connect(sensPUse.p, sinkPressure.in_p0)
    annotation (Line(points={{-18.8,-44.4},{-18.8,10},{-40.45,10},{-40.45,5.95}}, color={0,0,127}));
  connect(sensXTM.T, sourceMFPurge.in_T) annotation (Line(points={{-82.7,7.2},{-74,
          7.2},{-74,16},{50,16},{50,5}}, color={0,0,127}));
  connect(sensXTM.mf, sourceMFUse.in_w0) annotation (Line(points={{-84.1,10.7},{
          12,10.7},{12,-30},{5,-30}}, color={0,0,127}));
  connect(sensXTM.outlet, sinkPressure.flange) annotation (Line(points={{-84.8,0.2},
          {-65.4,0.2},{-65.4,0},{-44,0}}, color={159,159,223}));
  annotation (                    Icon(graphics={
        Line(
          points={{0,30},{0,0}},
          color={0,0,0},
          thickness=0.5),
        Polygon(
          points={{-80,40},{-80,-40},{0,0},{-80,40}},
          lineColor={128,128,128},
          lineThickness=0.5,
          fillColor={159,159,223},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{80,40},{0,0},{80,-40},{80,40}},
          lineColor={128,128,128},
          lineThickness=0.5,
          fillColor={159,159,223},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-20,50},{20,30}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-40,40},{-40,-40},{40,0},{-40,40}},
          lineColor={128,128,128},
          lineThickness=0.5,
          fillColor={159,159,223},
          fillPattern=FillPattern.Solid,
          origin={0,-40},
          rotation=90)}));
end ThreeWayValve;

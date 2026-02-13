within OpenTEMPEST.Flow;
model GasHeatedSplitter

  import SI = Modelica.SIunits;

  replaceable package medium =
      Modelica.Media.Interfaces.PartialMedium                          annotation(choicesAllMatching = true);

  parameter SI.Temperature TStart=1000 "cell Starting inlet temperature for simulation"
    annotation (Dialog(tab="Initialization"));
  parameter SI.AbsolutePressure pStart=101325 "cell Starting Pressure"
    annotation (Dialog(tab="Initialization"));
  parameter SI.MassFraction XStart[medium.nX]=medium.reference_X annotation (Dialog(tab="Initialization"));

  parameter Real outletAShare = 0.999;
  parameter Boolean[3] pipe={true,true,false}
    "= true if flow through a pipe, false if between two parallel plates";
  parameter Integer[3] nParallel={2,2,30} "number of parallel flow ducts";
  parameter SI.Length[3] l={0.09,0.09,0.09} "lenght of flow channel";
  parameter SI.Length[3] h={0.01,0.01,0.01} "height of flow channel between two parallel plates";
  parameter SI.Length[3] w={0.04,0.04,0.04} "width of flow channel between two parallel plates";
  parameter SI.Length[3] d={0.02,0.02,0.02} "diameter of flow channel in a pipe";
  parameter SI.MassFlowRate maximalFlow = -1 "maximal mass flow rate";

  ThermoPower.Gas.FlangeA inlet(redeclare package Medium = medium)
                                                                 annotation (
      Placement(transformation(extent={{-114,-14},{-86,14}}),iconTransformation(
          extent={{-80,-20},{-40,20}})));
  ThermoPower.Gas.FlangeB outletA(redeclare package Medium = medium)
    annotation (Placement(transformation(extent={{52,32},{80,60}}),
        iconTransformation(extent={{40,20},{80,60}})));
  ThermoPower.Gas.FlangeB outletB(redeclare package Medium = medium)
    annotation (Placement(transformation(extent={{52,-48},{80,-20}}),
        iconTransformation(extent={{40,-60},{80,-20}})));
  ThermoPower.Thermal.HT ht[3] annotation (Placement(transformation(extent={{-10,30},
            {10,50}}),     iconTransformation(extent={{-18,-20},{22,20}})));

  ThermoPower.Gas.FlowSplit flowSplit(redeclare package Medium = medium)
    annotation (Placement(transformation(extent={{-14,-10},{6,10}})));
  OpenTEMPEST.BOP.GasConditioning.SimpleGasHeater simpleGasHeaterA(
    redeclare package medium = medium,
    l=l[2],
    h=h[2],
    w=w[2],
    d=d[2],
    TStart=TStart,
    pStart=pStart,
    each xStart=XStart,
    maximalFlow=maximalFlow,
    pipe=pipe[2],
    nParallel=nParallel[2],
    useAlphaIn=false,
    nHT=1) annotation (Placement(transformation(extent={{18,14},{38,34}})));
  OpenTEMPEST.BOP.GasConditioning.SimpleGasHeater simpleGasHeaterB(
    redeclare package medium = medium,
    l=l[3],
    h=h[3],
    w=w[3],
    d=d[3],
    TStart=TStart,
    pStart=pStart,
    each xStart=XStart,
    maximalFlow=maximalFlow,
    pipe=pipe[3],
    nParallel=nParallel[3],
    useAlphaIn=false,
    nHT=1) annotation (Placement(transformation(extent={{26,-42},{46,-22}})));
  OpenTEMPEST.BOP.GasConditioning.SimpleGasHeater simpleGasHeaterA2(
    redeclare package medium = medium,
    l=l[1],
    h=h[1],
    w=w[1],
    d=d[1],
    TStart=TStart,
    pStart=pStart,
    each xStart=XStart,
    maximalFlow=maximalFlow,
    pipe=pipe[1],
    nParallel=nParallel[1],
    useAlphaIn=false,
    nHT=1) annotation (Placement(transformation(extent={{-64,-10},{-44,10}})));
equation

  connect(simpleGasHeaterA2.ht[1], ht[1]) annotation (Line(points={{-54,4},{-26,
          4},{-26,33.3333},{0,33.3333}}, color={191,0,0}));
  connect(simpleGasHeaterA2.inlet, inlet) annotation (Line(points={{-62,0},{-76,
          0},{-76,1.77636e-15},{-100,1.77636e-15}}, color={159,159,223}));
  connect(simpleGasHeaterA2.outlet, flowSplit.inlet)
    annotation (Line(points={{-46,0},{-10,0}}, color={159,159,223}));
  connect(flowSplit.outlet1, simpleGasHeaterA.inlet) annotation (Line(points={{2,
          4},{12,4},{12,24},{20,24}}, color={159,159,223}));
  connect(simpleGasHeaterA.outlet, outletA) annotation (Line(points={{36,24},{50,
          24},{50,46},{66,46}}, color={159,159,223}));
  connect(flowSplit.outlet2, simpleGasHeaterB.inlet) annotation (Line(points={{2,
          -4},{2,-18},{28,-18},{28,-32}}, color={159,159,223}));
  connect(simpleGasHeaterB.outlet, outletB) annotation (Line(points={{44,-32},{56,
          -32},{56,-34},{66,-34}}, color={159,159,223}));
  connect(simpleGasHeaterA.ht[1], ht[2]) annotation (Line(points={{28,28},{10,28},
          {10,40},{0,40}}, color={191,0,0}));
  connect(simpleGasHeaterB.ht[1], ht[3]) annotation (Line(points={{36,-28},{18,
          -28},{18,46.6667},{0,46.6667}},
                                     color={191,0,0}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                                                Polygon(
          points={{40,60},{0,20},{-40,20},{-40,-20},{0,-20},{40,-60},{40,-20},
              {22,0},{40,20},{40,60}},
          lineColor={128,128,128},
          fillColor={159,159,223},
          fillPattern=FillPattern.Solid)}),                      Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end GasHeatedSplitter;

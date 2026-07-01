within OpenTEMPEST.Flow;
model gasHeaterAndJoiner

  import SI = Modelica.SIunits;

  replaceable package medium =
      Modelica.Media.Interfaces.PartialMedium                          annotation(choicesAllMatching = true);

  parameter SI.Temperature TStart=1000 "cell Starting inlet temperature for simulation"
    annotation (Dialog(tab="Initialization"));
  parameter SI.AbsolutePressure pStart=101325 "cell Starting Pressure"
    annotation (Dialog(tab="Initialization"));
  parameter SI.MassFraction XStart[medium.nX]=medium.reference_X annotation (Dialog(tab="Initialization"));

  parameter Boolean[3] pipe={true,false,true}
    "= true if flow through a pipe, false if between two parallel plates";
  parameter Integer[3] nParallel={1,30,1} "number of parallel flow ducts";
  parameter SI.Length[3] l={0.09,0.09,0.09} "lenght of flow channel";
  parameter SI.Length[3] h={0.01,0.01,0.01} "height of flow channel between two parallel plates";
  parameter SI.Length[3] w={0.04,0.04,0.04} "width of flow channel between two parallel plates";
  parameter SI.Length[3] d={0.02,0.02,0.02} "diameter of flow channel in a pipe";
  parameter SI.MassFlowRate maximalFlow = -1 "maximal mass flow rate";

  ThermoPower.Gas.FlangeA inletA(redeclare package Medium = medium) annotation (
      Placement(transformation(extent={{-100,32},{-72,60}}), iconTransformation(
          extent={{-80,20},{-40,60}})));
  ThermoPower.Gas.FlangeB outlet(redeclare package Medium = medium) annotation (
      Placement(transformation(extent={{52,-8},{80,20}}), iconTransformation(
          extent={{40,-20},{80,20}})));
  ThermoPower.Thermal.HT ht[3] annotation (Placement(transformation(extent={{-10,30},
            {10,50}}),     iconTransformation(extent={{-18,-20},{22,20}})));

  ThermoPower.Gas.FlangeA inletB(redeclare package Medium = medium) annotation (
      Placement(transformation(extent={{-100,-46},{-72,-18}}), iconTransformation(
          extent={{-80,-60},{-40,-20}})));
  OpenTEMPEST.BOP.GasConditioning.SimpleGasHeater simpleGasHeaterA(
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
    nHT=1) annotation (Placement(transformation(extent={{-54,18},{-34,38}})));
  ThermoPower.Gas.FlowJoin flowJoin(redeclare package Medium = medium)
    annotation (Placement(transformation(extent={{-20,-4},{0,16}})));
  OpenTEMPEST.BOP.GasConditioning.SimpleGasHeater simpleGasHeaterB(
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
    nHT=1) annotation (Placement(transformation(extent={{-48,-40},{-28,-20}})));
  OpenTEMPEST.BOP.GasConditioning.SimpleGasHeater simpleGasHeaterOut(
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
    nHT=1) annotation (Placement(transformation(extent={{12,0},{32,20}})));
equation

  connect(simpleGasHeaterA.ht[1], ht[1]) annotation (Line(points={{-44,32},{-22,
          32},{-22,33.3333},{0,33.3333}}, color={191,0,0}));
  connect(simpleGasHeaterA.inlet, inletA) annotation (Line(points={{-52,28},{-64,
          28},{-64,46},{-86,46}}, color={159,159,223}));
  connect(flowJoin.inlet1, simpleGasHeaterA.outlet) annotation (Line(points={{-16,
          10},{-26,10},{-26,28},{-36,28}}, color={159,159,223}));
  connect(simpleGasHeaterB.inlet, inletB) annotation (Line(points={{-46,-30},{-62,
          -30},{-62,-32},{-86,-32}}, color={159,159,223}));
  connect(simpleGasHeaterB.outlet, flowJoin.inlet2) annotation (Line(points={{-30,
          -30},{-30,-19},{-16,-19},{-16,2}}, color={159,159,223}));
  connect(simpleGasHeaterB.ht[1], ht[2]) annotation (Line(points={{-38,-26},{-20,
          -26},{-20,40},{0,40}}, color={191,0,0}));
  connect(simpleGasHeaterOut.inlet, flowJoin.outlet) annotation (Line(points={{14,
          10},{6,10},{6,6},{-4,6}}, color={159,159,223}));
  connect(simpleGasHeaterOut.outlet, outlet) annotation (Line(points={{30,10},{52,
          10},{52,6},{66,6}}, color={159,159,223}));
  connect(simpleGasHeaterOut.ht[1], ht[3]) annotation (Line(points={{22,14},{12,
          14},{12,46.6667},{0,46.6667}}, color={191,0,0}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                                                Polygon(
          points={{-40,60},{0,20},{40,20},{40,-20},{0,-20},{-40,-60},{-40,-20},{-22,
              0},{-40,20},{-40,60}},
          lineColor={128,128,128},
          fillColor={159,159,223},
          fillPattern=FillPattern.Solid)}),                      Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end gasHeaterAndJoiner;

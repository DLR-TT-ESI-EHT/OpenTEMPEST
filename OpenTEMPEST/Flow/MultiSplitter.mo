within OpenTEMPEST.Flow;
model MultiSplitter

  parameter Integer nOutlets(min = 2);
  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
    annotation(choicesAllMatching = true);

  ThermoPower.Gas.FlangeB outlet[nOutlets](redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{60,0},{80,20}}),
        iconTransformation(extent={{40,-20},{80,20}})));
  ThermoPower.Gas.FlangeA inlet(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-60,0},{-40,20}}),
        iconTransformation(extent={{-80,-20},{-40,20}})));
  ThermoPower.Gas.FlowSplit flowSplit[nOutlets-1](redeclare package Medium = Medium,
    allowFlowReversal=true,
    checkFlowDirection=true,
    rev_inlet=true,
    rev_outlet1=true,
    rev_outlet2=true);
equation

  connect(flowSplit[1].inlet, inlet);
  connect(flowSplit[1].outlet1, outlet[1]);

  for i in 2:nOutlets-1 loop
    connect(flowSplit[i-1].outlet2, flowSplit[i].inlet);
    connect(flowSplit[i].outlet1, outlet[i]);
  end for;

  connect(flowSplit[nOutlets - 1].outlet2, outlet[nOutlets]);

    annotation (Placement(transformation(extent={{-44,-10},{-24,10}})), Icon(
        graphics={
          Rectangle(
          extent={{-40,32},{40,-30}},
          fillColor={170,213,255},
          fillPattern=FillPattern.Solid,
          pattern=LinePattern.None)}));
end MultiSplitter;

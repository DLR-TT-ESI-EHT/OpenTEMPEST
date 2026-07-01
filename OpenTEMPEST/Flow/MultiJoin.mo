within OpenTEMPEST.Flow;
model MultiJoin

  parameter Integer nInlets(min = 2);
  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
    annotation(choicesAllMatching = true);

  ThermoPower.Gas.FlangeB outlet(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{40,-10},{60,10}}),
        iconTransformation(extent={{40,-20},{80,20}})));
  ThermoPower.Gas.FlangeA inlet[nInlets](redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-60,0},{-40,20}}),
        iconTransformation(extent={{-80,-20},{-40,20}})));
  ThermoPower.Gas.FlowJoin flowJoin[nInlets-1](redeclare package Medium = Medium,
    allowFlowReversal=true,
    checkFlowDirection=true);
//     rev_inlet=true,
//     rev_outlet1=true,
//     rev_outlet2=true,
equation

  connect(flowJoin[1].inlet1, inlet[1]);

  for i in 2:nInlets-1 loop
    connect(flowJoin[i-1].inlet2, flowJoin[i].outlet);
    connect(flowJoin[i].inlet1, inlet[i]);
  end for;

  connect(flowJoin[nInlets - 1].inlet2, inlet[nInlets]);
  connect(flowJoin[1].outlet, outlet);

    annotation (Placement(transformation(extent={{-44,-10},{-24,10}})), Icon(
        graphics={
          Rectangle(
          extent={{-40,32},{40,-30}},
          fillColor={170,213,255},
          fillPattern=FillPattern.Solid,
          pattern=LinePattern.None), Text(
          extent={{-40,32},{40,-30}},
          lineColor={28,108,200},
          textString="MJ")}));
end MultiJoin;

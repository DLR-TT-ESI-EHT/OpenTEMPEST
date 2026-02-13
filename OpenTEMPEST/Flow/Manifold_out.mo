within OpenTEMPEST.Flow;
model Manifold_out
  function positiveMax
    extends Modelica.Icons.Function;
    input Real x;
    output Real y;
  algorithm
    y := max(x, 1e-10);
  end positiveMax;

  import Modelica.Constants;
  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium annotation (
    choicesAllMatching);
  // Ports
  parameter Integer nPorts_a = 0
    "Number of outlet ports (mass is distributed evenly between the outlet ports"                       annotation (
    Dialog(connectorSizing = true));
  Modelica.Fluid.Interfaces.FluidPorts_a ports_a[nPorts_a](redeclare package
      Medium =                                                                        Medium) annotation (
    Placement(transformation(extent = {{-36, -40}, {-16, 40}})));
  Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = Medium) annotation (
    Placement(transformation(extent = {{40, -10}, {60, 10}}), iconTransformation(extent = {{40, -10}, {60, 10}})));
  Medium.MassFraction ports_a_Xi_inStream[nPorts_a, Medium.nXi]
    "inStream mass fractions at ports_b";
  Medium.ExtraProperty ports_a_C_inStream[nPorts_a, Medium.nC]
    "inStream extra properties at ports_b";
equation
  // Only one connection allowed to a port to avoid unwanted ideal mixing
  for i in 1:nPorts_a loop
    assert(cardinality(ports_a[i]) <= 1, "
each ports_a[i] of boundary shall at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections, which is usually not the intention
of the modeller. Increase nPorts_b to add an additional port.
        ");
  end for;
  // mass and momentum balance
  0 = sum(ports_a.m_flow) + port_b.m_flow;
  ports_a.p = fill(port_b.p, nPorts_a);
  // mixing at port_b
  port_b.h_outflow = sum({positiveMax(ports_a[j].m_flow) * inStream(ports_a[j].h_outflow) for j in 1:nPorts_a}) / sum({positiveMax(ports_a[j].m_flow) for j in 1:nPorts_a});
  for j in 1:nPorts_a loop
    // expose stream values from port_b to ports_a
    ports_a[j].h_outflow = actualStream(port_b.h_outflow);
    ports_a[j].Xi_outflow = actualStream(port_b.Xi_outflow);
    ports_a[j].C_outflow = actualStream(port_b.C_outflow);
    ports_a_Xi_inStream[j, :] = inStream(ports_a[j].Xi_outflow);
    ports_a_C_inStream[j, :] = inStream(ports_a[j].C_outflow);
  end for;
  for i in 1:Medium.nXi loop
    port_b.Xi_outflow[i] = positiveMax(ports_a.m_flow) * ports_a_Xi_inStream[:, i] / sum(positiveMax(ports_a.m_flow));
  end for;
  for i in 1:Medium.nC loop
    port_b.C_outflow[i] = positiveMax(ports_a.m_flow) * ports_a_C_inStream[:, i] / sum(positiveMax(ports_a.m_flow));
  end for;
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics),
    Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Line(points = {{-26, 26}, {50, 0}}, color = {0, 128, 255}, thickness = 1), Line(points = {{-28, -28}, {52, 0}}, color = {0, 128, 255}, thickness = 1), Line(points = {{-28, 0}, {52, 0}}, color = {0, 128, 255}, thickness = 1)}));
end Manifold_out;

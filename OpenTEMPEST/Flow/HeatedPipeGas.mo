within OpenTEMPEST.Flow;
model HeatedPipeGas

  import SI = Modelica.SIunits;

  replaceable package Medium = OpenTEMPEST.Medium.Fuel_CH4
                                                        constrainedby
    Modelica.Media.Interfaces.PartialMedium annotation (choicesAllMatching=true);

  Medium.BaseProperties gasInH;
  Medium.BaseProperties gasInC;

  ThermoPower.Gas.FlangeA inlet(redeclare package Medium = Medium) annotation (
      Placement(transformation(extent={{-100,-10},{-80,10}}),
        iconTransformation(extent={{-100,-20},{-60,20}})));
  ThermoPower.Gas.FlangeB outlet(redeclare package Medium = Medium) annotation (
     Placement(transformation(extent={{80,-10},{100,10}}), iconTransformation(
          extent={{60,-30},{100,10}})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort
    annotation (Placement(transformation(extent={{-10,-110},{10,-90}})));

protected
  SI.MassFraction xIn[Medium.nXi] " MAss flow rates entering the reactors";
  SI.MassFraction xOut[Medium.nXi]
    "Mass flow rates from reactor to mixing gas node";
  SI.SpecificEnthalpy hIn;
  SI.SpecificEnthalpy hOut;

equation

  // Balances
  0 = inlet.m_flow + outlet.m_flow;
  //mass Balamce
  0 = inlet.m_flow*hIn + outlet.m_flow*hOut;
  xIn = gas.Xi;

  //connecting variables with connectors
  inlet.p = gas.p;
  outlet.p = gas.p;
  hIn = inStream(inlet.h_outflow);
  hOut = outlet.h_outflow;
  xIn = inStream(inlet.Xi_outflow);
  xOut = outlet.Xi_outflow;
  hOut = gas.h;
  xOut = gas.Xi;

  inlet.h_outflow = 0;
  inlet.Xi_outflow = Medium.reference_X[1:Medium.nXi];

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end HeatedPipeGas;

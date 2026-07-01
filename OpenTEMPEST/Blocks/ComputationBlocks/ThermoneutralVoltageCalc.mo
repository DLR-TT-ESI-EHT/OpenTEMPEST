within OpenTEMPEST.Blocks.ComputationBlocks;
block ThermoneutralVoltageCalc
  extends Modelica.Blocks.Icons.Block;

  Modelica.Blocks.Interfaces.RealInput in_T
    annotation (Placement(transformation(extent={{-140,0},{-100,40}})));

  Modelica.Blocks.Interfaces.RealOutput out_UTN
    annotation (Placement(transformation(extent={{100,30},{120,50}}),
        iconTransformation(extent={{100,30},{120,50}})));

  Modelica.Blocks.Interfaces.RealInput in_YH2O "Current input"
    annotation (Placement(transformation(extent={{-140,40},{-100,80}})));
  Modelica.Blocks.Interfaces.RealInput in_YCO2 "Current input"
    annotation (Placement(transformation(extent={{-140,-40},{-100,0}})));

  Modelica.Blocks.Interfaces.RealOutput out_UId
    annotation (Placement(transformation(extent={{100,-10},{120,10}}),
        iconTransformation(extent={{100,-10},{120,10}})));
  Modelica.Blocks.Interfaces.RealInput in_YH2 "Current input"
    annotation (Placement(transformation(extent={{-140,-80},{-100,-40}})));
  Modelica.Blocks.Interfaces.RealOutput out_ASR annotation (Placement(
        transformation(extent={{100,-50},{120,-30}}), iconTransformation(extent=
           {{100,-60},{120,-40}})));

equation

  out_UTN = -((in_YH2O/(in_YH2O + in_YCO2))*(0.003466165*in_T^2 - 13.28655501*
    in_T - 238033.4485) + (in_YCO2/(in_YH2O + in_YCO2))*(0.001884316*in_T^2 - 0.096511194
    *in_T - 284401.7366))/2/Modelica.Constants.F;

  out_UId = (-(0.005154497*in_T^2 + 45.41840517*in_T - 243122.6741) - Modelica.Constants.R
    *in_T*log(0.5*(in_YH2O + in_YH2O*(1 - 0.7))/0.5/(2*in_YH2 + in_YH2O*(0.7))/(
    0.21)^0.5/1.01325^0.5))/Modelica.Constants.F/2;

  out_ASR = 728*exp((in_T - 273.15)*(-0.00829));

  annotation(choicesAllMatching = true,
        Documentation(revisions="<html>
<ul>
<li><i>28 Okt 2021</i> by <a href=\"mailto:Marius.Tomberg@dlr.de\">Marius Tomberg</a>:<br>Created </li>
<li>12 Aug 2022 by Marius Tomberg:<br>Calculation of ideal voltage and ASR added</li>
</ul>
</html>"),
    Diagram(coordinateSystem(extent={{-100,-100},{100,100}})),
    Icon(coordinateSystem(extent={{-100,-100},{100,100}}), graphics={
        Text(
          extent={{-100,80},{-60,40}},
          lineColor={28,108,200},
          textString="YH2O"),
        Text(
          extent={{-100,40},{-60,0}},
          lineColor={28,108,200},
          textString="T"),
        Text(
          extent={{-100,0},{-60,-40}},
          lineColor={28,108,200},
          textString="YCO2"),
        Text(
          extent={{-100,-40},{-60,-80}},
          lineColor={28,108,200},
          textString="YH2"),
        Text(
          extent={{60,60},{100,20}},
          lineColor={28,108,200},
          textString="UTN"),
        Text(
          extent={{60,20},{100,-20}},
          lineColor={28,108,200},
          textString="UId"),
        Text(
          extent={{60,-20},{100,-60}},
          lineColor={28,108,200},
          textString="ASR")}));
end ThermoneutralVoltageCalc;

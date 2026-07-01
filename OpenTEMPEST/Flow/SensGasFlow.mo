within OpenTEMPEST.Flow;
model SensGasFlow "Sensor for flow parameters of gas flows"
  extends ThermoPower.Icons.Gas.SensThrough;

  parameter Boolean mfOutput = true "Sensor has output for mass flow";
  parameter Boolean pOutput = true "Sensor has output for pressure flow";
  parameter Boolean hOutput = true "Sensor has output for enthalpy flow";
  parameter Boolean XOutput = true "Sensor has output for mass fractions";
  parameter Boolean YOutput = true "Sensor has output for molar fractions";
  parameter Boolean HfOutput = true "Sensor has output for enthalpy flow";

  //Modelica.Media.Interfaces.PartialMedium
  replaceable package Medium = OpenTEMPEST.Medium.Fuel_CH4
                                                        constrainedby
    Modelica.Media.Interfaces.PartialMedium
    annotation(choicesAllMatching = true);
  parameter Boolean allowFlowReversal=system.allowFlowReversal
    "= true to allow flow reversal, false restricts to design direction"
    annotation(Evaluate=true);

  Modelica.SIunits.MolarMass M = 1/sum(outlet.Xi_outflow[:]./Medium.MMX[:]);

  outer ThermoPower.System system "System wide properties";
  ThermoPower.Gas.FlangeA inlet(redeclare package Medium = Medium, m_flow(min=
          if allowFlowReversal then -Modelica.Constants.inf else 0)) annotation (
     Placement(transformation(extent={{-80,-60},{-40,-20}}, rotation=0)));
  ThermoPower.Gas.FlangeB outlet(redeclare package Medium = Medium, m_flow(max=
          if allowFlowReversal then +Modelica.Constants.inf else 0)) annotation (
     Placement(transformation(extent={{40,-60},{80,-20}}, rotation=0)));
  Modelica.Blocks.Interfaces.RealOutput Hf if HfOutput annotation (Placement(
        transformation(extent={{56,50},{76,70}}, rotation=0),
        iconTransformation(extent={{60,60},{80,80}})));
  Modelica.Blocks.Interfaces.RealOutput mf if mfOutput annotation (Placement(
        transformation(extent={{56,30},{76,50}}, rotation=0),
        iconTransformation(extent={{60,100},{80,120}})));
  Modelica.Blocks.Interfaces.RealOutput p if pOutput annotation (Placement(
        transformation(extent={{56,68},{76,88}}, rotation=0),
        iconTransformation(extent={{60,40},{80,60}})));
  Modelica.Blocks.Interfaces.RealOutput h if hOutput annotation (Placement(
        transformation(extent={{56,14},{76,34}}, rotation=0),
        iconTransformation(extent={{60,80},{80,100}})));
  Modelica.Blocks.Interfaces.RealOutput x[Medium.nXi] if XOutput
                                                       annotation (Placement(
        transformation(extent={{56,-4},{76,16}}, rotation=0),
        iconTransformation(extent={{60,20},{80,40}})));
  Modelica.Blocks.Interfaces.RealOutput y[Medium.nXi] if YOutput annotation (
      Placement(transformation(extent={{56,-18},{76,2}}, rotation=0),
        iconTransformation(extent={{60,0},{80,20}})));
protected
  Modelica.Blocks.Sources.RealExpression realExpression_p(y=if pOutput then inlet.p else -1)
    annotation (Placement(transformation(extent={{-8,68},{12,88}})));
  Modelica.Blocks.Sources.RealExpression realExpression_Hf(y=if HfOutput then inlet.m_flow*outlet.h_outflow else -1)
    annotation (Placement(transformation(extent={{-10,50},{10,70}})));
  Modelica.Blocks.Sources.RealExpression realExpression_mf(y=if mfOutput then
        inlet.m_flow else -1)
    annotation (Placement(transformation(extent={{-10,30},{10,50}})));
  Modelica.Blocks.Sources.RealExpression realExpression_h(y=if hOutput then
        outlet.h_outflow else -1)
    annotation (Placement(transformation(extent={{-10,12},{10,32}})));
  Modelica.Blocks.Sources.RealExpression realExpression_X[Medium.nXi](y=outlet.Xi_outflow[
        :])
    annotation (Placement(transformation(extent={{-12,-6},{8,14}})));
  Modelica.Blocks.Sources.RealExpression realExpression_Y[Medium.nXi](y=outlet.Xi_outflow[
        :]./Medium.MMX[:]*M) annotation (Placement(transformation(extent={{-14,-20},{6,0}})));
equation
  inlet.m_flow + outlet.m_flow = 0 "Mass balance";
  inlet.p = outlet.p "Momentum balance";

  // Energy balance
  inlet.h_outflow = inStream(outlet.h_outflow);
  inStream(inlet.h_outflow) = outlet.h_outflow;

  // Independent composition mass balances
  inlet.Xi_outflow = inStream(outlet.Xi_outflow);
  inStream(inlet.Xi_outflow) = outlet.Xi_outflow;

  // Sensor output
  for i in 1:Medium.nXi loop
    connect(realExpression_X[i].y,x [i])
    annotation (Line(points={{9,4},{38,4},{38,6},{66,6}}, color={0,0,127}));
    connect(realExpression_Y[i].y,y [i]) annotation (Line(points={{7,-10},{34,-10},
          {34,-8},{66,-8}}, color={0,0,127}));
  end for;

  connect(realExpression_p.y, p)
    annotation (Line(points={{13,78},{66,78}}, color={0,0,127}));
  connect(realExpression_Hf.y, Hf)
    annotation (Line(points={{11,60},{66,60}}, color={0,0,127}));
  connect(realExpression_mf.y, mf) annotation (Line(points={{11,40},{34,40},{34,
          40},{66,40}}, color={0,0,127}));
  connect(realExpression_h.y, h) annotation (Line(points={{11,22},{36,22},{36,24},
          {66,24}}, color={0,0,127}));

  annotation (
    Documentation(revisions="<html>
<ul>
<li><i>20 Dec 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
      Adapted to Modelica.Media.</li>
<li><i>5 Mar 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</html>",
        info="<html>
<p>This component can be inserted in a hydraulic circuit to measure the flowrate of the fluid flowing through it.
<p>Flow reversal is supported.
</html>"),
    Icon(graphics={Text(
          extent={{-34,82},{34,40}},
          lineColor={0,0,0},
          textString="Sens")}));
end SensGasFlow;

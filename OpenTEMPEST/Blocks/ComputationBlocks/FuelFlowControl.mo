within OpenTEMPEST.Blocks.ComputationBlocks;
model FuelFlowControl
  "Flow rate source for gas flows in electrochemical processes"

  replaceable package Medium = OpenTEMPEST.Medium.Fuel_CH4;
  extends Modelica.Blocks.Icons.Block;

  parameter Real compositionNom[Medium.nX]=Medium.reference_X
    "Nominal gas composition";
  parameter Real FormingGasComposition[Medium.nX] = {0.00377,0.00001,0.00001,0.00001,0.00001,0.99623} "Secondary or Forming Gas composition (default 5/95)";
  parameter Medium.MassFraction Xnom[Medium.nX]=if molarInputComposition then Medium.moleToMassFractions(Ynom,Medium.MMX) else compositionNom
    "Nominal gas composition" annotation (Dialog(tab="Calculated parameter"));
  parameter Medium.MoleFraction Ynom[Medium.nX]=if molarInputComposition then compositionNom else Medium.massToMoleFractions(Xnom,Medium.MMX)
    "Nominal gas composition" annotation (Dialog(tab="Calculated parameter"));
  parameter Boolean molarInputComposition = false "Input compostions are molar instead of mass" annotation(choices(checkBox=true));
  parameter Modelica.SIunits.ElectricCurrent I0=0 "Nominal electrical current";
  parameter Real RC0 = 0.7 "Nominal reactant conversion";
  parameter Integer nCells "Number of cells that are fed by this source";
  parameter Modelica.SIunits.Time tSmooth(start=1) "Base time constant for PT1 output smoothing";
  parameter Modelica.SIunits.MassFlowRate mfStart "Initial value for the mass flow";
  parameter Modelica.SIunits.ElectricCurrent IStart "Initial value for the current";

  parameter Modelica.SIunits.ElectricCurrent minCurrentForDosing = 15
    "Minimal absolute current to calculate flows based on current" annotation(Dialog(group="Minimal flow"));
  parameter Boolean doseFGUnderIMin = false
    "Use forming gas (FormingGasComposition) instead of given compostion if flow is below minimal current" annotation(Dialog(group="Minimal flow"), choices(checkBox=true));
  parameter Modelica.SIunits.MassFlowRate mfFGPerCell = 6/1000/1440
    "Forming gas (FormingGasComposition) mass flow per cell at 0 A if activated" annotation(Dialog(group="Minimal flow"));
  parameter Boolean use_in_I = false
    "Use connector input for the electrical current" annotation(Dialog(group="External inputs"), choices(checkBox=true));
  parameter Boolean use_in_RC = false
    "Use connector input for the reactant conversion" annotation(Dialog(group="External inputs"), choices(checkBox=true));
  parameter Boolean use_in_composition = false
    "Use connector input for the composition"                                  annotation(Dialog(group="External inputs"), choices(checkBox=true));
  outer ThermoPower.System system "System wide properties";

  Medium.MassFlowRate w "Nominal mass flow rate";
  Medium.MoleFraction Y[Medium.nX];
  Medium.MoleFraction YFG[Medium.nX];
  Medium.MassFraction X[Medium.nX];
  Medium.MassFraction XFG[Medium.nX];
  Medium.MassFlowRate wFG "FG mass flow rate";
  Medium.MolarMass MM = sum(Y*Medium.MMX);
  Medium.MolarMass MMFG = sum(YFG*Medium.MMX);

  parameter Real zFC[Medium.nX] = {2,8,0,2,0,0};
  parameter Real zEC[Medium.nX] = {0,-4,2,0,2,0};
  Real out_z;

  Modelica.Blocks.Interfaces.RealInput in_I if use_in_I annotation (Placement(
        transformation(
        origin={-60,50},
        extent={{-10,-10},{10,10}},
        rotation=270), iconTransformation(extent={{-140,40},{-100,80}},
          rotation=0)));
  Modelica.Blocks.Interfaces.RealInput in_RC if use_in_RC annotation (Placement(
        transformation(
        origin={-40,50},
        extent={{-10,-10},{10,10}},
        rotation=270), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=0,
        origin={-120,0})));
  Modelica.Blocks.Interfaces.RealInput in_composition[Medium.nX]
                                                       if use_in_composition annotation (
      Placement(transformation(
        origin={28,50},
        extent={{-10,-10},{10,10}},
        rotation=270), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=0,
        origin={-120,-60})));
  Modelica.Blocks.Interfaces.RealOutput out_y[Medium.nX] annotation (Placement(
        transformation(extent={{10,-66},{30,-46}}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=0,
        origin={120,-80})));

  Modelica.Blocks.Interfaces.RealOutput out_I(start=IStart)
    annotation (Placement(transformation(extent={{100,60},{140,100}}),
        iconTransformation(extent={{100,60},{140,100}})));
  Modelica.Blocks.Interfaces.RealOutput out_mf(start=mfStart)
    annotation (Placement(transformation(extent={{100,20},{140,60}}),
        iconTransformation(extent={{100,20},{140,60}})));
  Modelica.Blocks.Interfaces.RealOutput out_x[Medium.nX] annotation (Placement(
        transformation(extent={{46,-46},{66,-26}}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=0,
        origin={120,-40})));
  Modelica.Blocks.Interfaces.RealOutput out_mfi[Medium.nX] annotation (
      Placement(transformation(extent={{100,-58},{128,-30}}),
        iconTransformation(extent={{100,-20},{140,20}})));
protected
  Modelica.Blocks.Interfaces.RealInput in_I_internal;
  Modelica.Blocks.Interfaces.RealInput in_RC_internal;
  Modelica.Blocks.Interfaces.RealInput in_composition_internal[Medium.nX];

  Modelica.SIunits.ElectricCurrent IForDosing;

equation

  if in_I_internal >= 0 then
    out_z = sum(zFC * Y);
    IForDosing = if in_I_internal < minCurrentForDosing and not doseFGUnderIMin then minCurrentForDosing else in_I_internal;
  else
    out_z = sum(zEC * Y);
    IForDosing = if abs(in_I_internal) < minCurrentForDosing and not doseFGUnderIMin then -minCurrentForDosing else in_I_internal;
  end if;

  wFG = if abs(in_I_internal) < minCurrentForDosing and doseFGUnderIMin then (1-abs(in_I_internal)/minCurrentForDosing) * mfFGPerCell*nCells else 0;
  w = abs(IForDosing)*nCells/out_z/Modelica.Constants.F/in_RC_internal * MM;

  if not use_in_I then
    in_I_internal = I0 "Flow rate set by parameter";
  end if;

  if not use_in_RC then
    in_RC_internal = RC0 "Flow rate set by parameter";
  end if;

  if molarInputComposition then
    Y = in_composition_internal[1:Medium.nXi];
    X = Medium.moleToMassFractions(Y, Medium.MMX);
    YFG = FormingGasComposition;
    XFG = Medium.moleToMassFractions(YFG, Medium.MMX);
  else
    X = in_composition_internal[1:Medium.nXi];
    Y = Medium.massToMoleFractions(X, Medium.MMX);
    XFG = FormingGasComposition;
    YFG = Medium.massToMoleFractions(XFG, Medium.MMX);
  end if;

  if not use_in_composition then
    in_composition_internal = compositionNom "Composition set by parameter";
  end if;

  out_y = (Y*w/MM+YFG*wFG/MMFG)/(w/MM+wFG/MMFG+Modelica.Constants.small) "outlet molar fraction";
  out_x = Medium.moleToMassFractions(out_y, Medium.MMX);
  if der(out_mf)>0 then
    der(out_I) = (in_I_internal - out_I)/2/tSmooth;
    der(out_mf) = (w+wFG - out_mf)/tSmooth;
  else
    der(out_I) = (in_I_internal - out_I)/tSmooth;
    der(out_mf) = (w+wFG - out_mf)/2/tSmooth;
  end if;
  out_mfi = out_mf*out_x;

  // Connect protected connectors to public conditional connectors
  connect(in_I, in_I_internal);
  connect(in_RC, in_RC_internal);
  connect(in_composition, in_composition_internal);

    annotation(choicesAllMatching = true,
              Documentation(info="<html>
<p><b>Modelling options</b></p>
<p>The actual gas used in the component is determined by the replaceable <tt>Medium</tt> package. In the case of multiple component, variable composition gases, the nominal gas composition is given by <tt>Xnom</tt>,whose default value is <tt>Medium.reference_X</tt> .
<p>If <tt>G</tt> is set to zero, the flowrate source is ideal; otherwise, the outgoing flowrate decreases proportionally to the outlet pressure.</p>
<p>If the <tt>in_w0</tt> connector is wired, then the source massflowrate is given by the corresponding signal, otherwise it is fixed to <tt>w0</tt>.</p>
<p>If the <tt>in_T</tt> connector is wired, then the source temperature is given by the corresponding signal, otherwise it is fixed to <tt>T</tt>.</p>
<p>If the <tt>in_X</tt> connector is wired, then the source massfraction is given by the corresponding signal, otherwise it is fixed to <tt>Xnom</tt>.</p>
</html>",
        revisions="<html>
        <ul>
<li><i>Jan 2023</i>
by <a href=\"mailto:Marius.Tomberg@dlr.de\">Marius Tomberg</a>:<br>
Created 

</ul>
</html>"),
    Diagram(coordinateSystem(extent={{-100,-100},{100,100}})),
    Icon(coordinateSystem(extent={{-100,-100},{100,100}}), graphics={
        Text(
          extent={{-98,80},{-58,40}},
          lineColor={28,108,200},
          textString="I"),
        Text(
          extent={{-98,20},{-58,-20}},
          lineColor={28,108,200},
          textString="RC"),
        Text(
          extent={{-98,-40},{-58,-80}},
          lineColor={28,108,200},
          textString="x/yi"),
        Text(
          extent={{58,100},{98,60}},
          lineColor={28,108,200},
          textString="I"),
        Text(
          extent={{58,60},{98,20}},
          lineColor={28,108,200},
          textString="mf"),
        Text(
          extent={{58,20},{98,-20}},
          lineColor={28,108,200},
          textString="mfi"),
        Text(
          extent={{58,-20},{98,-60}},
          lineColor={28,108,200},
          textString="xi"),
        Text(
          extent={{58,-60},{98,-100}},
          lineColor={28,108,200},
          textString="yi")}));
end FuelFlowControl;

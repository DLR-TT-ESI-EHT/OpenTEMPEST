within OpenTEMPEST.Flow;
model SourceMassFlowEC
  "Flow rate source for gas flows in electrochemical processes"
  extends ThermoPower.Icons.Gas.SourceW;
  replaceable package Medium = OpenTEMPEST.Medium.Fuel_CH4;
                                                         //Modelica.Media.Interfaces.PartialMedium
  Medium.BaseProperties gas(
    p(start=p0),
    T(start=T),
    Xi(start=compositionNom[1:Medium.nXi]));

  Medium.BaseProperties FormingGas(
    p(start=p0),
    T(start=T),
    Xi(start=FormingGasComposition));

  parameter Medium.AbsolutePressure p0=101325 "Nominal pressure";
  parameter Medium.Temperature T=300 "Nominal temperature";
  parameter Medium.Temperature TSteamMin=400 "Minimum temperature of the steam";
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
  parameter ThermoPower.Units.HydraulicConductance G=0 "HydraulicConductance";
  parameter Boolean allowFlowReversal=system.allowFlowReversal
    "= true to allow flow reversal, false restricts to design direction"
    annotation(Evaluate=true);
  parameter Integer nCells "Number of cells that are fed by this source";

  parameter Modelica.SIunits.ElectricCurrent minCurrentForDosing = 15
    "Minimal absolute current to calculate flows based on current" annotation(Dialog(group="Minimal flow"));
  parameter Boolean doseFGUnderIMin = false
    "Use forming gas (FormingGasComposition) instead of given compostion if flow is below minimal current" annotation(Dialog(group="Minimal flow"), choices(checkBox=true));
  parameter Modelica.SIunits.MassFlowRate mfFGPerCell = 6/1000/1440
    "Forming gas (FormingGasComposition) mass flow per cell at 0 A if activated" annotation(Dialog(group="Minimal flow"));
  parameter Modelica.SIunits.MassFraction minH2OXiforIout = 0.05 "Minimum H2/H2O (FC/EC) Massfraction for current signal to pass on";
  parameter Boolean use_in_I = false
    "Use connector input for the electrical current" annotation(Dialog(group="External inputs"), choices(checkBox=true));
  parameter Boolean use_in_RC = false
    "Use connector input for the reactant conversion" annotation(Dialog(group="External inputs"), choices(checkBox=true));
  parameter Boolean use_in_T = false
    "Use connector input for the temperature"                                  annotation(Dialog(group="External inputs"), choices(checkBox=true));
  parameter Boolean use_in_composition = false
    "Use connector input for the composition"                                  annotation(Dialog(group="External inputs"), choices(checkBox=true));
  parameter Boolean use_in_mfFGPerCell = false
    "Use connector input for mfFGPerCell"                                  annotation(Dialog(group="External inputs"), choices(checkBox=true));
  outer ThermoPower.System system "System wide properties";

  Medium.MassFlowRate w "Nominal mass flow rate";
  Medium.MoleFraction Y[Medium.nX];
  Medium.MoleFraction Yfg[Medium.nX];
  Medium.MassFlowRate wFG "FG mass flow rate";

  Real zFC[Medium.nX] = {2,8,0,2,0,0};
  Real zEC[Medium.nX] = {0,-2/in_RC_internal,2,0,2,0}; // Actual RC will vary based on extent of MSR reaction by -2/RC*Percent_Extent_of_MSR

  ThermoPower.Gas.FlangeB flange(redeclare package Medium = Medium, m_flow(max=
          if allowFlowReversal then +Modelica.Constants.inf else 0))
    annotation (Placement(transformation(extent={{80,-20},{120,20}}, rotation=0)));
  Modelica.Blocks.Interfaces.RealInput in_I if use_in_I annotation (Placement(
        transformation(
        origin={-60,50},
        extent={{-10,-10},{10,10}},
        rotation=270)));
  Modelica.Blocks.Interfaces.RealInput in_RC if use_in_RC annotation (Placement(
        transformation(
        origin={-40,50},
        extent={{-10,-10},{10,10}},
        rotation=270), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-34,50})));
  Modelica.Blocks.Interfaces.RealInput in_T if use_in_T annotation (Placement(
        transformation(
        origin={-8,50},
        extent={{10,-10},{-10,10}},
        rotation=90), iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={-8,50})));
  Modelica.Blocks.Interfaces.RealInput in_composition[Medium.nX] if use_in_composition annotation (
      Placement(transformation(
        origin={60,50},
        extent={{-10,-10},{10,10}},
        rotation=270)));
  Modelica.Blocks.Interfaces.RealInput in_mfFGPerCell if use_in_mfFGPerCell   annotation (Placement(
        transformation(
        origin={12,50},
        extent={{10,-10},{-10,10}},
        rotation=90), iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={22,50})));
  Modelica.Blocks.Interfaces.RealOutput out_z annotation (Placement(
        transformation(extent={{-46,-66},{-26,-46}}), iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={-60,-50})));
  Modelica.Blocks.Interfaces.RealOutput out_y[Medium.nX] annotation (Placement(
        transformation(extent={{10,-66},{30,-46}}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-20,-50})));
  Modelica.Blocks.Interfaces.RealOutput out_I  annotation (Placement(transformation(extent={{100,-56},{128,-28}})));

protected
  Modelica.Blocks.Interfaces.RealInput in_I_internal;
  Modelica.Blocks.Interfaces.RealInput in_RC_internal;
  Modelica.Blocks.Interfaces.RealInput in_T_internal;
  Modelica.Blocks.Interfaces.RealInput in_composition_internal[Medium.nX];
  Modelica.Blocks.Interfaces.RealInput in_mfFGPerCell_internal;

  Modelica.SIunits.ElectricCurrent IForDosing;

equation

  if G > 0 then
    flange.m_flow = -w-wFG + (flange.p - p0)*G;
  else
    flange.m_flow = -w-wFG;
  end if;

  if in_I_internal >= 0 then
    out_z = sum(zFC * Y);
    IForDosing = if in_I_internal < minCurrentForDosing and not doseFGUnderIMin then minCurrentForDosing else in_I_internal;
  else
    out_z = sum(zEC * Y);
    IForDosing = if abs(in_I_internal) < minCurrentForDosing and not doseFGUnderIMin then -minCurrentForDosing else in_I_internal;
  end if;

  wFG = if abs(in_I_internal) < minCurrentForDosing and doseFGUnderIMin then (1-abs(in_I_internal)/minCurrentForDosing) * in_mfFGPerCell_internal*nCells else 0;
  w = abs(IForDosing)*nCells/out_z/Modelica.Constants.F/in_RC_internal * gas.MM;

  if not use_in_I then
    in_I_internal = I0 "Flow rate set by parameter";
  end if;

  if not use_in_RC then
    in_RC_internal = RC0 "Flow rate set by parameter";
  end if;

  gas.T = Medium.T_hX((Medium.h_TX(in_T_internal,{gas.Xi[1],gas.Xi[2],gas.Xi[3],gas.Xi[4],0,gas.Xi[6]}/(1-gas.Xi[5]))*(1-gas.Xi[5])+
            Medium.h_TX(max(in_T_internal,TSteamMin),{0,0,0,0,1,0})*gas.Xi[5]),gas.Xi);//in_T_internal;
  FormingGas.T = in_T_internal;
  if not use_in_T then
    in_T_internal = T "Temperature set by parameter";
  end if;

   if not use_in_mfFGPerCell then
    in_mfFGPerCell_internal = mfFGPerCell "Temperature set by parameter";
  end if;

  if molarInputComposition then
    Y = in_composition_internal[1:Medium.nXi];
    gas.Xi = Medium.moleToMassFractions(Y, Medium.MMX);
    Yfg = FormingGasComposition;
    FormingGas.Xi = Medium.moleToMassFractions(Yfg, Medium.MMX);
  else
    gas.Xi = in_composition_internal[1:Medium.nXi];
    Y = Medium.massToMoleFractions(gas.Xi, Medium.MMX);
    FormingGas.Xi = FormingGasComposition;
    Yfg = Medium.massToMoleFractions(FormingGas.Xi, Medium.MMX);
  end if;

  if not use_in_composition then
    in_composition_internal = compositionNom "Composition set by parameter";
  end if;

  out_y = (Y*w/gas.MM+Yfg*wFG/FormingGas.MM)/(w/gas.MM+wFG/FormingGas.MM+Modelica.Constants.small) "outlet molar fraction";

  if in_I_internal > 0 then
    out_I = if flange.Xi_outflow[1] < minH2OXiforIout then 0.001 else in_I_internal;
  else
    out_I = if flange.Xi_outflow[5] < minH2OXiforIout then 0.001 else in_I_internal;
  end if;

  FormingGas.p = gas.p "Pressure is equal and fixed";
  flange.p = FormingGas.p;
  flange.h_outflow = (w*gas.h+wFG*FormingGas.h)/(w+wFG+Modelica.Constants.small) "outgoing enthalpy is ratio of fuel and forming gas";
  flange.Xi_outflow = (w*gas.Xi+wFG*FormingGas.Xi)/(w+wFG+Modelica.Constants.small) "outgoing composition is ratio of fuel and forming gas";

  // Connect protected connectors to public conditional connectors
  connect(in_I, in_I_internal);
  connect(in_RC, in_RC_internal);
  connect(in_T, in_T_internal);
  connect(in_composition, in_composition_internal);
  connect(in_mfFGPerCell, in_mfFGPerCell_internal);

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
<li><i>19 Jan 2022</i>
    by <a href=\"mailto:Marius.Tomberg@dlr.de\">Marius Tomberg</a>:<br>
    Changed z of CH4 in EC to -4 as it consumes water</li>        
<li><i>26 Okt 2021</i>
    by <a href=\"mailto:Marius.Tomberg@dlr.de\">Marius Tomberg</a>:<br>
    Several pathes fixed</li>
<li><i>Sep 2021</i>
by <a href=\"mailto:Marius.Tomberg@dlr.de\">Marius Tomberg</a>:<br>
Adapted from ThermoPower

</ul>
</html>"),
    Diagram(coordinateSystem(extent={{-100,-100},{100,100}})),
    Icon(coordinateSystem(extent={{-100,-100},{100,100}})));
end SourceMassFlowEC;

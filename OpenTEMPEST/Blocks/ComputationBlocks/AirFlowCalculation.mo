within OpenTEMPEST.Blocks.ComputationBlocks;
block AirFlowCalculation
  "Air flow calculations for solid oxide cell systems"
  extends Modelica.Blocks.Icons.Block;

  replaceable package Medium = OpenTEMPEST.Medium.Air_Medium;

  parameter Integer nCells "Number of cells that are fed by this source";

  parameter Modelica.SIunits.ElectricCurrent I0=0 "Nominal electrical current ";
  parameter Real OU0=0.25 "Nominal oxygen utilization";
  parameter Real yO2Out0=0.40 "Nominal outlet oxygen molar compostion";
  parameter Real lambdaSF0=2 "Nominal air number (lambda, sunfire style) ";
  parameter Real zFuel0=1.8 "Nominal zFuel ";
  parameter Modelica.SIunits.MassFlowRate mfMin = 0.005 "Minimum mass flow";

  Modelica.Blocks.Interfaces.RealInput in_I(unit="A")
                                            if use_in_I "Current input"
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
  Modelica.Blocks.Interfaces.RealInput in_OU if use_in_OU "oxygen utilization input"
    annotation (Placement(transformation(extent={{-140,60},{-100,100}})));
  Modelica.Blocks.Interfaces.RealInput in_yO2Out if use_in_yO2Out "oulet oxygen molar compostion"
    annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
  Modelica.Blocks.Interfaces.RealInput in_lambdaSF if use_in_lambdaSF "air number (lambda, sunfire style)"
    annotation (Placement(transformation(extent={{-140,-100},{-100,-60}})));
  Modelica.Blocks.Interfaces.RealInput in_zFuel if use_in_zFuel "zFuel"
    annotation (Placement(transformation(extent={{-140,-60},{-100,-20}})));

  parameter Boolean use_in_I=false
    "Use connector input for the electrical current"
    annotation (Dialog(group="External inputs"), choices(checkBox=true));
  parameter Boolean use_in_OU=false
    "Use connector input for the oxygen utilization"
    annotation (Dialog(group="External inputs"), choices(checkBox=true));
  parameter Boolean use_in_yO2Out=false
    "Use connector input for yO2Out"
    annotation (Dialog(group="External inputs"), choices(checkBox=true));
  parameter Boolean use_in_lambdaSF=false
    "Use connector input for lambdaSF"
    annotation (Dialog(group="External inputs"), choices(checkBox=true));
  parameter Boolean use_in_zFuel=false
    "Use connector input for zFuel"
    annotation (Dialog(group="External inputs"), choices(checkBox=true));

  Modelica.Blocks.Interfaces.RealOutput out_mfAirOU(unit="kg/s")
    "air mass flow based on oxygen utilization (OU)"
    annotation (Placement(transformation(extent={{100,40},{120,60}})));
  Modelica.Blocks.Interfaces.RealOutput out_mfAirOutletO2(unit="kg/s")
    "air mass flow based on oxygen outlet molar fraction "
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));
  Modelica.Blocks.Interfaces.RealOutput out_mfAirLambdaSF(unit="kg/s")
    "air mass flow based on air number (λsf, sunfire style) "
    annotation (Placement(transformation(extent={{100,-60},{120,-40}})));

protected
  Medium.MoleFraction YAir[Medium.nX]=Medium.massToMoleFractions(Medium.reference_X,
      Medium.MMX);
  Modelica.Blocks.Interfaces.RealInput in_I_internal;
  Modelica.Blocks.Interfaces.RealInput in_OU_internal;
  Modelica.Blocks.Interfaces.RealInput in_yO2Out_internal;
  Modelica.Blocks.Interfaces.RealInput in_lambdaSF_internal;
  Modelica.Blocks.Interfaces.RealInput in_zFuel_internal;

equation

  out_mfAirOU = max(mfMin,sum(Medium.MMX*YAir)*in_I_internal*nCells/4/Modelica.Constants.F/in_OU_internal/
    YAir[1]); //
  out_mfAirOutletO2 = max(mfMin,sum(Medium.MMX*YAir)*in_I_internal*nCells/4/Modelica.Constants.F * (in_yO2Out_internal-1)/(in_yO2Out_internal-YAir[1]));
  out_mfAirLambdaSF = max(mfMin,sum(Medium.MMX*YAir)*abs(in_I_internal)*nCells/in_zFuel_internal/Modelica.Constants.F * in_lambdaSF_internal);

  if not use_in_I then
    in_I_internal = I0 "Current set by parameter";
  end if;
  if not use_in_OU then
    in_OU_internal = OU0 "OU set by parameter";
  end if;
  if not use_in_yO2Out then
    in_yO2Out_internal = yO2Out0 "yO2Out set by parameter";
  end if;
  if not use_in_lambdaSF then
    in_lambdaSF_internal = lambdaSF0 "lambdaSF set by parameter";
  end if;
  if not use_in_zFuel then
    in_zFuel_internal = zFuel0 "zFuel set by parameter";
  end if;

  // Connect protected connectors to public conditional connectors
  connect(in_I, in_I_internal);
  connect(in_OU, in_OU_internal);
  connect(in_yO2Out, in_yO2Out_internal);
  connect(in_lambdaSF, in_lambdaSF_internal);
  connect(in_zFuel, in_zFuel_internal);

  annotation(choicesAllMatching = true,
        Documentation(revisions="<html>
<ul>
<li><i>28 Okt 2021</i>
by <a href=\"mailto:Marius.Tomberg@dlr.de\">Marius Tomberg</a>:<br>
Created
</html>"),
    Diagram(coordinateSystem(extent={{-100,-100},{100,100}})),
    Icon(coordinateSystem(extent={{-100,-100},{100,100}})));
end AirFlowCalculation;

within OpenTEMPEST.Flow.BaseClasses;
partial model simpleGasHeaterBase

  import SI = Modelica.SIunits;

  replaceable package medium = Modelica.Media.Interfaces.PartialMedium annotation(choicesAllMatching = true);

  parameter SI.Length l=0.09 "lenght of flow channel";
  parameter  SI.Length h=0.01
    "height of flow channel between two parallel plates";
  parameter  SI.Length w=0.04
    "width of flow channel between two parallel plates";
  parameter  SI.Length d=0.02
    "diameter of flow channel in a pipe";

  parameter  SI.Temperature TStart=25 + 273.15;
  parameter  SI.Pressure pStart=101325;
  parameter  SI.MassFraction xStart[medium.nX]=medium.reference_X;
  parameter  SI.MassFlowRate maximalFlow=-1
    "maximal mass flow rate";
    parameter Boolean allowFlowReversal=system.allowFlowReversal
    "= true to allow flow reversal, false restricts to design direction"
    annotation(Evaluate=true);
  outer ThermoPower.System system "System wide properties";

   SI.CoefficientOfHeatTransfer alpha;
  Real Nu;
  parameter Boolean pipe=false
    "= true if flow through a pipe, false if between two parallel plates";
  parameter Integer nParallel=30 "number of parallel flow ducts";
  parameter Boolean useAlphaIn=false "Do not colculate heat transfer but use inlet coefficient";

  medium.DynamicViscosity eta "Dynamic viscosity";
  medium.SpecificHeatCapacity cp "Heat capacity at constant pressure";
   SI.PerUnit Re "Reynolds number";
   SI.PerUnit Pr "Prandtl numbers";
  medium.ThermalConductivity lambda "Thermal conductivity";

  ThermoPower.Gas.FlangeA inlet(redeclare package Medium = medium) annotation (
      Placement(transformation(extent={{-100,-10},{-80,10}}),
        iconTransformation(extent={{-100,-20},{-60,20}})));
  ThermoPower.Gas.FlangeB outlet(redeclare package Medium = medium) annotation (
     Placement(transformation(extent={{80,0},{100,20}}), iconTransformation(
          extent={{60,-20},{100,20}})));
  Modelica.Blocks.Interfaces.RealInput alphaIn if useAlphaIn annotation (
      Placement(transformation(
        extent={{20,-20},{-20,20}},
        rotation=-90,
        origin={-60,28}), iconTransformation(
        extent={{20,-20},{-20,20}},
        rotation=-90,
        origin={-40,40})));
  Modelica.Blocks.Sources.Constant const(k=-1)
    annotation (Placement(transformation(extent={{-80,-40},{-60,-20}})));
protected
  Modelica.Blocks.Interfaces.RealOutput alphaInHelper
    annotation (Placement(transformation(extent={{-38,-16},{-18,4}})));

protected
  parameter  SI.Area A=if pipe then nParallel*l*Modelica.Constants.pi
      *d else nParallel*2*l*(w + h) "Surface relevant for the heat transfer";
  parameter  SI.Area ACS=if pipe then Modelica.Constants.pi*d^2/4
       else w*h "cross-sectional area of the flow duct";
  parameter  SI.Length dHyd=if pipe then d else 2*h
    "characteristic length";

equation

  if not useAlphaIn then

    if pipe then
      Nu = 3.66 + 0.0677 * abs(Re * Pr*d/l).^1.33 ./ ( 1+0.1*Pr.*((Re+1e-9)*d/l).^0.83);
    else
      Nu = 7.541; //0.332 * (Re+1e-9)^(0.5) * Pr^(1/3);
    end if;

    if maximalFlow>0 then
      Re = abs( max(min(inlet.m_flow,2*maximalFlow),-2*maximalFlow) /nParallel*dHyd/(ACS*eta));
    else
      Re = abs(inlet.m_flow/nParallel*dHyd/(ACS*eta));
    end if;
    Pr = cp*eta/lambda;

    alpha = max(min(Nu * lambda / dHyd,1000),0.01); // high number of TStar erros

  else
    Pr = -1;
    Re = -1;
    Nu = -1;

    alpha = alphaInHelper;
  end if;

  connect(alphaIn, alphaInHelper) annotation (Line(points={{-60,28},{-48,28},{-48,
          -6},{-28,-6}}, color={0,0,127}));
  if not useAlphaIn then
    connect(const.y, alphaInHelper) annotation (Line(points={{-59,-30},{-46,-30},{
          -46,-6},{-28,-6}}, color={0,0,127}));
  end if;
  annotation (Icon(graphics={Rectangle(
          extent={{-60,20},{60,-20}},
          lineColor={28,108,200},
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid),
        Text(extent={{-100,-30},{100,-54}},  textString="%name")}));
end simpleGasHeaterBase;

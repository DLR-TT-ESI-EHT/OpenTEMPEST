within OpenTEMPEST.BOP.GasConditioning;
model ExpansionTank "Header with metal walls for water/steam flows"
  extends ThermoPower.Icons.Water.Header;
  replaceable package Medium = ThermoPower.Water.StandardWater constrainedby
    Modelica.Media.Interfaces.PartialMedium "Medium model"
    annotation(choicesAllMatching = true);
  Medium.ThermodynamicState fluidState "Thermodynamic state of the fluid";
  parameter Modelica.SIunits.Volume V "Inner volume";
  parameter Modelica.SIunits.Area S=0 "Internal surface";
  parameter Modelica.SIunits.Position H=0 "Elevation of outlet over inlet"
    annotation (Evaluate=true);
  parameter Modelica.SIunits.CoefficientOfHeatTransfer gamma=0
    "Heat Transfer Coefficient" annotation (Evaluate=true);
  parameter Modelica.SIunits.HeatCapacity Cm=0 "Metal Heat Capacity"
    annotation (Evaluate=true);
  parameter Boolean allowFlowReversal=system.allowFlowReversal
    "= true to allow flow reversal, false restricts to design direction"
    annotation(Evaluate=true);
  outer ThermoPower.System system "System wide properties";
  parameter ThermoPower.Choices.FluidPhase.FluidPhases FluidPhaseStart=ThermoPower.Choices.FluidPhase.FluidPhases.Liquid
    "Fluid phase (only for initialization!)"
    annotation (Dialog(tab="Initialisation"));
  parameter Medium.AbsolutePressure pstart "Pressure start value"
    annotation (Dialog(tab="Initialisation"));
  parameter Medium.SpecificEnthalpy hstart=if FluidPhaseStart == ThermoPower.Choices.FluidPhase.FluidPhases.Liquid
       then 1e5 else if FluidPhaseStart == ThermoPower.Choices.FluidPhase.FluidPhases.Steam
       then 3e6 else 1e6 "Specific enthalpy start value"
    annotation (Dialog(tab="Initialisation"));
  parameter Medium.Temperature Tmstart=300
    "Metal wall temperature start value"
    annotation (Dialog(tab="Initialisation"));
  parameter ThermoPower.Choices.Init.Options initOpt=system.initOpt
    "Initialisation option"
    annotation (Dialog(tab="Initialisation"));
  parameter Boolean noInitialPressure=false
    "Remove initial equation on pressure"
    annotation (Dialog(tab="Initialisation"),choices(checkBox=true));
  parameter Boolean noInitialEnthalpy=false
    "Remove initial equation on enthalpy"
    annotation (Dialog(tab="Initialisation"),choices(checkBox=true));

  ThermoPower.Water.FlangeA inlet(
    h_outflow(start=hstart),
    redeclare package Medium = Medium,
    m_flow(min=if allowFlowReversal then -Modelica.Constants.inf else 0))
    annotation (Placement(transformation(extent={{-122,-20},{-80,20}}, rotation=
           0)));
  ThermoPower.Water.FlangeB outlet(
    h_outflow(start=hstart),
    redeclare package Medium = Medium,
    m_flow(max=if allowFlowReversal then +Modelica.Constants.inf else 0))
    annotation (Placement(transformation(extent={{80,-20},{120,20}}, rotation=0)));
  Medium.AbsolutePressure p(start=pstart, stateSelect=if Medium.singleState then StateSelect.avoid
         else StateSelect.prefer) "Fluid pressure at the outlet";
  Medium.SpecificEnthalpy h(start=hstart, stateSelect=StateSelect.prefer)
    "Fluid specific enthalpy";
  Medium.SpecificEnthalpy hi "Inlet specific enthalpy";
  Medium.SpecificEnthalpy ho "Outlet specific enthalpy";
  Modelica.SIunits.Mass M "Fluid mass";
  Modelica.SIunits.Energy E "Fluid energy";
  Medium.Temperature T "Fluid temperature";
  Medium.Temperature Tm(start=Tmstart) "Wall temperature";
  Modelica.SIunits.Time Tr "Residence time";
  Real dM_dt;
  Real dE_dt;

  Modelica.SIunits.Volume VG "gas volume";
  Modelica.SIunits.Volume VL "liquid volume";
  parameter Real gasShareStart = 0.1 "Start share gas"
    annotation (Dialog(tab="Initialisation"));

  replaceable ThermoPower.Thermal.HT thermalPort "Internal surface of metal wall"
    annotation (Dialog(enable=false), Placement(transformation(extent={{-24,
            50},{24,64}}, rotation=0)));
equation

  VG = (pstart*gasShareStart*V/Modelica.Constants.R/Tmstart)*Modelica.Constants.R*T/p;

  V = VG+VL;

  // Set fluid properties
  fluidState = Medium.setState_ph(p, h);
  T = Medium.temperature(fluidState);

  M = VL*Medium.density(fluidState) "Fluid mass";
  E = M*h - p*VL "Fluid energy";
  dM_dt = VL*(Medium.density_derp_h(fluidState)*der(p) + Medium.density_derh_p(
    fluidState)*der(h));
  dE_dt = h*dM_dt + M*der(h) - VL*der(p);
  dM_dt = inlet.m_flow + outlet.m_flow "Fluid mass balance";
  dE_dt = inlet.m_flow*hi + outlet.m_flow*ho + gamma*S*(Tm - T) + thermalPort.Q_flow
    "Fluid energy balance";
  if Cm > 0 and gamma > 0 then
    Cm*der(Tm) = gamma*S*(T - Tm) "Energy balance of the built-in wall model";
  else
    Tm = T "Trivial equation for metal temperature";
  end if;

  // Boundary conditions
  hi = homotopy(if not allowFlowReversal then inStream(inlet.h_outflow) else
    actualStream(inlet.h_outflow), inStream(inlet.h_outflow));
  ho = homotopy(if not allowFlowReversal then h else actualStream(outlet.h_outflow),
    h);
  inlet.h_outflow = h;
  outlet.h_outflow = h;
  inlet.p = p + Medium.density(fluidState)*Modelica.Constants.g_n*H;
  outlet.p = p;
  thermalPort.T = T;

  Tr = noEvent(M/max(abs(inlet.m_flow), Modelica.Constants.eps))
    "Residence time";
initial equation
  // Initial conditions
  if initOpt == ThermoPower.Choices.Init.Options.noInit then
    // do nothing
  elseif initOpt == ThermoPower.Choices.Init.Options.fixedState then
    if not noInitialPressure then
      p = pstart;
    end if;
    if not noInitialEnthalpy then
      h = hstart;
    end if;
    if (Cm > 0 and gamma > 0) then
      Tm = Tmstart;
    end if;
  elseif initOpt == ThermoPower.Choices.Init.Options.steadyState then
    if not noInitialEnthalpy then
      der(h) = 0;
    end if;
    if (not Medium.singleState and not noInitialPressure) then
      der(p) = 0;
    end if;
    if (Cm > 0 and gamma > 0) then
      der(Tm) = 0;
    end if;
  elseif initOpt == ThermoPower.Choices.Init.Options.steadyStateNoP then
    if not noInitialEnthalpy then
      der(h) = 0;
    end if;
    if (Cm > 0 and gamma > 0) then
      der(Tm) = 0;
    end if;
  else
    assert(false, "Unsupported initialisation option");
  end if;
  annotation (
    Icon(graphics),
    Documentation(info="<HTML>
<p>This model describes a constant volume header with metal walls. The fluid can be water, steam, or a two-phase mixture.
<p>It is possible to take into account the heat storage and transfer in the metal wall in two ways:
<ul>
<li>
  Leave <tt>InternalSurface</tt> unconnected, and set the appropriate
  values for the total wall heat capacity <tt>Cm</tt>, surface <tt>S</tt>
  and heat transfer coefficient <tt>gamma</tt>. In this case, the metal
  wall temperature is considered as uniform, and the wall is thermally
  insulated from the outside.
</li>
<li>
  Set <tt>Cm = 0</tt>, and connect a suitable thermal model of the the
  wall to the <tt>InternalSurface</tt> connector instead. This can be
  useful in case a more detailed thermal model is needed, e.g. for
  thermal stress studies.
</li>
</ul>
<p>The model can represent an actual header when connected to the model of a bank of tubes (e.g., <tt>Flow1D</tt> with <tt>Nt>1</tt>).</p>
</HTML>",
        revisions="<html>
<ul>
<li><i>30 May 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Initialisation support added.</li>
<li><i>12 Apr 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       <tt>InternalSurface</tt> connector added.</li>
<li><i>16 Dec 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Standard medium definition added.</li>
<li><i>28 Jul 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Added head between inlet and outlet.</li>
<li><i>7 Jul 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Changed name from <tt>Collector</tt> to <tt>Header</tt>.</li>
<li><i>18 Jun 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted to Modelica.Media.</li>
<li><i>1 Oct 2003</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</html>
"), Diagram(graphics));
end ExpansionTank;

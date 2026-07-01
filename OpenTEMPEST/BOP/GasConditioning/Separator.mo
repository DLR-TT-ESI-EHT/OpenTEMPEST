within OpenTEMPEST.BOP.GasConditioning;
model Separator "Separator with metal walls for gas flows"
  extends ThermoPower.Icons.Gas.Mixer;
  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
    annotation(choicesAllMatching = true);
  Medium.BaseProperties gas(
    p(start=pstart, stateSelect=StateSelect.prefer),
    T(start=Tstart, stateSelect=StateSelect.prefer),
    Xi(start=Xstart[1:Medium.nXi], each stateSelect=StateSelect.prefer));
  Medium.BaseProperties gasM(
    p(start=pstart, stateSelect=StateSelect.prefer),
    T(start=Tstart, stateSelect=StateSelect.prefer),
    Xi(start=Xstart[1:Medium.nXi], each stateSelect=StateSelect.prefer));
  Medium.BaseProperties gasS(
    p(start=pstart, stateSelect=StateSelect.prefer),
    T(start=Tstart, stateSelect=StateSelect.prefer),
    Xi(start=Xstart[1:Medium.nXi], each stateSelect=StateSelect.prefer));
  parameter Modelica.SIunits.Volume V "Inner volume";
  parameter Modelica.SIunits.Area S=0 "Inner surface";
  parameter Modelica.SIunits.CoefficientOfHeatTransfer gamma=0
    "Heat Transfer Coefficient" annotation (Evaluate=true);
  parameter Modelica.SIunits.HeatCapacity Cm=0 "Metal heat capacity"
    annotation (Evaluate=true);
  parameter Boolean allowFlowReversal=system.allowFlowReversal
    "= true to allow flow reversal, false restricts to design direction"
    annotation(Evaluate=true);
  outer ThermoPower.System system "System wide properties";
  parameter Medium.AbsolutePressure pstart=1e5 "Pressure start value"
    annotation (Dialog(tab="Initialisation"));
  parameter Medium.Temperature Tstart=300 "Temperature start value"
    annotation (Dialog(tab="Initialisation"));
  parameter Medium.MassFraction Xstart[Medium.nX]=Medium.reference_X
    "Start gas composition" annotation (Dialog(tab="Initialisation"));
  parameter Medium.Temperature Tmstart=Tstart "Metal wall start temperature"
    annotation (Dialog(tab="Initialisation"));
  parameter ThermoPower.Choices.Init.Options initOpt=system.initOpt
    "Initialisation option"
    annotation (Dialog(tab="Initialisation"));
  parameter Boolean noInitialPressure=false
    "Remove initial equation on pressure"
    annotation (Dialog(tab="Initialisation"),choices(checkBox=true));
  parameter Boolean noInitialTemperature=false
    "Remove initial equation on temperature"
    annotation (Dialog(tab="Initialisation"),choices(checkBox=true));
  parameter Real seperationFactor[Medium.nXi] "Seperated fraction per species";

  Modelica.SIunits.Mass M "Gas total mass";
  Modelica.SIunits.InternalEnergy E "Gas total energy";
  Modelica.SIunits.Temperature Tm(start=Tmstart) "Wall temperature";
  Medium.SpecificEnthalpy hi "Inlet specific enthalpy";
  Medium.SpecificEnthalpy hoM "Outlet main specific enthalpy";
  Medium.SpecificEnthalpy hoS "Outlet seperated specific enthalpy";
  Medium.MassFraction Xi[Medium.nXi] "Inlet composition";
  Medium.MassFraction XoM[Medium.nXi] "Outlet main composition";
  Medium.MassFraction XoS[Medium.nXi] "Outlet seperated composition";
  Modelica.SIunits.Time Tr "Residence time";

  ThermoPower.Gas.FlangeA inlet(redeclare package Medium = Medium, m_flow(min=
          if allowFlowReversal then -Modelica.Constants.inf else 0))
    annotation (Placement(transformation(extent={{-120,-20},{-80,20}}, rotation=
           0)));
  ThermoPower.Gas.FlangeB outletMain(redeclare package Medium = Medium, m_flow(
        max=if allowFlowReversal then +Modelica.Constants.inf else 0))
    annotation (Placement(transformation(extent={{60,40},{100,80}}, rotation=0)));
  ThermoPower.Gas.FlangeB outletSeperated(redeclare package Medium = Medium,
      m_flow(min=if allowFlowReversal then -Modelica.Constants.inf else 0))
    annotation (Placement(transformation(extent={{60,-80},{100,-40}}, rotation=0)));

  replaceable ThermoPower.Thermal.HT thermalPort annotation (Placement(transformation(
          extent={{-38,60},{42,80}}, rotation=0)));
equation
  M = gas.d*V "Gas mass";
  E = M*gas.u "Gas internal energy";
  der(M) = inlet.m_flow + outletMain.m_flow + outletSeperated.m_flow "Mass balance";
  der(E) = inlet.m_flow*hi + outletMain.m_flow*hoM + outletSeperated.m_flow*hoS - gamma*S*(gas.T
     - Tm) + thermalPort.Q_flow "Energy balance";
  for j in 1:Medium.nXi loop
    M*der(gas.Xi[j]) = inlet.m_flow*(Xi[j] - gas.Xi[j]) + outletMain.m_flow*(XoM[j]
       - gas.Xi[j]) + outletSeperated.m_flow*(XoS[j] - gas.Xi[j])
      "Independent component mass balance";
  end for;
  if Cm > 0 and gamma > 0 then
    Cm*der(Tm) = gamma*S*(gas.T - Tm) "Metal wall energy balance";
  else
    Tm = gas.T;
  end if;

  // Seperation
  outletMain.m_flow * XoM  =  (ones(Medium.nXi)-seperationFactor) .* gas.Xi * (outletMain.m_flow + outletSeperated.m_flow);
  outletSeperated.m_flow * XoS =  seperationFactor .* gas.Xi * (outletMain.m_flow + outletSeperated.m_flow);
  hoM = gasM.h; //Medium.specificEnthalpy(gas.state);
  hoS = gasS.h;//Medium.specificEnthalpy(gas.state);

  // Boundary conditions
  hi = homotopy(if not allowFlowReversal then inStream(inlet.h_outflow) else
    actualStream(inlet.h_outflow), inStream(inlet.h_outflow));
  Xi = homotopy(if not allowFlowReversal then inStream(inlet.Xi_outflow) else
    actualStream(inlet.Xi_outflow), inStream(inlet.Xi_outflow));
  hoM = homotopy(if not allowFlowReversal then gas.h else actualStream(outletMain.h_outflow),
    gas.h);
  XoM = homotopy(if not allowFlowReversal then gas.Xi else actualStream(outletMain.Xi_outflow),
    gas.Xi);
  hoS = homotopy(if not allowFlowReversal then gas.h else actualStream(outletSeperated.h_outflow),
    gas.h);
  XoS = homotopy(if not allowFlowReversal then gas.Xi else actualStream(outletSeperated.Xi_outflow),
    gas.Xi);
  inlet.p = gas.p;
  inlet.h_outflow = gas.h;
  inlet.Xi_outflow = gas.Xi;
  //out.p = gas.p;
  //out.h_outflow = gas.h;
  //out.Xi_outflow = gas.Xi;
  outletMain.p = gas.p;
  outletSeperated.p = gas.p;
  thermalPort.T = gas.T;

  gasM.p = gas.p;
  gasS.p = gas.p;
  gasM.T = gas.T;
  gasS.T = gas.T;
  gasS.Xi = XoS;
  gasM.Xi = XoM;

  Tr = noEvent(M/max(abs(-outletMain.m_flow-outletSeperated.m_flow), Modelica.Constants.eps))
    "Residence time";
initial equation
  // Initial conditions
  if initOpt == ThermoPower.Choices.Init.Options.noInit then
    // do nothing
  elseif initOpt == ThermoPower.Choices.Init.Options.fixedState then
    if not noInitialPressure then
      gas.p = pstart;
    end if;
    if not noInitialTemperature then
      gas.T = Tstart;
    end if;
    gas.Xi = Xstart[1:Medium.nXi];
    if (Cm > 0 and gamma > 0) then
      Tm  = Tmstart;
    end if;
  elseif initOpt == ThermoPower.Choices.Init.Options.steadyState then
    if not noInitialPressure then
      der(gas.p) = 0;
    end if;
    if not noInitialTemperature then
      der(gas.T) = 0;
    end if;
    der(gas.Xi) = zeros(Medium.nXi);
    if (Cm > 0 and gamma > 0) then
      der(Tm) = 0;
    end if;
  elseif initOpt == ThermoPower.Choices.Init.Options.steadyStateNoP then
    if not noInitialTemperature then
      der(gas.T) = 0;
    end if;
    der(gas.Xi) = zeros(Medium.nXi);
    if (Cm > 0 and gamma > 0) then
      der(Tm) = 0;
    end if;
  else
    assert(false, "Unsupported initialisation option");
  end if;

  annotation (
    Documentation(info="<html>
<p>This model describes a constant volume mixer with metal walls. The metal wall temperature and the heat transfer coefficient between the wall and the fluid are uniform. The wall is thermally insulated from the outside.</p>
<p><b>Modelling options</b></p>
<p>The actual gas used in the component is determined by the replaceable <tt>Medium</tt> package. In the case of multiple component, variable composition gases, the start composition is given by <tt>Xstart</tt>, whose default value is <tt>Medium.reference_X</tt>.
</html>",
        revisions="<html>
<ul>
<li><i>30 May 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Initialisation support added.</li>
<li><i>19 Nov 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted to Modelica.Media.</li>
<li><i>5 Mar 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</html>"));
end Separator;

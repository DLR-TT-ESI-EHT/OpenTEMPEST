within OpenTEMPEST.SOC.Cell.Cell1D;
partial model Channel1DBase

  import SI = Modelica.SIunits;
  replaceable package Medium = Modelica.Media.IdealGases.Common.MixtureGasNasa
    annotation(choicesAllMatching = true);
  parameter Integer N(min=3) = 5;
  constant Integer nSpecies = Medium.nXi;

  replaceable function fluxInterp =
      Flow.FluxInterpolators.UDSinterp              constrainedby
    Flow.FluxInterpolators.DifferencingSchemeInterpBase                                  annotation(choicesAllMatching = true);

  // Initial Values
  parameter SI.Temperature TStartIn = 773.15 annotation (Dialog(tab="Initialisation"));
  parameter SI.Temperature TStartOut = 773.15 annotation (Dialog(tab="Initialisation"));
  parameter SI.Temperature TStart[N] = linspace(TStartIn, TStartOut, N) annotation (Dialog(tab="Initialisation"));
  parameter SI.AbsolutePressure pStartIn = 101325 annotation (Dialog(tab="Initialisation"));
  parameter SI.AbsolutePressure pStartOut = pStartIn annotation (Dialog(tab="Initialisation"));
  parameter SI.AbsolutePressure pStart[N+1] = linspace(pStartIn,pStartOut,N+1) annotation (Dialog(tab="Initialisation"));
  parameter SI.MassFraction xStart[nSpecies] = Medium.reference_X annotation (Dialog(tab="Initialisation"));

  // Dimensions
  parameter SI.Length lY "width of channel" annotation (Dialog(tab="Dimensions"));
  parameter SI.Length lZ "Height of channel"
                                            annotation (Dialog(tab="Dimensions"));
  parameter SI.Length lX = 1 "Length of channel /m"
                                                   annotation (Dialog(tab="Dimensions"));
  parameter Real por "channel porosity"
                                       annotation (Dialog(tab="Dimensions"));

  parameter Real Nu_PEN "Nusselt number on PEN side";
  parameter Real Nu_IC "Nusselt Number on IC side";

  // Control Volume sizing
  SI.Volume dV = dx*lY*lZ;
  SI.Length dx = lX/N;

  // Gas Object
  Medium.BaseProperties Gas[N](
    p(start=pStart[2:N+1], each stateSelect=StateSelect.prefer),
    T(start=TStart[1:N], each stateSelect=StateSelect.prefer),
    Xi(each stateSelect=StateSelect.prefer))
    "Gas volume properties";

  // Centre cell values
  SI.MassFlowRate mf[N] "Mass flow rate in and leaving CV";
  SI.EnergyFlowRate QgasExt[N] "Gas phase heat flows in CV";
  Real Ycell[N, nSpecies];

  // Vertex Values
  SI.MassFlowRate mfv[N+1] "Mass flow rate CV Vertices/Nodes";
  SI.SpecificEnthalpy hv[N+1] "Enthalpy at Vertices";
  SI.MassFraction xiv[N+1, nSpecies] "Mass fractions at Vertices";

  // Kinetics
  SI.MassFlowRate R[N, nSpecies]
    "net Rate of production and consumption of products and reactants in the reactor - from thermochemical AND electrochemical reactions";

  SI.Energy Emg[N] "Energy density in control volume";
  SI.MassFlowRate massTransfer[N] "ion transfer rate";

  SI.Area crossSectionalArea = lY*lZ;
  SI.Length wettedPerimeter = 2*(lY+lZ);
  SI.Length Dh = 4*crossSectionalArea/wettedPerimeter;

  Real Re[N] = (mf[:]/crossSectionalArea * Dh)./(eta);
  Real Pr[N] = Medium.prandtlNumber(Gas.state);
  SI.DynamicViscosity eta[N];

  SI.EnergyFlowRate deltaH = infl.m_flow*infl.h_outflow + outfl.m_flow*outfl.h_outflow;
  SI.EnergyFlowRate EB = deltaH + sum(QgasExt);

  ThermoPower.Gas.FlangeA infl(redeclare package Medium = Medium)  annotation (Placement(transformation(extent={{-100,-8},{-80,12}}),
        iconTransformation(extent={{-108,-16},{-80,12}})));
  ThermoPower.Gas.FlangeB outfl(redeclare package Medium = Medium)  annotation (Placement(transformation(extent={{80,-10},{100,10}}),
        iconTransformation(extent={{80,-14},{108,14}})));
  ThermoPower.Thermal.DHTVolumes Q_IC(N=N) "HT Between Channel and IC"  annotation (Placement(transformation(extent={{-62,34},{-40,52}}),
        iconTransformation(extent={{-26,34},{26,46}})));
  ThermoPower.Thermal.DHTVolumes Q_PEN(N=N) "HT Between Channel and PEN"  annotation (Placement(transformation(extent={{-60,-32},{-38,-14}}),
        iconTransformation(extent={{-26,-40},{24,-28}})));
  Electrochem.Interfaces.VariablesStream PEN_in[N](each nspecies=nSpecies)
    annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=-90,
        origin={-14,-22}), iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=-90,
        origin={-52,-32})));

initial equation
  Gas[:].T = TStart[:];
  Gas[:].X[:] = fill(xStart[:], N);
  //Gas[1:N-1].p = pStart[2:N];
equation

  // Total Mass Balance
   dV.*por.*der(Gas[:].d) = (mfv[1:N] .- mfv[2:N+1]) .+ massTransfer[:];

   // Species Mass Balance
   for i in 1:nSpecies loop
      dV.*por*der(Gas[:].d.*Gas[:].Xi[i]) = mfv[1:N].*xiv[1:N, i] .- mfv[2:N+1].*xiv[2:N+1, i].+ R[:, i];
   end for;

  for i in 1:N loop
    Ycell[i,:] = Gas[i].Xi[:]./Medium.MMX[:]/sum(Gas[i].Xi[:]./Medium.MMX[:]);
    eta[i] = Medium.dynamicViscosity(Gas[i].state);
  end for;

  // Energy Balance
  der(Emg[:]) = mfv[1:N].*hv[1:N] .- mfv[2:N+1].*hv[2:N+1] .+ QgasExt[:];

  Q_PEN.Q[:] = Nu_PEN*Medium.thermalConductivity(Gas[:].state)./(2*lZ)*dx*lY.*(Q_PEN.T[:] .- Gas[:].T);
  Q_IC.Q[:]  = Nu_IC *Medium.thermalConductivity(Gas[:].state)./(2*lZ)*dx*lY.*(Q_IC.T[:]  .- Gas[:].T); // First simplifying assumption (for the time being): negligible convective HT from the sides of the IC

  // Boundary conditions
  outfl.m_flow = -mfv[N+1];
  outfl.Xi_outflow[:] = xiv[N+1, :];
  outfl.h_outflow = hv[N+1];

  // Interface interpolation
    // At x=0
  mfv[1] = infl.m_flow;
  hv[1]  = infl.h_outflow;
  hv[1]  = inStream(infl.h_outflow);
  xiv[1, :] = infl.Xi_outflow[:];
  xiv[1, :] = inStream(infl.Xi_outflow[:]);

    // Interior and x=L
  mfv[2:N+1] = fluxInterp(N, mf, mfv[1]);
  hv[2:N+1] = fluxInterp(N, Gas.h, hv[1]);
  for i in 1:nSpecies loop
    xiv[2:N+1, i] = fluxInterp(N, Gas.Xi[i], xiv[1, i]);
  end for;

  // To PEN
  for i in 1:N loop
    PEN_in[i].P = Gas[i].p;
    PEN_in[i].Y[:] = Ycell[i, :];
  end for;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
          extent={{-84,34},{82,-28}},
          lineColor={28,108,200},
          fillColor={28,108,200},
          fillPattern=FillPattern.HorizontalCylinder)}),         Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<h2>Channel1DBase</h2>

<p>
Partial finite-volume base model for one-dimensional gas flow channels. The channel is discretized into N control volumes 
(CVs) along its length.
</p>

<h3>Governing Equations</h3>
<ul>
<li>Total mass conservation (finite-volume form)</li>
<li>Species mass conservation for each gas component</li>
<li>Gas-phase energy conservation (internal energy formulation)</li>
<li>Convective fluxes computed via replaceable interpolation scheme</li>
<li>Heat exchange with PEN and interconnect via distributed heat ports</li>
</ul>
</html>"));
end Channel1DBase;

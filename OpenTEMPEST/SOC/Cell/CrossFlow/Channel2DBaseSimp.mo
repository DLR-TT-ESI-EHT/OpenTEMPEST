within OpenTEMPEST.SOC.Cell.CrossFlow;
partial model Channel2DBaseSimp
  "Base model for SOC cell channels 1D discretised with 2D interfaces - need to give pressure drop at top level"

  import SI = Modelica.SIunits;
  replaceable package Medium = Modelica.Media.IdealGases.Common.MixtureGasNasa
    annotation(choicesAllMatching = true);
  parameter Integer nX(min=3) = 5 "Number of control volumes in first direction";
  parameter Integer nY(min=3) = 5 "Number of control volumes in second direction for 2D interface";
  constant Integer nSpecies = Medium.nXi "Number of different species in the channel";
  parameter Real alfa(min=0, max=1) = 0.01 "Weight factor for local vs. averaged 2D interface temperature in convective heat transfer (0=fully averaged, 1=fully local)";

  // Selections
  parameter Boolean heatTransferCorrelationFormDuct = true "True for Nusselt correlation duct geometry with characteristic length=2*lZ (default), false for plate geometry with characteristic length=lX";
  replaceable function fluxInterp =
      Flow.FluxInterpolators.UDSinterp              constrainedby
    Flow.FluxInterpolators.DifferencingSchemeInterpBase                                  annotation(choicesAllMatching = true);

  // Initial Values
  parameter SI.Temperature TStart=773.15 annotation (Dialog(tab="Initialisation"));
  parameter SI.AbsolutePressure pStart=101325 annotation (Dialog(tab="Initialisation"));
  parameter SI.MassFraction xStart[nSpecies] = Medium.reference_X annotation (Dialog(tab="Initialisation"));

  // Dimensions
  parameter SI.Length lX = 1 "Length of channel" annotation (Dialog(tab="Dimensions"));
  parameter SI.Length lY = 1 "Width of channel" annotation (Dialog(tab="Dimensions"));
  parameter SI.Length lZ = 1 "Height of channel" annotation (Dialog(tab="Dimensions"));
  parameter Real por = 0.7 "Channel porosity"  annotation (Dialog(tab="Dimensions"));

  parameter Real Nu_PEN "Nusselt number for convective heat transfer at PEN interface";
  parameter Real Nu_IC "Nusselt number for convective heat transfer at IC interface";

  // Ribs thermophysical parameters - CFY (94.9% Cr, 5% Fe, 0.1% Y)
  parameter Modelica.SIunits.ThermalConductivity kRibs=40 "Thermal conductivity of CFY ribs in W/mK (35-45 W/mK for 20-900 °C)"  annotation(Dialog(tab="Ribs"));
  parameter Modelica.SIunits.SpecificHeatCapacity cpRibs = 451.8 "Specific heat capacity of ribs material"  annotation(Dialog(tab="Ribs"));
  parameter Modelica.SIunits.Density rhoRibs = 7233 "Calculated CFY ribs density in kg/m3"  annotation(Dialog(tab="Ribs"));

  // Control Volume sizing
  SI.Volume dV = dx*lY*lZ "Control volume size";
  SI.Length dx = lX/nX "Axial length of one control volume";

  // Gas Object
  SI.Temperature T[nX](each start=TStart) "Average Temperature in CV";

  Medium.BaseProperties Gas[nX](
    p(each start=pStart, each stateSelect=StateSelect.prefer),
    T(each start=TStart, each stateSelect=StateSelect.prefer),
    Xi(start=fill(xStart, nX), each stateSelect=StateSelect.prefer))
    "Gas volume properties";
  SI.ThermalConductivity[nX] lambdaGas = Medium.thermalConductivity(Gas[:].state) "Thermal conductivity of the gas in the CV";

  // Centre cell values
  SI.MassFlowRate mf[nX] "Cell-centered mass flow rate (used for flux interpolation)";
  SI.EnergyFlowRate QgasExt[nX] "Net external heat flow rate added to gas control volume";
  Modelica.Media.Interfaces.Types.MoleFraction Ycell[nX, nSpecies] "Mole fractions in control volumes (derived from mass fractions)";

  // Vertex Values
  SI.MassFlowRate mfv[nX+1] "Mass flow rate CV Vertices/Nodes";
  SI.SpecificEnthalpy hv[nX+1] "Enthalpy at Vertices";
  SI.MassFraction xiv[nX+1, nSpecies] "Mass fractions at Vertices";

  // Kinetics
  SI.MassFlowRate R[nX, nSpecies]
    "Net rate of production and consumption of products and reactants in the reactor - from thermochemical AND electrochemical reactions";

  SI.Energy Emg[nX] "Total gas internal energy in control volume";
  SI.MassFlowRate massTransfer[nX] "Mass source term due to electrochemical ion transport";

  // Inlet and outlet flanges
  ThermoPower.Gas.FlangeA infl(redeclare package Medium = Medium)  annotation (Placement(transformation(extent={{-106,0},
            {-86,20}}),
        iconTransformation(extent={{-114,-8},{-86,20}})));
  ThermoPower.Gas.FlangeB outfl(redeclare package Medium = Medium)  annotation (Placement(transformation(extent={{80,-10},{100,10}}),
        iconTransformation(extent={{86,-8},{114,20}})));

  // Heat ports PEN and IC
  Heat.DHTVolumes2D Q_IC(i=nX, j=nY) "2D HT Between Channel and IC" annotation (
     Placement(transformation(extent={{-28,34},{28,46}}), iconTransformation(
          extent={{-28,34},{28,46}})));
  Heat.DHTVolumes2D Q_PEN(i=nX, j=nY) "HT between Channel and PEN" annotation (
      Placement(transformation(extent={{-28,-40},{28,-28}}), iconTransformation(
          extent={{-28,-40},{28,-28}})));

  // Variables stream objects
  Electrochem.Interfaces.VariablesStream PEN_in[nX,nY](each nspecies=nSpecies)
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-50,-30}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-50,-30})));

  SI.Length Dhth = if heatTransferCorrelationFormDuct then (2*lZ) else (dx);

initial equation
  Gas.T = fill(TStart, nX);
  Gas.X = fill(xStart, nX);
  Gas[1:nX-1].p = fill(pStart, nX-1);
equation

  // Total Mass Balance 2D
  for i in 1:nX loop
    dV*por*der(Gas[i].d) = (mfv[i] .- mfv[i+1]) .+ massTransfer[i];
    Ycell[i,:] = Gas[i].Xi[:]./Medium.MMX[:]./sum(Gas[i].Xi[:]./Medium.MMX[:]);
  end for;

  // Species Mass Balance
  for i in 1:nSpecies loop
    dV.*por.*der(Gas[:].d.*Gas[:].Xi[i]) = mfv[1:nX].*xiv[1:nX, i] .- mfv[2:nX+1].*xiv[2:nX+1, i] .+ R[:, i];
  end for;

  // Energy Balance
  der(Emg[:]) = mfv[1:nX].*hv[1:nX] .- mfv[2:nX+1].*hv[2:nX+1] .+ QgasExt[:];


  for i in 1:nX loop
    for j in 1:nY loop
       Q_PEN.Q[i,j] = Nu_PEN*lambdaGas[i]/Dhth*por*dx*lY*(alfa*(Q_PEN.T[i,j] - Gas[i].T) + (1-alfa)*((sum(Q_PEN.T[i,:])/nY) - Gas[i].T));
       Q_IC.Q[i,j] = Nu_IC*lambdaGas[i]/Dhth*por*dx*lY*(alfa*(Q_IC.T[i,j] - Gas[i].T) + (1-alfa)*((sum(Q_IC.T[i,:])/nY) - Gas[i].T));
    end for;
  end for;

  // Boundary conditions
  outfl.m_flow = -mfv[nX+1];
  outfl.Xi_outflow[:] = xiv[nX+1, :];
  outfl.h_outflow = hv[nX+1];

  // Interface interpolation
  // At x=0
  mfv[1] = infl.m_flow;
  hv[1]  = infl.h_outflow;
  hv[1]  = inStream(infl.h_outflow);
  xiv[1, :] = infl.Xi_outflow[:];
  xiv[1, :] = inStream(infl.Xi_outflow[:]);

  // Interior and x=L
  mfv[2:nX+1] = fluxInterp(nX, mf[1:nX], mfv[1]);
  hv[2:nX+1] = fluxInterp(nX, Gas[1:nX].h, hv[1]);
  for i in 1:nSpecies loop
    xiv[2:nX+1, i] = fluxInterp(nX, Gas[1:nX].Xi[i], xiv[1, i]);
  end for;

  // To PEN 2D interface
  for i in 1:nX loop
    for j in 1:nY loop
      PEN_in[i,j].P = Gas[i].p; // Imposing the same pressure in the other CVs
      for k in 1:nSpecies loop
        PEN_in[i,j].Y[k] = Ycell[i,k];  // Imposing the same mole fractions in the other CVs
      end for;
    end for;
  end for;

  annotation (Documentation(revisions="<html>
</html>", info="<html>
<h2>Channel2DBaseSimp</h2>

<p>
Simplified channel base model for Solid Oxide Cell (SOC) applications.
The model resolves the flow and thermodynamic state in the streamwise
direction (X) using a 1D finite-volume discretisation, while maintaining
a 2D interface in the transverse direction (Y) for thermal and
electrochemical coupling.
</p>
<h3>Discretisation Concept</h3>
<ul>
<li>1D finite-volume formulation in X-direction (<code>nX</code> control volumes)</li>
<li>2D heat-transfer and electrochemical interface in Y-direction (<code>nY</code>)</li>
<li>Uniform gas properties assumed across Y within each X control volume</li>
<li>Convective fluxes computed via replaceable interpolation scheme</li>
</ul>
</html>"), Icon(graphics={
          Rectangle(
          extent={{-100,34},{100,-28}},
          lineColor={28,108,200},
          fillColor={28,108,200},
          fillPattern=FillPattern.HorizontalCylinder)}));
end Channel2DBaseSimp;

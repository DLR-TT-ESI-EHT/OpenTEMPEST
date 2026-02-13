within OpenTEMPEST.SOC.Cell.CrossFlow;
partial model Channel2DBase
  "Base model for 2D SOC cell channels - need to give energy balance and pressure drop at top level"

  replaceable package Medium = Modelica.Media.IdealGases.Common.MixtureGasNasa
    annotation(choicesAllMatching = true);
  parameter Integer nX(min=3) = 5;
  parameter Integer nY(min=1) = 5;
  constant Integer nSpecies = Medium.nXi;

  replaceable function fluxInterp =
      Flow.FluxInterpolators.UDSinterp              constrainedby
    Flow.FluxInterpolators.DifferencingSchemeInterpBase                                  annotation(choicesAllMatching = true);

  parameter Boolean heatTransferCorrelationFormDuct=true "true for Nusselt correlation duct geometry with characteristic length=2*lZ (default), false for plate geometry with characteristic length=lX";

  // Initial Values
  parameter Modelica.SIunits.Temperature TStart=773.15 annotation (Dialog(tab="Initialisation"));
  parameter Modelica.SIunits.AbsolutePressure pStart=101325 annotation (Dialog(tab="Initialisation"));
  parameter Modelica.SIunits.MassFraction xStart[nSpecies]=Medium.reference_X annotation (Dialog(tab="Initialisation"));

  // Dimensions
  parameter Modelica.SIunits.Length lX=1 "Length of channel /m" annotation (Dialog(tab="Dimensions"));
  parameter Modelica.SIunits.Length lY=1 "width of channel" annotation (Dialog(tab="Dimensions"));
  parameter Modelica.SIunits.Length lZ=1 "Height of channel" annotation (Dialog(tab="Dimensions"));
  parameter Real por = 0.7 "channel porosity" annotation (Dialog(tab="Dimensions"));

  parameter Real Nu_PEN "Nusselt number on PEN side";
  parameter Real Nu_IC "Nusselt Number on IC side";

  // Control Volume sizing
  parameter Modelica.SIunits.Volume dV=dx*dy*lZ annotation (Dialog(tab="Dimensions"));
  parameter Modelica.SIunits.Length dx=lX/nX annotation (Dialog(tab="Dimensions"));
  parameter Modelica.SIunits.Length dy=lY/nY annotation (Dialog(tab="Dimensions"));

  // Gas Object
  Modelica.SIunits.Temperature T[nX,nY](each start=TStart) "Average Temperature in CV";

  Medium.BaseProperties Gas[nX, nY](
    p(each start=pStart, each stateSelect=StateSelect.prefer),
    T(each start=TStart, each stateSelect=StateSelect.prefer),
    Xi(start=fill(xStart, nX, nY), each stateSelect=StateSelect.prefer))
    "Gas volume properties";
  Modelica.SIunits.ThermalConductivity[nX,nY] lambdaGas = Medium.thermalConductivity(Gas[:,:].state) "Thermal conductivity of the gas in the CV";

  // Centre cell values
  Modelica.SIunits.MassFlowRate mf[nX,nY] "Mass flow rate in and leaving CV";
  Modelica.SIunits.EnergyFlowRate QgasExt[nX,nY] "Gas phase heat flows in CV";
  Real Ycell[nX, nY, nSpecies];
  Modelica.SIunits.DynamicViscosity eta[nX,nY]=Medium.dynamicViscosity(Gas.state);

  // Vertex Values
  Modelica.SIunits.MassFlowRate mfv[nX + 1,nY] "Mass flow rate CV Vertices/Nodes";
  Modelica.SIunits.SpecificEnthalpy hv[nX + 1,nY] "Enthalpy at Vertices";
  Modelica.SIunits.MassFraction xiv[nX + 1,nY,nSpecies] "Mass fractions at Vertices";

  // Kinetics
  Modelica.SIunits.MassFlowRate R[nX,nY,nSpecies] "net Rate of production and consumption of products and reactants in the reactor - from thermochemical AND electrochemical reactions";

  ThermoPower.Gas.FlangeA infl[nY](redeclare package Medium = Medium)  annotation (Placement(transformation(extent={{-100,-8},{-80,12}}),
        iconTransformation(extent={{-108,-16},{-80,12}})));
  ThermoPower.Gas.FlangeB outfl[nY](redeclare package Medium = Medium)  annotation (Placement(transformation(extent={{80,-10},{100,10}}),
        iconTransformation(extent={{80,-14},{108,14}})));
  Heat.DHTVolumes2D Q_IC(i=nX, j=nY) "HT Between Channel and IC" annotation (
      Placement(transformation(extent={{-62,34},{-40,52}}), iconTransformation(
          extent={{-26,34},{26,46}})));
  Heat.DHTVolumes2D Q_PEN(i=nX, j=nY) "HT Between Channel and PEN" annotation (
      Placement(transformation(extent={{-60,-32},{-38,-14}}),
        iconTransformation(extent={{-26,-40},{24,-28}})));
  Electrochem.Interfaces.VariablesStream PEN_in[nX,nY](each nspecies=nSpecies)
    annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=-90,
        origin={-14,-22}), iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=-90,
        origin={-52,-32})));
  // Differentiating between air and fuel side
  parameter Integer channelFac "set to +1 for fuel channel, -1 for air channel";

  Modelica.SIunits.Length Dhth = if heatTransferCorrelationFormDuct then (2*lZ) else (dx);

initial equation
  Gas.T = fill(TStart, nX, nY);
  Gas.X = fill(xStart, nX, nY);
  for i in 1:nX-1 loop
    for j in 1:nY loop
      Gas[i,j].p = pStart;
    end for;
  end for;

equation

  // Total Mass Balance
  dV.*por.*der(Gas[:,:].d) = (mfv[1:nX, :] .- mfv[2:nX+1, :]) .+ channelFac*PEN_in.I./(4*Modelica.Constants.F).*Modelica.Media.IdealGases.SingleGases.O2.data.MM;

  // Species Mass Balance
  for i in 1:nSpecies loop
    dV.*por*der(Gas[:,:].d.*Gas[:,:].Xi[i]) = mfv[1:nX,:].*xiv[1:nX,:, i] .- mfv[2:nX+1,:].*xiv[2:nX+1,:, i] .+ R[:,:, i];
  end for;

  for i in 1:nX loop
    for k in 1:nY loop
      Ycell[i,k,:] = Gas[i,k].Xi[:]./Medium.MMX[:]./sum(Gas[i,k].Xi[:]./Medium.MMX[:]);
    end for;
  end for;

  // Energy Balance - Thermal Equilibrium between Gas and Solid phase
  dV.*der(por.*Gas[1,:].d.*Gas[1,:].u)                 = mfv[1,:].*hv[1,:]           .- mfv[2,:].*hv[2,:]       .+ QgasExt[1,:];
  dV.*der(por.*Gas[2:nX-1,:].d.*Gas[2:nX-1,:].u)  = mfv[2:nX-1,:].*hv[2:nX-1,:] .- mfv[3:nX,:].*hv[3:nX,:] .+ QgasExt[2:nX-1,:];
  dV.*der(por.*Gas[nX,:].d.*Gas[nX,:].u)              = mfv[nX,:].*hv[nX,:]         .- mfv[nX+1,:].*hv[nX+1,:] .+ QgasExt[nX,:];

  Gas.T = T;

  Q_PEN.Q[:,:] = Nu_PEN*lambdaGas./Dhth*por*dx*dy.*(Q_PEN.T .- Gas.T);
  Q_IC.Q[:,:]  = Nu_IC *lambdaGas./Dhth*por*dx*dy.*(Q_IC.T .-  Gas.T); // First simplifying assumption (for the time being): neglegible convective HT from the sides of the IC

  // Boundary conditions
  outfl.m_flow = -mfv[nX+1,:];
  outfl.Xi_outflow[:] = xiv[nX+1,:, :];
  outfl.h_outflow = hv[nX+1,:];

  // Interface interpolation
    // At x=0
  mfv[1,:] = infl.m_flow;
  hv[1,:]  = infl.h_outflow;
  hv[1,:]  = inStream(infl.h_outflow);
  xiv[1,:, :] = infl.Xi_outflow[:];
  xiv[1,:, :] = inStream(infl.Xi_outflow[:]);

    // Interior and x=L
  for j in 1:nY loop
    mfv[2:nX+1, j] = fluxInterp(nX, mf[:,j], mfv[1, j]);
    hv[2:nX+1, j] = fluxInterp(nX, Gas[:,j].h, hv[1, j]);
    for i in 1:nSpecies loop
      xiv[2:nX+1, j, i] = fluxInterp(nX, Gas[:,j].Xi[i], xiv[1, j, i]);
    end for;
  end for;

  // To PEN
  PEN_in[:,:].P = Gas.p;
  PEN_in[:,:].Y[:] = Ycell[:,:, :];

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
          extent={{-84,34},{82,-28}},
          lineColor={28,108,200},
          fillColor={28,108,200},
          fillPattern=FillPattern.HorizontalCylinder)}),         Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end Channel2DBase;

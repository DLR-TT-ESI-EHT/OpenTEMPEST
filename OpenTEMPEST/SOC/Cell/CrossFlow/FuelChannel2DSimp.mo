within OpenTEMPEST.SOC.Cell.CrossFlow;
model FuelChannel2DSimp "Fuel channel model 1D discretised with 2D interfaces"

  import SI = Modelica.SIunits;

  extends Channel2DBaseSimp(
    heatTransferCorrelationFormDuct=true,
    redeclare package Medium = OpenTEMPEST.Medium.Fuel_CH4);

  // Kinetic parameters
  constant Real A_WGSf = 0.0171  annotation(Dialog(tab="Kinetics"));
  constant Real Ea_WGSf = 103191  annotation(Dialog(tab="Kinetics"));
  constant Real A_MSRf = 2395  annotation(Dialog(tab="Kinetics"));
  constant Real Ea_MSRf = 231266*0.875  annotation(Dialog(tab="Kinetics"));

  // Pressure Drop
  SI.PressureDifference dp[nX] = Blocks.Functions.pressureDropDarcy(mfv[1:nX], eta, Gas.d, dx, (lZ - hAbsetzen)*lY, 3e-9) "pressure loss";
  SI.DynamicViscosity eta[nX] = Medium.dynamicViscosity(Gas.state);
  SI.Length hAbsetzen = 0;

  // Kinetics
  Real kfMSR[nX], kfWGS[nX];
  SI.MolarFlowRate rMSR[nX], rWGS[nX], rHEL[nX], rCEL[nX] "Reaction rates (mol/s) of MSR, WGS and steam and CO2 electrolysis";
  SI.MassFlowRate R[nX, nSpecies] "Net rate of production and consumption of products and reactants in the reactor";
  Real Kmsr[nX], Kwgs[nX] "Reaction equillibrium constant of reverse reforming and WGS reaction respectively";
  SI.MolarEnergy DeltaG_wgs[nX], DeltaG_msr[nX];

  SI.EnergyFlowRate Qcond[nX];
  SI.EnergyFlowRate q_electrochem[nX] =  rHEL.*(H_h2 .- H_h2o) .+ rCEL.*(H_co .- H_co2);

protected
    constant Real a[4, nSpecies]={{3,-1,0,1,-1,0}, {1,0,1,-1,-1,0}, {-1,0,0,0,1,0}, {0,0,1,-1,0,0}} "Stoichiometric coefficients for 1) MSR, 2) WGS, 3) H2 Fuel Cell, 4) CO Fuel Cell reactions";
  SI.MolarEnthalpy H_h2o[nX]=Medium.MMX[5]*Modelica.Media.IdealGases.Common.Functions.h_T(
      T=Gas[:].T,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.H2O);
  SI.MolarEnthalpy H_h2[nX]=Medium.MMX[1]*Modelica.Media.IdealGases.Common.Functions.h_T(
      T=Gas[:].T,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.H2);
  SI.MolarEnthalpy H_co2[nX]=Medium.MMX[3]*Modelica.Media.IdealGases.Common.Functions.h_T(
      T=Gas[:].T,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.CO2);
  SI.MolarEnthalpy H_co[nX]=Medium.MMX[4]*Modelica.Media.IdealGases.Common.Functions.h_T(
      T=Gas[:].T,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.CO);
equation

  // Energy balance
  Emg = por*dV * Gas.d .* Gas.u;

  Qcond[1] = (kRibs/dx)*(1-por)*lY*lZ.*(T[2] .- T[1]);
  Qcond[2:nX-1] = (kRibs/dx)*(1-por)*lY*lZ.*(T[3:nX] .- 2.*T[2:nX-1] .+ T[1:nX-2]);
  Qcond[nX] = (kRibs/dx)*(1-por)*lY*lZ.*(T[nX-1] .- T[nX]);

  // External heat flow 2D
  for i in 1:nX loop
    QgasExt[i] = sum(Q_PEN.Q[i,:]) .+ sum(Q_IC.Q[i,:]) .+ Qcond[i] .- q_electrochem[i];
  end for;

  Gas.T = T;

  // Reaction kinetics
  for i in 1:nX loop
    massTransfer[i] = sum(PEN_in[i,:].I)./(4*Modelica.Constants.F).*Modelica.Media.IdealGases.SingleGases.O2.data.MM;
    rMSR[i] = kfMSR[i]*(Gas[i].p^2)*(Ycell[i, 2]*Ycell[i, 5] - ((Gas[i].p^2)/Kmsr[i])*(Ycell[i, 1]^3)*(Ycell[i,4]))*dV*por;
    rWGS[i] = kfWGS[i]*(Gas[i].p^2)*(Ycell[i, 4]*Ycell[i, 5] - Ycell[i, 1]*Ycell[i, 3]/Kwgs[i])*dV*por;

    // 2D variables transport
    rHEL[i] = sum(PEN_in[i,:].I_H)./(2*Modelica.Constants.F);
    rCEL[i] = sum(PEN_in[i,:].I_C)./(2*Modelica.Constants.F);

    DeltaG_msr[i] = -252.642810968035.*Gas[i].T .+ 225215.698063031;
    Kmsr[i] = 1.01325^2*1e10*Modelica.Math.exp(-DeltaG_msr[i]./(Modelica.Constants.R.*Gas[i].T));
    kfMSR[i] = A_MSRf*Modelica.Math.exp(-Ea_MSRf./(Modelica.Constants.R.*Gas[i].T));

    DeltaG_wgs[i] = 32.1153.*Gas[i].T .- 3.5211E4;
    Kwgs[i] = Modelica.Math.exp(-DeltaG_wgs[i]./(Modelica.Constants.R.*Gas[i].T));
    kfWGS[i] = A_WGSf.*Modelica.Math.exp(-Ea_WGSf./(Modelica.Constants.R.*Gas[i].T));

    // Mass flow rate calculation
    for j in 1:nSpecies loop
      R[i,j] = Medium.MMX[j].*(a[1,j].*rMSR[i] + a[2,j].*rWGS[i] + a[3,j].*rHEL[i] + a[4,j].*rCEL[i]);
    end for;
  end for;

  // Pressure Drop
  Gas[1].p = infl.p .- dp[1];
  Gas[2:nX].p = Gas[1:nX-1].p .- dp[2:nX];
  outfl.p = Gas[nX].p; // Pressure is upwinded for the last control volume

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                                                       Rectangle(
          extent={{-100,30},{100,14}},
          lineColor={28,108,200},
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={134,134,134})}),                            Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    Documentation(revisions="<html>
</html>", info="<html>
<h2>FuelChannel2DSimp</h2>

<p>
Top-level simplified fuel channel model.
The model extends <code>Channel2DBaseSimp</code> by adding internal
reforming chemistry, electrochemical source terms, and an axial
Darcy-based pressure drop formulation.
</p>

<h3>Spatial Discretisation</h3>
<ul>
<li>1D finite-volume discretisation in streamwise direction (X)</li>
<li>2D thermal and electrochemical coupling to PEN and interconnect (Y)</li>
<li>Uniform gas properties across Y within each X control volume</li>
</ul>
</html>"));
end FuelChannel2DSimp;

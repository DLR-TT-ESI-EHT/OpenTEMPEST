within OpenTEMPEST.SOC.Cell.CrossFlow;
model FuelChannel2D
  extends Channel2DBase(
    channelFac=1,
    redeclare package Medium = OpenTEMPEST.Medium.Fuel_CH4);

  // Kinetic parameters
  parameter Real A_WGSf annotation(Dialog(tab="Kinetics"));
  parameter Real Ea_WGSf annotation(Dialog(tab="Kinetics"));
  parameter Real A_MSRf annotation(Dialog(tab="Kinetics"));
  parameter Real Ea_MSRf annotation(Dialog(tab="Kinetics"));

  // Pressure Drop
  Modelica.SIunits.PressureDifference dp[nX,nY]=Blocks.Functions.pressureDropDarcy(
      mfv[1:nX, :],
      eta,
      Gas.d,
      dx,
      (lZ - hAbsetzen)*dy,
      2.5e-9) "pressure loss";
  Modelica.SIunits.Length hAbsetzen = 0;

  // Kinetics
  constant Real a[4, nSpecies]={{3,-1,0,1,-1,0}, {1,0,1,-1,-1,0}, {-1,0,0,0,1,0}, {0,0,1,-1,0,0}} "Stoichiometric coefficients for 1) MSR, 2) WGS, 3) H2 Fuel Cell, 4) CO Fuel Cell reactions";
  Real kfMSR[nX, nY], kfWGS[nX, nY] "Rate constant for MSR and WGS forward reaction";
  Modelica.SIunits.MolarFlowRate rMSR[nX, nY], rWGS[nX, nY], rEl[nX, nY] "Reaction rates (mol/s) of MSR, WGS and steam electrolysis";
  Modelica.SIunits.MassFlowRate R[nX,nY,nSpecies] "net Rate of production and consumption of products and reactants in the reactor";
  Real Kmsr[nX, nY], Kwgs[nX, nY] "Reaction equillibrium constant of MSR and WGS reaction respectively";
  Modelica.SIunits.MolarEnergy DeltaG_wgs[nX, nY], DeltaG_msr[nX, nY];

  Modelica.SIunits.EnergyFlowRate q_electrochem[nX, nY] =  rEl.*(H_h2 .- H_h2o);

protected
  Modelica.SIunits.MolarEnthalpy H_h2o[nX,nY]=Medium.MMX[5]*Modelica.Media.IdealGases.Common.Functions.h_T(
      T=Gas[:, :].T,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.H2O);
  Modelica.SIunits.MolarEnthalpy H_h2[nX,nY]=Medium.MMX[1]*Modelica.Media.IdealGases.Common.Functions.h_T(
      T=Gas[:, :].T,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.H2);

equation

  // Energy Balance - Thermal Equilibrium between Gas and Solid phase
  QgasExt = Q_PEN.Q .+ Q_IC.Q .- q_electrochem;

  // Reaction kinetics
  rMSR = kfMSR.*(Gas.p.^2).*(Ycell[:,:, 2].*Ycell[:,:, 5] .- ((Gas.p.^2)./Kmsr).*(Ycell[:,:, 1].^3).*(Ycell[:,:,4])) .*dV.*por;
  rWGS = kfWGS.*(Gas.p.^2).*(Ycell[:,:, 4].*Ycell[:,:, 5] .- Ycell[:,:, 1].*Ycell[:,:, 3]./Kwgs) .*dV.*por;
  rEl = PEN_in.I./(2*Modelica.Constants.F);

  DeltaG_msr = -252.642810968035*Gas.T .+ 225215.698063031;
  Kmsr = 1.01325^2*1e10*Modelica.Math.exp(-DeltaG_msr./(Modelica.Constants.R*Gas.T));
  kfMSR = A_MSRf*Modelica.Math.exp(-Ea_MSRf./(Modelica.Constants.R*Gas.T));

  DeltaG_wgs = 32.1153*Gas.T .- 3.5211E4;
  Kwgs = Modelica.Math.exp(-DeltaG_wgs./(Modelica.Constants.R*Gas.T));
  kfWGS = A_WGSf*Modelica.Math.exp(-Ea_WGSf./(Modelica.Constants.R*Gas.T));

  for j in 1:nSpecies loop
    R[:,:, j] =  Medium.MMX[j]*(a[3, j]*rEl);
  end for;

  // Pressure Drop
  // flow channel area reduction due to degradation
 // time>720000 and ASROhmDegrationRate>0 then

  Gas[1, :].p = infl[:].p .- dp[1, :];
  Gas[2:nX, :].p = Gas[1:nX-1, :].p .- dp[2:nX, :];
  outfl[:].p = Gas[nX, :].p; // Pressure is upwinded for the last control volume

end FuelChannel2D;

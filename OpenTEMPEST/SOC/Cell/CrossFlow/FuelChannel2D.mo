within OpenTEMPEST.SOC.Cell.CrossFlow;
model FuelChannel2D
  extends Channel2DBase(
    redeclare package Medium = OpenTEMPEST.Medium.Fuel_CH4);

  import SI = Modelica.SIunits;

  // Kinetic parameters
  constant Real A_WGSf = 0.0171 annotation(Dialog(tab="Kinetics"));
  constant Real Ea_WGSf = 103191 annotation(Dialog(tab="Kinetics"));
  constant Real A_MSRf = 2395 annotation(Dialog(tab="Kinetics"));
  constant Real Ea_MSRf = 231266*0.875 annotation(Dialog(tab="Kinetics"));

  // Pressure Drop
  SI.PressureDifference dp[nX,nY]=Blocks.Functions.pressureDropDarcy(
      mfv[1:nX, :],
      eta,
      Gas.d,
      dx,
      (lZ - hAbsetzen)*dy,
      3e-9) "pressure loss";
  SI.Length hAbsetzen = 0;

  // Kinetics
  Real kfMSR[nX, nY], kfWGS[nX, nY] "Rate constant for MSR and WGS forward reaction";
  SI.MolarFlowRate rMSR[nX, nY], rWGS[nX, nY], rHEL[nX, nY], rCEL[nX, nY] "Reaction rates (mol/s) of MSR, WGS and steam electrolysis";
  SI.MassFlowRate R[nX,nY,nSpecies] "net rate of production and consumption of products and reactants in the reactor";
  Real Kmsr[nX, nY], Kwgs[nX, nY] "Reaction equillibrium constant of MSR and WGS reaction respectively";
  SI.MolarEnergy DeltaG_wgs[nX, nY], DeltaG_msr[nX, nY];

  SI.EnergyFlowRate Qcond[nX, nY];
  SI.EnergyFlowRate q_electrochem[nX, nY] =  rHEL.*(H_h2 .- H_h2o) .+ rCEL.*(H_co .- H_co2);

protected
  constant Real a[4, nSpecies]={{3,-1,0,1,-1,0}, {1,0,1,-1,-1,0}, {-1,0,0,0,1,0}, {0,0,1,-1,0,0}} "Stoichiometric coefficients for 1) MSR, 2) WGS, 3) H2 Fuel Cell, 4) CO Fuel Cell reactions";
  SI.MolarEnthalpy H_h2o[nX,nY]=Medium.MMX[5]*Modelica.Media.IdealGases.Common.Functions.h_T(
      T=Gas[:, :].T,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.H2O);
  SI.MolarEnthalpy H_h2[nX,nY]=Medium.MMX[1]*Modelica.Media.IdealGases.Common.Functions.h_T(
      T=Gas[:, :].T,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.H2);
  SI.MolarEnthalpy H_co2[nX,nY]=Medium.MMX[3]*Modelica.Media.IdealGases.Common.Functions.h_T(
      T=Gas[:, :].T,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.CO2);
  SI.MolarEnthalpy H_co[nX,nY]=Medium.MMX[4]*Modelica.Media.IdealGases.Common.Functions.h_T(
      T=Gas[:, :].T,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.CO);

equation
  // Energy balance
  Emg = por*dV * Gas.d .* Gas.u;

  Qcond[1,:] = (kRibs/dx)*(1-por)*dy*lZ.*(T[2,:] .- T[1,:]);
  Qcond[2:nX-1,:] = (kRibs/dx)*(1-por)*dy*lZ.*(T[3:nX,:] .- 2.*T[2:nX-1,:] .+ T[1:nX-2,:]);
  Qcond[nX,:] = (kRibs/dx)*(1-por)*dy*lZ.*(T[nX-1,:] .- T[nX,:]);

  // Energy Balance
  QgasExt = Q_PEN.Q .+ Q_IC.Q .+ Qcond .- q_electrochem;

  Gas.T = T;

  // Reaction kinetics
  massTransfer = PEN_in.I./(4*Modelica.Constants.F).*Modelica.Media.IdealGases.SingleGases.O2.data.MM;
  rMSR = kfMSR.*(Gas.p.^2).*(Ycell[:,:, 2].*Ycell[:,:, 5] .- ((Gas.p.^2)./Kmsr).*(Ycell[:,:, 1].^3).*(Ycell[:,:,4])) .*dV.*por;
  rWGS = kfWGS.*(Gas.p.^2).*(Ycell[:,:, 4].*Ycell[:,:, 5] .- Ycell[:,:, 1].*Ycell[:,:, 3]./Kwgs) .*dV.*por;
  rHEL = PEN_in.I_H./(2*Modelica.Constants.F);
  rCEL = PEN_in.I_C./(2*Modelica.Constants.F);

  DeltaG_msr = -252.642810968035*Gas.T .+ 225215.698063031;
  Kmsr = 1.01325^2*1e10*Modelica.Math.exp(-DeltaG_msr./(Modelica.Constants.R*Gas.T));
  kfMSR = A_MSRf*Modelica.Math.exp(-Ea_MSRf./(Modelica.Constants.R*Gas.T));

  DeltaG_wgs = 32.1153*Gas.T .- 3.5211E4;
  Kwgs = Modelica.Math.exp(-DeltaG_wgs./(Modelica.Constants.R*Gas.T));
  kfWGS = A_WGSf*Modelica.Math.exp(-Ea_WGSf./(Modelica.Constants.R*Gas.T));

  for i in 1:nX loop
    for j in 1:nY loop
      for k in 1:nSpecies loop
        R[i,j, k] =  Medium.MMX[k]*(a[1, k]*rMSR[i,j] + a[2, k]*rWGS[i,j] + a[3, k]*rHEL[i,j] + a[4, k]*rCEL[i,j]);
      end for;
    end for;
  end for;

  // Pressure drop
  Gas[1, :].p = infl[:].p .- dp[1, :];
  Gas[2:nX, :].p = Gas[1:nX-1, :].p .- dp[2:nX, :];
  outfl[:].p = Gas[nX, :].p; // Pressure is upwinded for the last control volume

  annotation (Documentation(info="<html>
<h2>FuelChannel2D</h2>

<p>
Two-dimensional finite-volume fuel channel model extending 
<code>Channel2DBase</code>. The model describes coupled mass, energy, 
and momentum balances in a porous fuel channel discretized in the 
streamwise (x) and transverse (y) directions, including reforming 
and electrochemical reactions.
</p>

<h3>Modeling Assumptions</h3>
<ul>
<li>2D discretization in x–y plane</li>
<li>Primary flow in x-direction with distributed pressure losses</li>
<li>Convective heat exchange with PEN and interconnect</li>
<li>Axial solid rib conduction in x-direction</li>
</ul>

<h3>Included Phenomena</h3>
<ul>
<li>Methane Steam Reforming (MSR)</li>
<li>Water-Gas Shift (WGS)</li>
<li>Hydrogen electrochemical reaction</li>
<li>Carbon dioxide electrochemical reaction</li>
<li>Darcy-based pressure drop formulation</li>
</html>"));
end FuelChannel2D;

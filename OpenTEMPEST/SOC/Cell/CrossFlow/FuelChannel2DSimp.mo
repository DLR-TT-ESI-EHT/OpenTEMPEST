within OpenTEMPEST.SOC.Cell.CrossFlow;
model FuelChannel2DSimp "Fuel channel model 1D discretised with 2D interfaces"

  import SI = Modelica.SIunits;

  extends Channel2DBaseSimp(
    Nu_IC=9.86,
    Nu_PEN=12,
    channelFac=1,
    redeclare package Medium = OpenTEMPEST.Medium.Fuel_CH4);

  // Selections
  parameter Boolean LUDS=false "Set true if Linear Upwind difference wanted (more accuracy), false for Upwind difference scheme (more speed and stability)";
  parameter Boolean heatTransferCorrelationFormDuct=true "true for Nusselt correlation duct geometry with characteristic length=2*lZ (default), false for plate geometry with characteristic length=lX";

  // Kinetic parameters
  constant Real A_sf = 0.0171;
  constant Real Ea_sf = 103191;
  constant Real A_rf = 2395;
  constant Real Ea_rf = 231266*0.875;

  // Ni Foam Thermophysical parameters
  parameter SI.ThermalConductivity k_eff_Ni=7.971  "Nickel foam effective thermal conductivity";
  parameter SI.SpecificHeatCapacity cp_Ni=455  "Nickel heat capacity";
  parameter SI.Density rho_Ni=4066/(1-0.5437)  "Nickel density";

  // Pressure Drop
  SI.PressureDifference dp[nX] = Blocks.Functions.pressureDropDarcy(         mfv[1:nX], eta, Gas.d, dx, (lZ - hAbsetzen)*lY, 2.5e-9) "pressure loss";
//   parameter SI.PressureDifference dpNom=0.01 "nominal pressure loss for initialization";
  SI.DynamicViscosity eta[nX] = Medium.dynamicViscosity(Gas.state);
  SI.Length hAbsetzen = 0;

  // Kinetics
  Real kfMSR[nX], kfWGS[nX];
  SI.MolarFlowRate rMSR[nX], rWGS[nX], rHEL[nX], rCEL[nX] "Reaction rates (mol/s) of MSR, WGS and steam and CO2 electrolysis";
  SI.MassFlowRate R[nX, nSpecies] "Net rate of production and consumption of products and reactants in the reactor";
  Real Kmsr[nX], Kwgs[nX] "Reaction equillibrium constant of reverse reforming and WGS reaction respectively";
  SI.MolarEnergy DeltaG_wgs[nX], DeltaG_msr[nX];

  SI.EnergyFlowRate q_electrochem[nX] =  rHEL.*(H_h2 .- H_h2o) .+ rCEL.*(H_co .- H_co2);

protected
    constant Real a[4, nSpecies]={{3,-1,0,1,-1,0}, {1,0,1,-1,-1,0}, {-1,0,0,0,1,0}, {0,0,1,-1,0,0}} "Stoichiometric coefficients for 1) MSR, 2) WGS, 3) H2 Fuel Cell, 4) CO Fuel Cell reactions";
  Modelica.SIunits.MolarEnthalpy H_h2o[nX]=Medium.MMX[5]*Modelica.Media.IdealGases.Common.Functions.h_T(
      T=Gas[:].T,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.H2O);
  Modelica.SIunits.MolarEnthalpy H_h2[nX]=Medium.MMX[1]*Modelica.Media.IdealGases.Common.Functions.h_T(
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

  // External heat flow 2D
  for i in 1:nX loop
    QgasExt[i] = sum(Q_PEN.Q[i,:]) .+ sum(Q_IC.Q[i,:]) .- q_electrochem[i];
  end for;

  // Reaction kinetics
  for i in 1:nX loop
    rMSR[i] = kfMSR[i]*(Gas[i].p^2)*(Ycell[i, 2]*Ycell[i, 5] - ((Gas[i].p^2)/Kmsr[i])*(Ycell[i, 1]^3)*(Ycell[i,4]))*dV*por;
    rWGS[i] = kfWGS[i]*(Gas[i].p^2)*(Ycell[i, 4]*Ycell[i, 5] - Ycell[i, 1]*Ycell[i, 3]/Kwgs[i])*dV*por;

    // 2D variables transport
    rHEL[i] = sum(PEN_in[i,:].I_H)./(2*Modelica.Constants.F);
    rCEL[i] = sum(PEN_in[i,:].I_C)./(2*Modelica.Constants.F);

    DeltaG_msr[i] = -252.642810968035.*Gas[i].T .+ 225215.698063031;
    Kmsr[i] = 1.01325^2*1e10*Modelica.Math.exp(-DeltaG_msr[i]./(Modelica.Constants.R.*Gas[i].T));
    kfMSR[i] = A_rf*Modelica.Math.exp(-Ea_rf./(Modelica.Constants.R.*Gas[i].T));

    DeltaG_wgs[i] = 32.1153.*Gas[i].T .- 3.5211E4;
    Kwgs[i] = Modelica.Math.exp(-DeltaG_wgs[i]./(Modelica.Constants.R.*Gas[i].T));
    kfWGS[i] = A_sf.*Modelica.Math.exp(-Ea_sf./(Modelica.Constants.R.*Gas[i].T));

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
<ul>
<li><i>2 Feb 2025</i>
    by Anis Taissir</a>:<br>
       First release.</li>
</ul>
</html>"));
end FuelChannel2DSimp;

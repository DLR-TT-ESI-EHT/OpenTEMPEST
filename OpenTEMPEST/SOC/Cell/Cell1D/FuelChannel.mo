within OpenTEMPEST.SOC.Cell.Cell1D;
model FuelChannel
  extends OpenTEMPEST.SOC.Cell.Cell1D.Channel1DBase(
    redeclare package Medium = OpenTEMPEST.Medium.Fuel_CH4);

  import SI = Modelica.SIunits;
  // Thermal
  SI.Temperature T[N](start = TStart) "Average Temperature in CV";

  // Ni Foam Thermophysical parameters
  parameter SI.ThermalConductivity k_eff_Ni  "Nickel foam effective thermal conductivity";
  parameter SI.SpecificHeatCapacity cp_Ni  "Nickel heat capacity";
  parameter SI.Density rho_Ni  "Nickel density";

  // Kinetic parameters
  constant Real A_WGSf = 0.0171  annotation(Dialog(tab="Kinetics"));
  constant Real Ea_WGSf = 103191  annotation(Dialog(tab="Kinetics"));
  constant Real A_MSRf = 2395  annotation(Dialog(tab="Kinetics"));
  constant Real Ea_MSRf = 231266*0.875  annotation(Dialog(tab="Kinetics"));

  // Pressure Drop
  SI.PressureDifference dp[N] "pressure loss";
  parameter SI.PressureDifference dpNom=0.01 "nominal pressure loss for initialization";
  SI.Length hAbsetzen;

  // Kinetics
  Real kfMSR[N], kfWGS[N] "Rate constant for MSR and WGS forward reaction";
  SI.MolarFlowRate rMSR[N], rWGS[N], rHEL[N], rCEL[N] "Reaction rates (mol/s) of MSR, WGS and steam and CO2 electrolysis";
  SI.MassFlowRate R[N, nSpecies]
    "net Rate of production and consumption of products and reactants in the reactor";
  Real Kmsr[N], Kwgs[N]
    "Reaction equillibrium constant of reverse reforming and WGS reaction respectively";
  SI.MolarEnergy DeltaG_wgs[N], DeltaG_msr[N];
  SI.EnergyFlowRate Qcond[N];
  SI.EnergyFlowRate q_electrochem[N] =  rHEL.*(H_h2 .- H_h2o) .+ rCEL.*(H_co .- H_co2);

  ThermoPower.Thermal.HT hTNi_x0 "conductive HT through Ni at x=0" annotation (Placement(transformation(extent={{-104,
            18},{-84,38}}), iconTransformation(extent={{-104,18},{-84,38}})));
  ThermoPower.Thermal.HT hTNi_xN "conductive HT through Ni at x=L" annotation (Placement(transformation(extent={{82,
            20},{102,40}}), iconTransformation(extent={{82,20},{102,40}})));

protected
  constant Real a[4, nSpecies]={{3,-1,0,1,-1,0}, {1,0,1,-1,-1,0}, {-1,0,0,0,1,0}, {0,0,1,-1,0,0}} "Stoichiometric coefficients for 1) MSR, 2) WGS, 3) H2 Fuel Cell, 4) CO Fuel Cell reactions";

  SI.MolarEnthalpy H_h2o[N]=Medium.MMX[5]*Modelica.Media.IdealGases.Common.Functions.h_T(
      T=Q_PEN.T,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.H2O);
  SI.MolarEnthalpy H_h2[N]=Medium.MMX[1]*Modelica.Media.IdealGases.Common.Functions.h_T(
      T=Q_PEN.T,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.H2);
  SI.MolarEnthalpy H_co2[N]=Medium.MMX[3]*Modelica.Media.IdealGases.Common.Functions.h_T(
      T=Q_PEN.T,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.CO2);
  SI.MolarEnthalpy H_co[N]=Medium.MMX[4]*Modelica.Media.IdealGases.Common.Functions.h_T(
      T=Q_PEN.T,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.CO);

equation

// Energy Balance - Thermal Equilibrium between Gas and Solid phase
  Emg[:] = Gas[:].d.*Gas[:].u.*por*dV .+ rho_Ni.*(1-por)*dV.*cp_Ni.*T[:];

  Qcond[1] = (k_eff_Ni/dx)*lY*lZ*(T[2] - T[1]) + hTNi_x0.Q_flow;
  Qcond[2:N-1] = (k_eff_Ni/dx)*lY*lZ.*(T[3:N] .- 2.*T[2:N-1] .+ T[1:N-2]);
  Qcond[N] = (k_eff_Ni/dx)*lY*lZ*(T[N-1] - T[N]) + hTNi_xN.Q_flow;

   Gas[:].T = T[:];

   QgasExt[:] = Q_PEN.Q[:] .+ Q_IC.Q[:] .+ Qcond[:] .- q_electrochem[:];
   hTNi_x0.Q_flow = (k_eff_Ni/(0.5*dx))*lY*lZ*(hTNi_x0.T - T[1]);
   hTNi_xN.Q_flow = (k_eff_Ni/(0.5*dx))*lY*lZ*(hTNi_xN.T - T[N]);

  // Reaction kinetics
  massTransfer[:] = PEN_in[:].I/(4*Modelica.Constants.F)*Modelica.Media.IdealGases.SingleGases.O2.data.MM;
  for i in 1:N loop
    rMSR[i] = kfMSR[i]*(Gas[i].p^2)*(Ycell[i, 2]*Ycell[i, 5] - ((Gas[i].p^2)/Kmsr[i])*(Ycell[i, 1]^3)*(Ycell[i,4]))*dV*por;
    rWGS[i] = kfWGS[i]*(Gas[i].p^2)*(Ycell[i, 4]*Ycell[i, 5] - Ycell[i, 1]*Ycell[i, 3]/Kwgs[i])*dV*por;
    rHEL[i] = PEN_in[i].I_H/(2*Modelica.Constants.F);
    rCEL[i] = PEN_in[i].I_C/(2*Modelica.Constants.F);

    DeltaG_msr[i] = -252.642810968035.*Gas[i].T + 225215.698063031;
    Kmsr[i] = 1.01325^2*1e10*Modelica.Math.exp(-DeltaG_msr[i]/Modelica.Constants.R/Gas[i].T);
    kfMSR[i] = A_MSRf*Modelica.Math.exp(-Ea_MSRf/(Modelica.Constants.R*Gas[i].T));

    DeltaG_wgs[i] = 32.1153*(Gas[i].T) - 3.5211E4; // Marius' shortcut from NASA
    Kwgs[i] = Modelica.Math.exp(-DeltaG_wgs[i]/Modelica.Constants.R/Gas[i].T);
    kfWGS[i] = A_WGSf*Modelica.Math.exp(-Ea_WGSf/(Modelica.Constants.R*Gas[i].T));

    for j in 1:nSpecies loop
      R[i, j] =  Medium.MMX[j]*(a[1, j]*rMSR[i] .+ a[2, j]*rWGS[i] .+ a[3, j]*rHEL[i] .+ a[4, j]*rCEL[i]);
    end for;
  end for;

  // flow channel area reduction due to degradation
  if false then // time>720000 and ASROhmDegrationRate>0 then
    hAbsetzen=0.1*lZ;
  else
    hAbsetzen=0;
  end if;

  // Pressure Drop
  for i in 1:N loop
    dp[i] =OpenTEMPEST.Blocks.Functions.pressureDropDarcy(
      mfv[i],
      eta[i],
      Gas[i].d,
      dx,
      (lZ - hAbsetzen)*lY,
      3e-9);
  end for;

  Gas[1].p = infl.p - dp[1];
  Gas[2:N].p = Gas[1:N-1].p - dp[2:N];
  outfl.p = Gas[N].p; // Pressure is upwinded for the last control volume

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
          extent={{-84,34},{82,-28}},
          lineColor={28,108,200},
          fillColor={28,108,200},
          fillPattern=FillPattern.HorizontalCylinder), Rectangle(
          extent={{-84,28},{82,12}},
          lineColor={28,108,200},
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={134,134,134})}),                            Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<h2>FuelChannel1D</h2>

<p>
One-dimensional finite-volume fuel channel model extending 
<code>Channel1DBase</code>. The model describes coupled mass, energy, 
and momentum balances in a porous nickel foam channel including 
heterogeneous reforming reactions and electrochemical source terms.
</p>

<h3>Modeling Assumptions</h3>
<ul>
<li>1D discretization in flow direction</li>
<li>Thermal equilibrium between gas and nickel foam solid phase</li>
<li>Effective axial heat conduction in the solid matrix</li>
<li>Darcy-based pressure drop formulation</li>
</ul>

<h3>Included Phenomena</h3>
<ul>
<li>Methane Steam Reforming (MSR)</li>
<li>Water-Gas Shift (WGS)</li>
<li>Hydrogen electrochemical reaction</li>
<li>Carbon dioxide electrochemical reaction</li>
<li>Axial solid conduction through nickel foam</li>
<li>Distributed pressure losses</li>
</html>"));
end FuelChannel;

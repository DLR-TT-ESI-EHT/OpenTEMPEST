within OpenTEMPEST.SOC.Cell.CrossFlow;
model AirChannel2DSimp "Air channel model 1D discretised with 2D interfaces"

  import SI = Modelica.SIunits;

  extends Channel2DBaseSimp(
    redeclare package Medium = OpenTEMPEST.Medium.Air_Medium);

  // Pressure Drop
  SI.PressureDifference dp[nX]=Blocks.Functions.pressureDropDarcy(mfv[1:nX], eta, Gas.d, dx, lZ*lY, 3e-9) "pressure loss";
  parameter Real pDrop(max=0.99) "pressure drop as a factor of inlet pressure (between 0 and 0.99)";
  SI.DynamicViscosity eta[nX]=Medium.dynamicViscosity(Gas.state);

  // Kinetics
  SI.MolarFlowRate rEl[nX];

  SI.EnergyFlowRate Qcond[nX];
  SI.EnergyFlowRate q_electrochem[nX] = rEl .* (0.5*H_o2);

protected
  Real a[2] = {-0.5, 0};
  SI.MolarEnthalpy H_o2[nX]=Medium.MMX[1]*Modelica.Media.IdealGases.Common.Functions.h_T(
    T=Gas[:].T,
    exclEnthForm=false,
    refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
    data=Modelica.Media.IdealGases.Common.SingleGasesData.O2);

equation
  // Energy Balance
  Emg[:] = Gas[:].d.*Gas[:].u.*por*dV;

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
    massTransfer[i] = -sum(PEN_in[i,:].I)./(4*Modelica.Constants.F).*Modelica.Media.IdealGases.SingleGases.O2.data.MM;
    rEl[i] = sum(PEN_in[i,:].I)./(2*Modelica.Constants.F); // 2D variables stream interface

    // Mass flow rate calculation
    for j in 1:nSpecies loop
      R[i, j] = a[j].*rEl[i].*Medium.MMX[j];
    end for;
  end for;

  // Pressure Drop
  //dp~momentum Balance
  Gas[1].p = infl.p - dp[1];
  Gas[2:nX].p = Gas[1:nX-1].p .- dp[2:nX];
  outfl.p = Gas[nX].p; // Pressure is upwinded for the last control volume

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                                                       Rectangle(
          extent={{-100,28},{100,12}},
          lineColor={28,108,200},
          fillPattern=FillPattern.CrossDiag,
          fillColor={134,134,134})}),                            Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    Documentation(revisions="<html>
</html>", info="<html>
<h2>AirChannel2DSimp</h2>

<p>
Top-level air channel simplified model.
The model extends <code>Channel2DBaseSimp</code> and represents the
oxygen electrode side including electrochemical oxygen
reduction/evolution and axial pressure losses.
</p>

<h3>Spatial Discretisation</h3>
<ul>
<li>1D finite-volume discretisation in streamwise direction (X)</li>
<li>2D thermal and electrochemical coupling to PEN and interconnect (Y)</li>
<li>Uniform gas properties across Y within each X control volume</li>
</ul>
</html>"));
end AirChannel2DSimp;

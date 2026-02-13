within OpenTEMPEST.SOC.Cell.CrossFlow;
model AirChannel2D
  extends Channel2DBase(
    redeclare package Medium = OpenTEMPEST.Medium.Air_Medium,
    channelFac=-1);

  // Pressure Drop
  Modelica.SIunits.PressureDifference dp[nX,nY]=Blocks.Functions.pressureDropDarcy(
      mfv[1:nX, :],
      eta,
      Gas.d,
      dx,
      lZ*dy,
      2.5e-9) "pressure loss";
  parameter Real pDrop(max=0.99) "pressure drop as a factor of inlet pressure (between 0 and 0.99)";

  // Kinetics
  Modelica.SIunits.MolarFlowRate rEl[nX, nY];
  Modelica.SIunits.EnergyFlowRate q_electrochem[nX, nY] = rEl .* (0.5*H_o2);

protected
  constant Real a[2] = {-0.5, 0};
  Modelica.SIunits.MolarEnthalpy H_o2[nX,nY]=Medium.MMX[1]*Modelica.Media.IdealGases.Common.Functions.h_T(
      T=Gas[:, :].T,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.O2);

equation

  // Energy Balance
  QgasExt = Q_PEN.Q .+ Q_IC.Q .- q_electrochem;

  // Reaction kinetics
  rEl = PEN_in.I/(2*Modelica.Constants.F);
  for j in 1:nSpecies loop
    R[:, :, j] = a[j]*rEl*Medium.MMX[j];
  end for;

  // Pressure Drop
  //dp~momentum Balance
  Gas[1,:].p = infl[:].p .- dp[1,:];
  Gas[2:nX,:].p = Gas[1:nX-1,:].p .- dp[2:nX,:];
  outfl.p = Gas[nX,:].p; // Pressure is upwinded for the last control volume

end AirChannel2D;

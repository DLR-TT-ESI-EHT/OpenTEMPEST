within OpenTEMPEST.SOC.Cell.CrossFlow;
model AirChannel2DSimp "Air channel model 1D discretised with 2D interfaces"

  import SI = Modelica.SIunits;

  extends Channel2DBaseSimp(
    Nu_IC=7.54,
    Nu_PEN=8.235,
    redeclare package Medium = OpenTEMPEST.Medium.Air_Medium,
    channelFac=-1);

  // Selections
  parameter Boolean LUDS=false "Set true if Linear Upwind difference wanted (more accuracy), false for Upwind difference scheme (more speed and stability)";
  parameter Boolean heatTransferCorrelationFormDuct=true "true for Nusselt correlation duct geometry with characteristic length=2*lZ (default), false for plate geometry with characteristic length=lX";

  // Pressure Drop
  SI.PressureDifference dp[nX]=Blocks.Functions.pressureDropDarcy(        mfv[1:nX], eta, Gas.d, dx, lZ*lY, 2.5e-9) "pressure loss";
  parameter Real pDrop(max=0.99) "pressure drop as a factor of inlet pressure (between 0 and 0.99)";
  SI.DynamicViscosity eta[nX]=Medium.dynamicViscosity(Gas.state);

  // Kinetics
  SI.MolarFlowRate rEl[nX];
  SI.EnergyFlowRate q_electrochem[nX] = rEl .* (0.5*H_o2);
protected
  Real a[2] = {-0.5, 0};
  SI.MolarEnthalpy H_o2[nX]=Medium.MMX[1]*Modelica.Media.IdealGases.Common.Functions.h_T(
    T=Gas[:].T,
    exclEnthForm=false,
    refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
    data=Modelica.Media.IdealGases.Common.SingleGasesData.O2);

equation

   // External heat flow 2D
   for i in 1:nX loop
     QgasExt[i] = sum(Q_PEN.Q[i,:]) .+ sum(Q_IC.Q[i,:]) .- q_electrochem[i];
   end for;

  // Reaction kinetics
  for i in 1:nX loop
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
<ul>
<li><i>2 Feb 2025</i> by Anis Taissir</a>:<br>
       First release.</li>
</ul>
</html>"));
end AirChannel2DSimp;

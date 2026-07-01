within OpenTEMPEST.SOC.Cell.Cell1D;
model AirChannel
  extends OpenTEMPEST.SOC.Cell.Cell1D.Channel1DBase(
    pStartOut=(1 - pDrop)*pStartIn,
    redeclare package Medium = OpenTEMPEST.Medium.Air_Medium);

  import SI = Modelica.SIunits;

  // Pressure Drop
  SI.PressureDifference dp[N] "pressure loss";
  parameter Real pDrop(max=0.99) "pressure drop as a factor of inlet pressure (between 0 and 0.99)";

  // Kinetics
  SI.MolarFlowRate rEl[N];
  SI.EnergyFlowRate q_electrochem[N] = rEl .* (0.5*H_o2);

protected
  Real a[nSpecies] = {-0.5, 0};

  SI.MolarEnthalpy H_o2[N]=Medium.MMX[1]*Modelica.Media.IdealGases.Common.Functions.h_T(
      T=Q_PEN.T,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.O2);

equation

  // Energy Balance
  Emg[:] = Gas[:].d.*Gas[:].u.*por*dV;
  QgasExt[:] = Q_PEN.Q[:] + Q_IC.Q[:] .- q_electrochem[:];

  // Reaction kinetics
  massTransfer[:] = -PEN_in[:].I/(4*Modelica.Constants.F)*Modelica.Media.IdealGases.SingleGases.O2.data.MM;
  for i in 1:N loop
    rEl[i] = PEN_in[i].I/(2*Modelica.Constants.F);
    for j in 1:nSpecies loop
      R[i, j] = a[j]*rEl[i]*Medium.MMX[j];
    end for;
  end for;

  // Pressure Drop
  //dp~momentum Balance
  for i in 1:N loop
    dp[i] = pDrop/N*infl.p;
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
          fillPattern=FillPattern.CrossDiag,
          fillColor={134,134,134})}),                            Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<h2>AirChannel</h2>

<p>
One-dimensional finite-volume air channel model extending 
<code>Channel1DBase</code>. The model describes coupled mass, energy, 
and momentum balances in the channel including oxygen electrochemical reaction 
and axial pressure drop.
</p>

<h3>Included Phenomena</h3>
<ul>
<li>Oxygen electrochemical reaction via <code>PEN_in</code></li>
<li>Electrochemical heat source term</li>
<li>Darcy-based axial pressure drop</li>
<li>Convective heat exchange with PEN and interconnect</li>
</ul>
</html>"));
end AirChannel;

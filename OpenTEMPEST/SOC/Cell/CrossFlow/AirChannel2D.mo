within OpenTEMPEST.SOC.Cell.CrossFlow;
model AirChannel2D
  extends Channel2DBase(
    redeclare package Medium = OpenTEMPEST.Medium.Air_Medium);

  import SI = Modelica.SIunits;

  // Pressure Drop
  SI.PressureDifference dp[nX,nY]=Blocks.Functions.pressureDropDarcy(
      mfv[1:nX, :],
      eta,
      Gas.d,
      dx,
      lZ*dy,
      3e-9) "pressure loss";
  parameter Real pDrop(max=0.99) "pressure drop as a factor of inlet pressure (between 0 and 0.99)";

  // Kinetics
  SI.MolarFlowRate rEl[nX, nY];

  SI.EnergyFlowRate Qcond[nX, nY];
  SI.EnergyFlowRate q_electrochem[nX, nY] = rEl .* (0.5*H_o2);

protected
  constant Real a[2] = {-0.5, 0};
  SI.MolarEnthalpy H_o2[nX,nY]=Medium.MMX[1]*Modelica.Media.IdealGases.Common.Functions.h_T(
      T=Gas[:, :].T,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.O2);

equation
  // Energy balance
  Emg = por*dV * Gas.d .* Gas.u;

  Qcond[1,:] = (kRibs/dx)*(1-por)*dy*lZ.*(T[2,:] .- T[1,:]);
  Qcond[2:nX-1,:] = (kRibs/dx)*(1-por)*dy*lZ.*(T[3:nX,:] .- 2.*T[2:nX-1,:] .+ T[1:nX-2,:]);
  Qcond[nX,:] = (kRibs/dx)*(1-por)*dy*lZ.*(T[nX-1,:] .- T[nX,:]);

  QgasExt = Q_PEN.Q .+ Q_IC.Q .+ Qcond .- q_electrochem;

  Gas.T = T;

  // Reaction kinetics
  massTransfer = -PEN_in.I ./ (4*Modelica.Constants.F).*Modelica.Media.IdealGases.SingleGases.O2.data.MM;
  rEl = PEN_in.I/(2*Modelica.Constants.F);
  for j in 1:nSpecies loop
    R[:, :, j] = a[j]*rEl[:,:]*Medium.MMX[j];
  end for;

  // Pressure Drop
  //dp~momentum Balance
  Gas[1,:].p = infl[:].p .- dp[1,:];
  Gas[2:nX,:].p = Gas[1:nX-1,:].p .- dp[2:nX,:];
  outfl.p = Gas[nX,:].p; // Pressure is upwinded for the last control volume

  annotation (Documentation(info="<html>
<h2>AirChannel2D</h2>


<p>
Two-dimensional finite-volume air channel model extending 
<code>Channel2DBase</code>. The model describes coupled mass, energy, 
and momentum balances in an air channel discretized in the 
streamwise (x) and transverse (y) directions, including oxygen 
electrochemical reaction and distributed pressure losses.
</p>

<h3>Modeling Assumptions</h3>
<ul>
<li>2D discretization in x–y plane</li>
<li>Primary flow in x-direction</li>
<li>Convective heat exchange with PEN and interconnect</li>
<li>Axial solid rib conduction in x-direction</li>
</ul>

<h3>Included Phenomena</h3>
<ul>
<li>Oxygen electrochemical reaction coupled via <code>PEN_in</code></li>
<li>Electrochemical heat source term</li>
<li>Darcy-based pressure drop formulation</li>
</ul>
</html>"));
end AirChannel2D;

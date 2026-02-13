within OpenTEMPEST.Flow;
model Flow1D "1D flow model FVM, uses UDS and LUDS"
  extends Modelica.Icons.UnderConstruction; // use flux interpolator and remove comments
  extends OpenTEMPEST.Flow.BaseClasses.Flow1DBaseClass;
  import SI = Modelica.SIunits;

  parameter Real dpdxNom = 100 "Pressure drop in pipe Pa/m";
  parameter Boolean LUDS=false "Set true for linear upwind differencing";
  parameter Boolean internalHT=true "set true if you want to model HT phenomena inside the channel, false for HT modelled outside";

  parameter SI.CoefficientOfHeatTransfer h_conv=1 "Heat Transfer Coefficient /W/m2/K" annotation(dialog(enable=internalHT));

  // Vertex Values
  SI.MassFlowRate mfv[N+1] "Mass flow rate in and leaving CV";
  SI.SpecificEnthalpy hv[N+1] "specific enthalpy at CV boundary";
  SI.MassFraction xiv[N+1, nSpecies] "Mass fractions in CV";

  // External Heat Transfer
  Real Qext[N] = wall.Q./Ntubes;
  Real Qtot;

//   SI.SpecificHeatCapacity cp[N] = Medium.heatCapacity_cp(Gas[1:N].state);
//   SI.SpecificHeatCapacity cpm = sum(cp)/N;

//   SI.ReynoldsNumber Re[N];
//   SI.DynamicViscosity eta[N];
//   SI.PrandtlNumber Pr[N];

  ThermoPower.Thermal.DHTVolumes wall(N=N)
                                       "HT through cylindrical wall" annotation (Placement(
        transformation(extent={{-80,30},{80,50}}),iconTransformation(extent={{-80,34},
            {80,50}})));

equation

  // Mass balance - single channel
  dV.*theta[:,1].*der(Gas[:].d) = (mfv[1:N] .- mfv[2:N+1])./Ntubes;

//   for i in 1:N loop
//     Re[i] = (2*Rh*mf[i]/Ntubes)/(Acs*eta[i]); // Re in single channel CV
//     eta[i] = Medium.dynamicViscosity(Gas[i].state);
//     Pr[i] = Medium.prandtlNumber(Gas[i].state);
//   end for;

  // Species Mass balance - single channel
  for i in 1:nSpecies loop
    dV.*theta[:,1].*der(Gas[:].d.*Gas[:].Xi[i]) = (mfv[1:N].*xiv[1:N, i] .- mfv[2:N+1].*xiv[2:N+1, i])./Ntubes;
  end for;

  // Energy Balance - single channel
  Ntubes*dV.*theta[:,1].*der(Gas[:].d.*Gas[:].u) = (mfv[1:N].*hv[1:N] .- mfv[2:N+1].*hv[2:N+1]) .+ wall.Q[:];

  if internalHT then
    wall.Q[:] = h_conv.*Ntubes.*Aw.*(wall.T[:] - Gas[:].T);
  else
    wall.T[:] = Gas[:].T;
  end if;
  Qtot = sum(wall.Q);

  // Momentum Balance
  0 = (infl.p - Gas[1].p) - dpdxNom*dx;
  fill(0, N-1) = (Gas[1:N-1].p .- Gas[2:N].p) .- dpdxNom*dx;

  // Interface interpolation
    // At x=0
  mfv[1] = infl.m_flow;
  xiv[1, :] = infl.Xi_outflow[:];
  xiv[1, :] = inStream(infl.Xi_outflow[:]);
  hv[1] = infl.h_outflow;
  hv[1] = inStream(infl.h_outflow);

  if LUDS then
    // 1st CV
    mfv[2] = (3*mf[1] - mfv[1])/2;
    hv[2] = (3*Gas[1].h - hv[1])/2;
    // Remaining CVs
    mfv[3:N+1] = (3*mf[2:N] - mf[1:N-1])/2; // Linear Upwind Differencing
    hv[3:N+1] = (3*Gas[2:N].h - Gas[1:N-1].h)/2;
    for i in 1:nSpecies loop
      xiv[2, i] = (3*Gas[1].Xi[i] - xiv[1, i])/2;
      xiv[3:N+1, i] = (3*Gas[2:N].Xi[i] - Gas[1:N-1].Xi[i])/2;
    end for;
  else
    mfv[2:N+1] = mf[:]; // Upwind Differencing
    hv[2:N+1] = Gas[:].h;
    for i in 1:nSpecies loop
      xiv[2:N+1, i] = Gas[1:N].Xi[i];
    end for;
  end if;

  // Boundary conditions
  outfl.m_flow = -mfv[N+1];
  outfl.Xi_outflow[:] = xiv[N+1, :];
  outfl.h_outflow = hv[N+1];
  outfl.p = Gas[N].p; // Pressure is upwinded for the last control volume, alternative options given below
  //outfl.p = infl.p - dpdx*L ;
  //outfl.p = (1-dpFac)*infl.p ;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
          extent={{-80,34},{80,-32}},
          lineColor={28,108,200},
          fillColor={28,108,200},
          fillPattern=FillPattern.HorizontalCylinder)}),         Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end Flow1D;

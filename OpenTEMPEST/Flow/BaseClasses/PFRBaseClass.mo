within OpenTEMPEST.Flow.BaseClasses;
partial model PFRBaseClass
  "Base Class for a 1 phase Plug flow reactor with no external mass transfer"
  extends Flow1DBaseClass;
  import SI = Modelica.SIunits;

  // Vertex Values
  SI.MassFlowRate mfv[N+1] "Mass flow rate in and leaving CV";
  SI.SpecificEnthalpy hv[N+1] "Average enthalpy in CV";
  SI.MassFraction xiv[N+1, nSpecies] "Mass fractions in CV";

  // Catalyst parameters
  parameter Real por = 1 "Reactor porosity, leave as 1 if no heterogenous catalyst";
  parameter SI.Density rhoCat=1;
  parameter SI.SpecificHeatCapacity cpCat=1;

  // External Heat Transfer - total HT to surroundings
  Real Qext[N];

  // Total Species Generation Rate - User will need to define to enter the reaction rates at top level
  SI.MassFlowRate R[N, nSpecies]; // for single channel

//   SI.MassFlowRate D_mf[N+1] "Diffusive Mass Transport";
//   SI.MassFlowRate D_Xi[N+1, nSpecies] "Diffusive Species Transport";

equation

  // Mass balance
  por*dV.*theta[:,1].*der(Gas[:].d) = (mfv[1:N] .- mfv[2:N+1])./Ntubes;

  // Species Mass balance
  for i in 1:nSpecies loop
    por.*dV.*theta[:,1].*der(Gas[:].d.*Gas[:].Xi[i]) = (mfv[1:N].*xiv[1:N, i] .- mfv[2:N+1].*xiv[2:N+1, i])/Ntubes .+ R[:, i];
  end for;

  // Energy Balance - Heat generation from reaction already accounted for
  Ntubes*dV.*theta[:,1].*der(por.*Gas[:].d.*Gas[:].u .+ (1-por).*rhoCat.*cpCat.*Gas[:].T) = mfv[1:N].*hv[1:N] .- mfv[2:N+1].*hv[2:N+1] .+ Qext[:];

  // Momentum Balance Needed at top level

  // Interface interpolation
    // At x=0
  mfv[1] = infl.m_flow;
  hv[1] = infl.h_outflow;
  hv[1] = inStream(infl.h_outflow);
  xiv[1, :] = infl.Xi_outflow[:];
  xiv[1, :] = inStream(infl.Xi_outflow[:]);

//   D_mf[1] = 0;
//   D_mf[2:N] = 1e-3*(Gas[2:N].d .- Gas[1:N-1].d)*por*Acs/dx;
//   D_mf[N+1] = 0;
//
//   D_Xi[1,:] = zeros(nSpecies);
//   for j in 2:N loop
//     D_Xi[j, :] = 1e-3*(Gas[j].d.*Gas[j].Xi[1:nSpecies] .- Gas[j-1].d.*Gas[j-1].Xi[1:nSpecies])*por*Acs/dx;
//   end for;
//   D_Xi[N+1, :] = zeros(nSpecies);

    // First CV (x=1*dx)
//   mfv[2] = (3*mf[1] - mfv[1])/2;//mf[1]; //homotopy((3*mf[1] - mfv[1])/2, mf[1]);
//   hv[2] = (3*Gas[1].h - hv[1])/2; //homotopy((3*Gas[1].h - hv[1])/2, Gas[1].h);
//   xiv[2, 1:nSpecies] = (3*Gas[1].Xi[1:nSpecies] - xiv[1, 1:nSpecies])/2; //homotopy((3*Gas[1].Xi[1:nSpecies] - xiv[1, 1:nSpecies])/2, Gas[1].Xi[1:nSpecies]);
//   for j in 2:N loop // Remaining CVs (x=j*dx)
//     mfv[j+1] = (3*mf[j] - mf[j-1])/2;//mf[j];// + 1e-2*(Gas[j+1].d - Gas[j].d);  //homotopy((3*mf[j] - mf[j-1])/2, mf[j]);
//     hv[j+1] = (3*Gas[j].h - Gas[j-1].h)/2; //homotopy((3*Gas[j].h - Gas[j-1].h)/2, Gas[j].h);
//     xiv[j+1, 1:nSpecies] = (3*Gas[j].Xi[1:nSpecies] - Gas[j-1].Xi[1:nSpecies])/2; // .+ 1e-3*(Gas[j+1].Xi[1:nSpecies] .- Gas[j].Xi[1:nSpecies]);   //homotopy((3*Gas[j].Xi[1:nSpecies] - Gas[j-1].Xi[1:nSpecies])/2, Gas[j].Xi[1:nSpecies]);
//   end for;
  mfv[2] = (3*mf[1] - mfv[1])/2;
  hv[2] = (3*Gas[1].h - hv[1])/2;
  // Remaining CVs
  mfv[3:N+1] = (3*mf[2:N] - mf[1:N-1])/2; // Linear Upwind Differencing
  hv[3:N+1] = (3*Gas[2:N].h - Gas[1:N-1].h)/2;
  for i in 1:nSpecies loop
    xiv[2, i] = (3*Gas[1].Xi[i] - xiv[1, i])/2;
    xiv[3:N+1, i] = (3*Gas[2:N].Xi[i] - Gas[1:N-1].Xi[i])/2;
  end for;
  //mfv[N+1] = (3*mf[N] - mf[N-1])/2;//mf[N];
  //hv[N+1] = Gas[N].h;
  //xiv[N+1, 1:nSpecies] = Gas[N].Xi[1:nSpecies];

  // Boundary conditions
  outfl.m_flow = -mfv[N+1];
  outfl.Xi_outflow[:] = xiv[N+1, :];
  outfl.h_outflow = hv[N+1];
  outfl.p = Gas[N].p; // Pressure is upwinded for the last control volume, alternative options given below

end PFRBaseClass;

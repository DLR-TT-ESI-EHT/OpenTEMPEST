within OpenTEMPEST.BOP.Reactors;
model Reformer1D
  extends Flow.BaseClasses.PFRBaseClass(redeclare package Medium =
        OpenTEMPEST.Medium.Fuel_CH4);
import SI = Modelica.SIunits;

  parameter Real Fac_ref = 1;
  parameter Real Fac_wgs = 0.75;
  // Kinetic parameters
  parameter Real A_sf = 0.0171;
  parameter Real Ea_sf = 103191*Fac_wgs;
  parameter Real A_rf = 2395;
  parameter Real Ea_rf = 231266*Fac_ref;

  // Kinetics
  Real a[nSpecies,2]={{3,1},{-1,0},{0,1},{1,-1},{-1,-1},{0,0}};
  Real krf[N];
  Real ksf[N];
  Real r1[N];
  Real r2[N];

  Real Kref[N];
  Real Kwgs[N]
    "Reaction equillibrium constant of reverse reforming and WGS reaction respectively";
  Real DeltaG_wgs[N];
  Real DeltaG_ref[N];
  Real Yi[N, Medium.nXi];

  Real Qref[N];
  Real Eq_ref[N];

  parameter Real dpdx=100 " Constant pressure gradient Pa/m";
  parameter SI.CoefficientOfHeatTransfer Uconv=100 "Overall HTC /W/m2/K";

  ThermoPower.Thermal.DHTVolumes wall(N=N)
                                       "HT through cylindrical wall" annotation (Placement(
        transformation(extent={{-80,30},{80,50}}),iconTransformation(extent={{-80,34},
            {80,50}})));

equation

// Reaction kinetics
  for i in 1:N loop
    r1[i] = krf[i]*(Gas[i].p*Yi[i, 2]*Gas[i].p*Yi[i, 5] - (Gas[i].p*Yi[i, 1])^3*(Gas[i].p*Yi[i,4])/Kref[i]); //mol/m3.s
    r2[i] = ksf[i]*(Gas[i].p*Yi[i, 4]*Gas[i].p*Yi[i, 5] - Gas[i].p*Yi[i, 1]*Gas[i].p*Yi[i, 3]/Kwgs[i]); //mol/m3.s
    Qref[i] = ((Gas[i].p*Yi[i, 1])^3*(Gas[i].p*Yi[i,4]))/(Gas[i].p*Yi[i, 2]*Gas[i].p*Yi[i, 5]);
    Eq_ref[i] = Qref[i]/Kref[i];

    DeltaG_ref[i] = -252.642810968035.*Gas[i].T + 225215.698063031;
    Kref[i] = 1e10*exp(-DeltaG_ref[i]/Modelica.Constants.R/Gas[i].T);
    krf[i] = A_rf*Modelica.Math.exp(-Ea_rf/(8.314*Gas[i].T));

    DeltaG_wgs[i] = 32.1153*(Gas[i].T) - 3.5211E4; // Marius' shortcut from NASA
    Kwgs[i] = exp(-DeltaG_wgs[i]/Modelica.Constants.R/Gas[i].T);
    ksf[i] = A_sf*Modelica.Math.exp(-Ea_sf/(8.314*Gas[i].T));

    Yi[i,:] = Gas[i].Xi[:]./Medium.MMX[:]/sum(Gas[i].Xi[:]./Medium.MMX[:]); // this line is able to be compiled by openmodelica (1.14), enabling FMU use.

    for j in 1:nSpecies loop
      R[i, j] =  Medium.MMX[j]*(dV*por*a[j, 1]*r1[i] .+ dV*por*a[j, 2]*r2[i]);
    end for;
  end for;

// Momentum Balance
  0 = (infl.p - Gas[1].p) - dpdx*dx;
  fill(0, N-1) = (Gas[1:N-1].p .- Gas[2:N].p) .- dpdx*dx;

  Qext = wall.Q;
  wall.Q = -Uconv*Aw*(Gas[:].T .- wall.T[:]);

end Reformer1D;

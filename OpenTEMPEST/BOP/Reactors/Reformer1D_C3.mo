within OpenTEMPEST.BOP.Reactors;
model Reformer1D_C3

  extends Flow.BaseClasses.PFRBaseClass(
    cpCat=710,
    rhoCat=7130,
    por=0.75,
    redeclare package Medium = OpenTEMPEST.Medium.Fuel_C3);
  import SI = Modelica.SIunits;

  ThermoPower.Thermal.DHTVolumes wall(N=N)
                                       "HT through cylindrical wall" annotation (Placement(
        transformation(extent={{-80,30},{80,50}}),iconTransformation(extent={{-80,34},
            {80,50}})));

  SI.TemperatureSlope dTcond[N] "left and right temperature gradient for conductive heat flux";
  parameter SI.ThermalConductivity lambda=22;

  Real r1[N];
  Real r2[N];
  Real r3[N];

  constant Integer a1[nSpecies]={1,0,1,-1,-1,0,0};
  constant Integer a2[nSpecies]={3,-1,0,1,-1,0,0};
  constant Integer a3[nSpecies]={10,0,3,0,-6,0,-1};

  parameter Real A1=174;
  parameter Real A2=239.5;
  parameter Real A3=50;
  parameter Real Ea1=103191;
  parameter Real Ea2=231266*0.82;
  parameter Real Ea3=159241;
  parameter Real alpha3=2.8;
  parameter Real beta3=-1;

  parameter Real UA=0.005;

  Real K_WGS[N], K_MSR[N], Q_MSR[N], Q_WGS[N];
  Real Eq_MSR, Eq_WGS;
  Real Xtot[nSpecies];

  parameter Real dpdxNom=1;

  SI.MoleFraction Yi[N, nSpecies];
  SI.MoleFraction Ydry[nSpecies];

equation

  // Reaction rate
  for j in 1:nSpecies loop
    R[:, j] = dV*por*Medium.MMX[j]*(a1[j]*r1
               .+ a2[j]*r2
               .+ a3[j]*r3);
  end for;

  for i in 1:N loop
    Yi[i,:] = Medium.massToMoleFractions(Gas[i].Xi[:], Medium.MMX);
    r1[i] = A1*Modelica.Math.exp(-Ea1/Modelica.Constants.R/Gas[i].T)*(Yi[i,4]*Yi[i,5]*Gas[i].p^2 - (1/K_WGS[i])*Yi[i,3]*Yi[i,1]*(Gas[i].p^2));
    r2[i] = A2*Modelica.Math.exp(-Ea2/Modelica.Constants.R/Gas[i].T)*(Yi[i,2]*Yi[i,5]*Gas[i].p^2 - (1/K_MSR[i])*Yi[i,4]*(Yi[i,1]^3)*(Gas[i].p^4));
    r3[i] = A3*abs(Gas[i].T/1000)^beta3*Modelica.Math.exp(-Ea3/Modelica.Constants.R/Gas[i].T)*abs(Yi[i,7]*Gas[i].p)^alpha3;
    K_MSR[i] = 101325^2*Modelica.Math.exp(-(-252.642810968035*Gas[i].T + 225215.698063031)/Modelica.Constants.R/Gas[i].T);
    K_WGS[i] = Modelica.Math.exp(-(32.1153*(Gas[i].T) - 3.5211E4)/Modelica.Constants.R/Gas[i].T);
    Q_MSR[i] = Yi[i,4]*(Yi[i,1]^3)*(Gas[i].p^2)/(Yi[i,2]*Yi[i,5]);
    Q_WGS[i] = Yi[i,3]*Yi[i,1]/(Yi[i,4]*Yi[i,5]);
  end for;
  Eq_MSR = Q_MSR[N]/K_MSR[N];
  Eq_WGS = Q_WGS[N]/K_WGS[N];

  Ydry[5] = 0;
  Ydry[1] = Yi[N,1]./(Yi[N,1] .+ Yi[N,2] .+ Yi[N,3] .+ Yi[N,4] .+ Yi[N,6] .+ Yi[N,7]);
  Ydry[2] = Yi[N,2]./(Yi[N,1] .+ Yi[N,2] .+ Yi[N,3] .+ Yi[N,4] .+ Yi[N,6] .+ Yi[N,7]);
  Ydry[3] = Yi[N,3]./(Yi[N,1] .+ Yi[N,2] .+ Yi[N,3] .+ Yi[N,4] .+ Yi[N,6] .+ Yi[N,7]);
  Ydry[4] = Yi[N,4]./(Yi[N,1] .+ Yi[N,2] .+ Yi[N,3] .+ Yi[N,4] .+ Yi[N,6] .+ Yi[N,7]);
  Ydry[6] = Yi[N,6]./(Yi[N,1] .+ Yi[N,2] .+ Yi[N,3] .+ Yi[N,4] .+ Yi[N,6] .+ Yi[N,7]);
  Ydry[7] = Yi[N,7]./(Yi[N,1] .+ Yi[N,2] .+ Yi[N,3] .+ Yi[N,4] .+ Yi[N,6] .+ Yi[N,7]);

  for j in 1:nSpecies loop
    Xtot[j] = (mfv[N]*xiv[N,j] - mfv[1]*xiv[1,j])/(mfv[1]*xiv[1,j] + Modelica.Constants.eps);
  end for;

  // Heat transfer
  Qext= wall.Q .+ lambda*((1-por)*Acs)*dTcond;
    // wall to reactor
  wall.Q = UA*(wall.T.-Gas.T);
    // Internal conductivity
  dTcond[1] = (Gas[2].T .- Gas[1].T)/dx;
  dTcond[2:N-1] = (Gas[1:N-2].T .- 2*Gas[2:N-1].T + Gas[3:N].T)/dx;
  dTcond[N] = (Gas[N-1].T .- Gas[N].T)/dx;

  // Momentum Balance
  0 = (infl.p - Gas[1].p) - dpdxNom*dx;
  fill(0, N-1) = (Gas[1:N-1].p .- Gas[2:N].p) .- dpdxNom*dx;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end Reformer1D_C3;

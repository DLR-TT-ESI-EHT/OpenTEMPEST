within OpenTEMPEST.BOP.Reactors;
model Combustor0D
  "0D combustor assuming complete combustion of fuels. Needs lambda >=1"

  import SI = Modelica.SIunits;
  package Air = Medium.Air_Medium;
  package Fuel = Medium.Fuel_CH4;
  package Exhaust = Medium.SimpleExhaust;

  parameter SI.Volume V "Inner volume";
  parameter SI.Area S=0 "Inner surface";
  parameter SI.CoefficientOfHeatTransfer gamma=0  "Heat Transfer Coefficient"    annotation (Evaluate=true);

  parameter SI.HeatCapacity Cm=0 "Metal Heat Capacity" annotation (Evaluate=true);
  parameter SI.Temperature Tmstart=300 "Metal wall start temperature"  annotation (Dialog(tab="Initialisation"));

  parameter Air.AbsolutePressure pstart=101325 "Pressure start value"    annotation (Dialog(tab="Initialisation"));
  parameter Air.Temperature Tstart=300 "Temperature start value"    annotation (Dialog(tab="Initialisation"));
  parameter Air.MassFraction Xstart[Exhaust.nX]=Exhaust.reference_X   "Start flue gas composition" annotation (Dialog(tab="Initialisation"));

  Exhaust.BaseProperties fluegas(
    p(start=pstart),
    T(start=Tstart),
    Xi(start=Xstart[1:Exhaust.nXi]));
  SI.Mass M "Gas total mass";
  SI.Mass MX[Exhaust.nXi] "Partial flue gas masses";
  SI.InternalEnergy E "Gas total energy";
  SI.Temperature Tm(start=Tmstart) "Wall temperature";
  Air.SpecificEnthalpy hia "Air specific enthalpy";
  Fuel.SpecificEnthalpy hif "Fuel specific enthalpy";

  SI.PerUnit lambda    "Stoichiometric ratio (>1 if air flow is greater than stoichiometric)";

  SI.Time Tr "Residence time";
  ThermoPower.Gas.FlangeA ina(redeclare package Medium = Air) "inlet air"
    annotation (Placement(transformation(extent={{-120,-20},{-80,20}}, rotation=
           0)));
  ThermoPower.Gas.FlangeA inf(redeclare package Medium = Fuel) "inlet fuel"
    annotation (Placement(transformation(extent={{-20,80},{20,120}}, rotation=0)));
  ThermoPower.Gas.FlangeB out(redeclare package Medium = Exhaust) "flue gas"
    annotation (Placement(transformation(extent={{80,-20},{120,20}}, rotation=0)));

  constant Real[3, Exhaust.nXi] a_i= {
          {-0.5, 0, 1, 0},
          {-2,   1, 2, 0},
          {-0.5, 1, 0, 0}};

  Real[3] R = {
        inf_X[1]/Fuel.data[1].MM,
        inf_X[2]/Fuel.data[2].MM,
        inf_X[4]/Fuel.data[4].MM}
   *inf.m_flow;

  Air.MassFraction ina_X[Air.nXi]=inStream(ina.Xi_outflow);
  Fuel.MassFraction inf_X[Fuel.nXi]=inStream(inf.Xi_outflow);

equation
  M = fluegas.d*V "Gas mass";
  E = fluegas.u*M "Gas energy";
  MX = fluegas.Xi*M "Component masses";

  der(M) = ina.m_flow + inf.m_flow + out.m_flow "Gas mass balance";
  der(E) = ina.m_flow*hia + inf.m_flow*hif + out.m_flow*fluegas.h - gamma*S*(fluegas.T - Tm) "Gas energy balance";
  if Cm > 0 and gamma > 0 then
    Cm*der(Tm) = gamma*S*(fluegas.T - Tm) "Metal wall energy balance";
  else
    Tm = fluegas.T;
  end if;

  if inf.m_flow <=1e-8 then
     lambda = 100;
  else
    lambda  = (ina.m_flow*ina_X[1]/Air.data[1].MM)/(-a_i[:,1]*R);
  end if;
  //assert(lambda >= 1, "Not enough oxygen flow");

  der(MX[1]) = ina.m_flow*ina_X[1] + out.m_flow*fluegas.X[1] + (a_i[:,1]*R[:])*Exhaust.data[1].MM "oxygen";
  der(MX[2]) = inf.m_flow*inf_X[3] + out.m_flow*fluegas.X[2] + (a_i[:,2]*R[:])*Exhaust.data[2].MM "carbondioxide";
  der(MX[3]) = inf.m_flow*inf_X[5] + out.m_flow*fluegas.X[3] + (a_i[:,3]*R[:])*Exhaust.data[3].MM "water";
  der(MX[4]) = ina.m_flow*ina_X[2] + out.m_flow*fluegas.X[4] + inf.m_flow*inf_X[6] "nitrogen";

  // Set gas properties
  out.p = fluegas.p;
  out.h_outflow = fluegas.h;
  out.Xi_outflow = fluegas.Xi;

  // Boundary conditions
  ina.p = fluegas.p;
  ina.h_outflow = 0;
  ina.Xi_outflow = Air.reference_X[1:Air.nXi];
  inf.p = fluegas.p;
  inf.h_outflow = 0;
  inf.Xi_outflow = Fuel.reference_X[1:Fuel.nXi];
  hia = inStream(ina.h_outflow);
  hif = inStream(inf.h_outflow);

  Tr = noEvent(M/max(abs(out.m_flow), Modelica.Constants.eps));

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Ellipse(
          extent={{-90,82},{90,-100}},
          lineColor={28,108,200},
          fillColor={28,108,200},
          fillPattern=FillPattern.Solid)}),                      Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end Combustor0D;

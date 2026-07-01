within OpenTEMPEST.Flow.BaseClasses;
partial model Flow1DBaseClass
  "Base class model for 1D flow model based on FVM or FDM. Can be used to build PFR models or heterogeneous catalytic reactor models"

  import SI = Modelica.SIunits;
  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
    annotation(choicesAllMatching = true);
  parameter Integer N(min=3) = 5;
  parameter Integer Ntubes = 1 "number of individual tubes";

  // Initial Values
  parameter SI.Temperature TStartIn=773.15   annotation (Dialog(tab="Initialisation"));
  parameter SI.Temperature TStartOut=773.15   annotation (Dialog(tab="Initialisation"));
  parameter SI.Temperature TStart[N] = linspace(TStartIn, TStartOut, N) annotation (Dialog(tab="Initialisation"));
  parameter SI.AbsolutePressure pStart=101325   annotation (Dialog(tab="Initialisation"));
  parameter SI.MassFraction xStart[nSpecies] = Medium.reference_X annotation (Dialog(tab="Initialisation"));

  // Dimensions
  parameter SI.Length Rh = 2.5e-2 "Hydraulic Radius single channel (Dh/2) /m" annotation (Dialog(tab="Dimensions"));
  parameter SI.Length omega = 2*Modelica.Constants.pi*Rh "Wetted Perimeter single channel /m" annotation (Dialog(tab="Dimensions"));
  parameter Real theta0_1 = 1 "Proportion of entrance and exit volumes to centre volumes /-" annotation (Dialog(tab="Dimensions"));
  parameter SI.Length L = 1 "Length of pipe /m" annotation (Dialog(tab="Dimensions"));

  parameter SI.Area Acs = Rh*omega/2 "channel cross sectional area" annotation (Dialog(tab="Dimensions"));
  parameter SI.Area Aw = omega*dx "HT area for one CV single channel" annotation (Dialog(tab="Dimensions"));
  parameter Real theta[N,1] = [theta0_1; fill(1, N-2); theta0_1] "coefficient for CV size" annotation (Dialog(tab="Dimensions"));

  parameter SI.Volume dV = Acs*dx "Size of nominal control volume /m3" annotation (Dialog(tab="Dimensions"));
  parameter SI.Length dx = L/(N - 2 + 2*theta0_1) "length of nominal 1D control volume /m" annotation (Dialog(tab="Dimensions"));

  Medium.BaseProperties Gas[N](
    p(each start=pStart, each stateSelect=StateSelect.prefer),
    T(start=TStart[1:N], each stateSelect=StateSelect.prefer),
    Xi(each stateSelect=StateSelect.prefer))
    "Gas volume properties";

  // Centre cell values
  SI.MassFlowRate mf[N] "Mass flow rate in and leaving CV";

  ThermoPower.Gas.FlangeA infl(redeclare package Medium = Medium)
      annotation (Placement(transformation(extent={{-100,-8},{-80,12}}),
        iconTransformation(extent={{-108,-16},{-80,12}})));
  ThermoPower.Gas.FlangeB outfl(redeclare package Medium = Medium)
      annotation (Placement(transformation(extent={{80,-10},{100,10}}),
        iconTransformation(extent={{80,-14},{108,14}})));

protected
  constant Integer nSpecies = Medium.nXi;

initial equation
  //Gas[:].T = TStart[:];
  Gas[:].Xi = fill(xStart, N);

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
          extent={{-80,34},{80,-32}},
          lineColor={28,108,200},
          fillColor={28,108,200},
          fillPattern=FillPattern.HorizontalCylinder)}),         Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end Flow1DBaseClass;

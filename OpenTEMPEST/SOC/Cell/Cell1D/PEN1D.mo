within OpenTEMPEST.SOC.Cell.Cell1D;
model PEN1D "1D FV PEN model for extending Solid1D"
  extends Heat.BaseClasses.Solid1DBase;

  import SI = Modelica.SIunits;

  parameter SI.CurrentDensity Jstart = 0;

  constant SI.AbsolutePressure p0 = 1e5 "Standard Pressure";

  replaceable model Electrochem =
      OpenTEMPEST.SOC.Electrochem.Components.BV_Steam          constrainedby
    OpenTEMPEST.SOC.Electrochem.Components.ElectrochemBase
      annotation (Placement(transformation(extent={{78,50},{98,70}})), choicesAllMatching=true);

  Electrochem electrochem[N](Tpen=T, J=J, P_A=PEN_ina.P, P_F=PEN_in.P, p0=fill(p0,N), yA=PEN_ina.Y, yF=PEN_in.Y);

  SI.Voltage Uop[N] "Operating Voltage of Cell";
  SI.CurrentDensity J[N](each start=Jstart) "Current Density of Control Volume";
  SI.Current Icell "Total current into Cell";
  SI.Current Iv[N];
  SI.EnergyFlowRate q_electrochem[N] = electrochem.q_electroChem*dx*lY;
  SI.MolarFlowRate rEl[N] = electrochem.r*dx*lY;

  OpenTEMPEST.SOC.Electrochem.Interfaces.VariablesStream PEN_ina[N](each
      nspecies=Medium.Air_Medium.nXi) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=-90,
        origin={-14,-30}), iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=-90,
        origin={0,-70})));
  OpenTEMPEST.SOC.Electrochem.Interfaces.VariablesStream PEN_in[N](each
      nspecies=Medium.Fuel_CH4.nXi) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={-8,36}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={0,70})));
  Modelica.Electrical.Analog.Interfaces.NegativePin pin_n annotation (
    Placement(transformation(extent={{-104,20},{-84,40}}),      iconTransformation(extent={{-106,34},{-96,44}})));
  Modelica.Electrical.Analog.Interfaces.PositivePin pin_p annotation (
    Placement(transformation(extent={{-116,-56},{-96,-36}}),      iconTransformation(extent={{-106,-46},{-96,-36}})));

  ThermoPower.Thermal.DHTVolumes Qrad_FI(N=N)
    annotation (Placement(transformation(extent={{44,10},{64,30}}), iconTransformation(extent={{44,10},{64,30}})));
  ThermoPower.Thermal.DHTVolumes Qrad_AI(N=N)
    annotation (Placement(transformation(extent={{46,-30},{66,-10}}), iconTransformation(extent={{46,-30},{66,-10}})));

equation
  // Thermal
  Qext = Qrad_FI.Q .+ Qrad_AI.Q .- dx*lY*(Uop.*J) .+ q_electrochem;
  Qrad_FI.T = T; // Radiative Flow
  Qrad_AI.T = T; // Radiative Flow

  // Electrochemistry
  electrochem.Uop   = Uop;

  // Electrical connections
  Uop = fill(pin_p.v - pin_n.v, N);
  0 = pin_p.i + pin_n.i;
  Icell = pin_p.i;
  J*dx*lY = Iv;
  sum(Iv) = Icell;

  // Mass Transfer
  PEN_in.I_H = electrochem.J_H*dx*lY;
  PEN_in.I_C = electrochem.J_C*dx*lY;
  PEN_ina.I_H = zeros(N);
  PEN_ina.I_C = zeros(N);
  PEN_in.I = 2*Modelica.Constants.F*rEl;  // in metal supported cell J is external circuit only and not whole reaction rate
  PEN_ina.I = 2*Modelica.Constants.F*rEl;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end PEN1D;

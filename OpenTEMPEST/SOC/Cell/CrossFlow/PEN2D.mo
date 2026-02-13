within OpenTEMPEST.SOC.Cell.CrossFlow;
model PEN2D
  extends Heat.BaseClasses.Solid2DBase(        redeclare package SolidMat =
        TEMPEST.Solid.Material.Custom,
        cpCustom=500,
    rhoCustom(displayUnit="kg/m3") = 5900,
    kCustom_trans=2,
    lX=9.0e-2,
    lY=1.42e-1,
    lZ=1.51e-4);

  parameter Modelica.SIunits.CurrentDensity Jstart=0;
  constant Modelica.SIunits.AbsolutePressure p0=1e5 "Standard Pressure";

  replaceable model Electrochem =
       OpenTEMPEST.SOC.Electrochem.Components.Crossflow_Electrochem     constrainedby
    OpenTEMPEST.SOC.Electrochem.Components.ElectrochemBase
       annotation (Placement(transformation(extent={{78,50},{98,70}})), choicesAllMatching=true);

   Electrochem electrochem[nX, nY](Tpen=T, J=J, P_A=PEN_ina.P, P_F=PEN_in.P, p0=fill(p0, nX, nY), yA=PEN_ina.Y, yF=PEN_in.Y);

   Modelica.SIunits.Voltage Uop[nX,nY] "Operating Voltage of Cell";
   Modelica.SIunits.CurrentDensity J[nX,nY](each start=Jstart) "Current Density of Control Volume";
   Modelica.SIunits.Current Icell "Total current into Cell";
   Modelica.SIunits.Current Iv[nX, nY];
   Modelica.SIunits.EnergyFlowRate q_electrochem[nX, nY] = electrochem.q_electroChem*dx*lY;
   Modelica.SIunits.MolarFlowRate rEl[nX, nY] = electrochem.r*dx*lY;

  OpenTEMPEST.SOC.Electrochem.Interfaces.VariablesStream PEN_ina[nX,nY](each
      nspecies=Medium.Air_Medium.nXi) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=-90,
        origin={0,-70}), iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=-90,
        origin={0,-70})));
  OpenTEMPEST.SOC.Electrochem.Interfaces.VariablesStream PEN_in[nX,nY](each
      nspecies=Medium.Fuel_CH4.nXi) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={0,70})));
  Modelica.Electrical.Analog.Interfaces.NegativePin pin_n annotation (
    Placement(transformation(extent={{-104,20},{-84,40}}),      iconTransformation(extent={{-106,30},{-96,40}})));
  Modelica.Electrical.Analog.Interfaces.PositivePin pin_p annotation (
    Placement(transformation(extent={{-116,-54},{-96,-34}}),      iconTransformation(extent={{-106,-44},{-96,-34}})));

  Heat.DHTVolumes2D Qrad_FI(i=nX, j=nY) annotation (Placement(transformation(
          extent={{40,12},{60,32}}), iconTransformation(extent={{40,12},{60,32}})));
  Heat.DHTVolumes2D Qrad_AI(i=nX, j=nY) annotation (Placement(transformation(
          extent={{40,-28},{60,-8}}), iconTransformation(extent={{40,-28},{60,-8}})));

equation

  // Electrochemistry
   electrochem.Uop   = Uop;

  // Thermal
  Qext[:,:] = Qrad_FI.Q .+ Qrad_AI.Q .+ q_electrochem .- dx*dy*(Uop.*J);
  Qrad_FI.T[:,:] = T; // Radiative Flow
  Qrad_AI.T[:,:] = T; // Radiative Flow

  // Electrical connections
  Uop[:, :] = fill(pin_p.v - pin_n.v, nX, nY);
  0 = pin_p.i + pin_n.i;
  Icell = pin_p.i;
  J*dx*dy = Iv;
  sum(Iv[:,:]) = Icell;

  // Mass Transfer
  PEN_in.I = 2*Modelica.Constants.F*rEl; // PEN_in.I_h = i_h*L*B
  PEN_ina.I = 2*Modelica.Constants.F*rEl; //PEN_ina.I_h = 0
  PEN_in.I_H = electrochem.J_H*dx*dy;
  PEN_in.I_C = electrochem.J_C*dx*dy;
  PEN_ina.I_H = zeros(nX,nY);
  PEN_ina.I_C = zeros(nX,nY);

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end PEN2D;

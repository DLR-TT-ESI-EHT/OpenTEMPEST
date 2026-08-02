within OpenTEMPEST.SOC.Cell.CrossFlow;
model CellGeneral "Crossflow generalized cell model with possibility to select 2D or 1D channels"

  import SI = Modelica.SIunits;

  parameter Integer nX=5 "Number of control volumes in the x-direction";
  parameter Integer nY=5 "Number of control volumes in the y-direction";
  parameter Boolean simplification=false  "Choice between simplified 1D channels or 2D channels" annotation(Dialog(group="Selections"), choices(__Dymola_radioButtons=true, choice=true "1D Channels", choice=false "2D Channels"));
  parameter Boolean LUDS=false "Linear Upwind difference (more accuracy), Upwind difference scheme (more speed and stability)" annotation(Dialog(group="Selections"), choices(__Dymola_radioButtons=true, choice=true "LUDS", choice=false "UDS"));
  parameter Boolean heatTransferCorrelationFormDuct=true "Nusselt correlation duct geometry with characteristic length=2*lZ (default), plate geometry with characteristic length=lX" annotation(Dialog(group="Selections"), choices(__Dymola_radioButtons=true, choice=true "Duct", choice=false "Plate"));

  // Dimensions from DOI 10.1149/05701.0089ecst
  parameter SI.Length lX = 0.115 "Length of active area" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lY = 0.110 "Width of active area" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lXpen = lX "Length of pen" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lYpen = lY "Width of pen" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lZpen = 1.65e-4 "Thickness of pen" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lXac = lX "Length of air channel" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lYac = lY "Width of air channel" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lZac = 0.53e-3 "Height of air channel" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lXfc = lX "Length of fuel channel" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lYfc = lY "Width of fuel channel" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lZfc = 0.4e-3 "Height of fuel channel" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lZsolid = 1.6e-3 "Height of interconnector"  annotation(Dialog(tab="Dimensions"));
  parameter SI.Length Bac = porAC*lYac "Width of air channel without ribs" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length Bfc = porFC*lYfc "Width of fuel channel without ribs" annotation(Dialog(tab="Dimensions"));

  // Initialization
  parameter SI.Temperature TStart=1023.15 "Uniform start temperature" annotation(Dialog(tab="Initialization"));
  parameter SI.AbsolutePressure pStart=101325 "Starting pressure" annotation(Dialog(tab="Initialization"));
  parameter SI.CurrentDensity Jstart = 0 "Starting current density in PEN" annotation(Dialog(tab="Initialization"));
  parameter SI.MassFraction xStartFC[FCMedium.nXi] = FCMedium.reference_X "Starting mass fraction in fuel channel" annotation(Dialog(tab="Initialization"));
  parameter SI.MassFraction xStartAC[ACMedium.nX] = ACMedium.reference_X "Starting mass fraction in air channel" annotation(Dialog(tab="Initialization"));

  // PEN parameters - 10Sc1CeSZ
  replaceable model Electrochem =
      OpenTEMPEST.SOC.Electrochem.Components.Crossflow_Electrochem
      constrainedby OpenTEMPEST.SOC.Electrochem.Components.ElectrochemBase annotation (
  Dialog(tab="PEN"),
  Placement(transformation(extent={{78,50},{98,70}})),
  choicesAllMatching=true);
  parameter SI.ThermalConductivity kCustom_trans = 2.3 "Thermal Conductivity across layers of PEN in W/mK" annotation(Dialog(tab="PEN"));
  parameter SI.ThermalConductivity kCustom_long = 2.3 "Thermal Conductivity in plane of layers of PEN in W/mK (=k_trans for homogeneous materials)" annotation(Dialog(tab="PEN"));
  parameter SI.Density rhoPEN=5740 "Density of PEN in kg/m3 from https://doi.org/10.1016/j.matpr.2018.10.071" annotation(Dialog(tab="PEN"));
  parameter SI.SpecificHeatCapacity cpPEN = 500 "Specific heat capacity of PEN in J/kgK" annotation(Dialog(tab="PEN"));
  parameter SI.SpectralEmissivity epsilonPEN = 0.8 "Emissivity of Anode-Electrolyte-Cathode unit" annotation(Dialog(tab="PEN"));

  // FC parameters
  replaceable package FCMedium = OpenTEMPEST.Medium.Fuel_CH4 annotation(Dialog(tab="Fuel Channel"));
  parameter Real porFC = 0.6 "Porosity in fuel channel" annotation(Dialog(tab="Fuel Channel"));
  parameter Real NuFCPEN = 12 "Nusselt number fuel channel on PEN side" annotation(Dialog(tab="Fuel Channel"));
  parameter Real NuFCIC = 9.86 "Nusselt number fuel channel on IC side" annotation(Dialog(tab="Fuel Channel"));
  parameter Real alfafc=0.01528182
                               "Weight for 2D temperature in z-direction convection for fuel channel" annotation(Dialog(tab="Fuel Channel"));

  // AC parameters
  replaceable package ACMedium = OpenTEMPEST.Medium.Air_Medium annotation(Dialog(tab="Air Channel"));
  parameter Real porAC = 0.5 "Porosity in air channel" annotation(Dialog(tab="Air Channel"));
  parameter Real pDrop(max=0.99) = 0.04 "Pressure drop as a factor of inlet pressure (between 0 and 0.99)" annotation(Dialog(tab="Air Channel"));
  parameter Real NuACPEN = 8.235 "Nusselt number air channel on PEN side" annotation(Dialog(tab="Air Channel"));
  parameter Real NuACIC = 7.54 "Nusselt number air channel on IC side" annotation(Dialog(tab="Air Channel"));
  parameter Real alfaac=0.11955708
                               "Weight for 2D temperature in z-direction convection for air channel" annotation(Dialog(tab="Air Channel"));

  // IC parameters - CFY (94.9% Cr, 5% Fe, 0.1% Y)
  parameter SI.ThermalConductivity kIC = 40 "Thermal conductivity of CFY interconnects in W/mK (35-45 W/mK for 20-900 °C)" annotation(Dialog(tab="Interconnects"));
  parameter SI.SpecificHeatCapacity cpIC = 451.8 "Ribs heat capacity in W/mK" annotation(Dialog(tab="Interconnects"));
  parameter SI.Density rhoIC=7233 "Ribs density in kg/m3" annotation(Dialog(tab="Interconnects"));
  parameter SI.SpectralEmissivity epsilonIC = 0.1 "Emissivity of interconnects" annotation(Dialog(tab="Interconnects"));

  // PEN
  Components.PEN2D pen(
    nX=nX,
    nY=nY,
    Tstartbar=TStart,
    lX=lXpen,
    lY=lYpen,
    lZ=lZpen,
    kCustom_trans=kCustom_trans,
    kCustom_long=kCustom_long,
    rhoCustom=rhoPEN,
    cpCustom=cpPEN,
    Jstart=Jstart,
    redeclare model Electrochem = Electrochem) annotation (Placement(transformation(extent={{-34,-20},{34,44}})));

  // Simplified 1D Air Channel
  Components.AirChannel2DSimp airChannel1D(
    nX=nY,
    nY=nX,
    TStart=TStart,
    pStart=pStart,
    xStart=xStartAC,
    lX=lYac,
    lY=lXac,
    lZ=lZac,
    por=porAC,
    pDrop=pDrop,
    Nu_PEN=NuACPEN,
    Nu_IC=NuACIC,
    alfa=alfaac,
    heatTransferCorrelationFormDuct=heatTransferCorrelationFormDuct) if simplification;

  // 2D Air Channel
  Components.AirChannel2D airChannel2D(
    nX=nY,
    nY=nX,
    TStart=TStart,
    pStart=pStart,
    xStart=xStartAC,
    lX=lYac,
    lY=lXac,
    lZ=lZac,
    por=porAC,
    pDrop=pDrop,
    Nu_PEN=NuACPEN,
    Nu_IC=NuACIC,
    heatTransferCorrelationFormDuct=heatTransferCorrelationFormDuct) if not simplification
    annotation (Placement(transformation(extent={{-38,-24},{44,-92}})));

  // Simplified 1D Fuel Channel
  Components.FuelChannel2DSimp fuelChannel1D(
    nX=nX,
    nY=nY,
    TStart=TStart,
    pStart=pStart,
    xStart=xStartFC,
    lX=lXfc,
    lY=lYfc,
    lZ=lZfc,
    por=porFC,
    kRibs=kIC,
    cpRibs=cpIC,
    rhoRibs=rhoIC,
    Nu_PEN=NuFCPEN,
    Nu_IC=NuFCIC,
    alfa=alfafc,
    heatTransferCorrelationFormDuct=heatTransferCorrelationFormDuct) if simplification;

  // 2D Fuel Channel
  Components.FuelChannel2D fuelChannel2D(
    nX=nX,
    nY=nY,
    TStart=TStart,
    pStart=pStart,
    xStart=xStartFC,
    lX=lXfc,
    lY=lYfc,
    lZ=lZfc,
    por=porFC,
    kRibs=kIC,
    cpRibs=cpIC,
    rhoRibs=rhoIC,
    Nu_PEN=NuFCPEN,
    Nu_IC=NuFCIC,
    heatTransferCorrelationFormDuct=heatTransferCorrelationFormDuct) if not simplification annotation (Placement(transformation(extent={{-34,32},{40,86}})));

  // Interconnects
  Heat.Solid2D iCAir(
    redeclare package SolidMat = Solid.Material.Crofer22APU,
    nX=nX,
    nY=nY,
    Tstartbar=TStart,
    lX=lX,
    lY=lY,
    lZ=lZsolid)
    annotation (Placement(transformation(extent={{-8,-124},{12,-104}})));
  Heat.Solid2D iCFuel(
    redeclare package SolidMat = Solid.Material.Crofer22APU,
    nX=nX,
    nY=nY,
    Tstartbar=TStart,
    lX=lX,
    lY=lY,
    lZ=lZsolid)
    annotation (Placement(transformation(extent={{-4,86},{16,106}})));

  // Rotation between PEN and AC
  Components.CrossFlowTopology crossFlowTopology(
    nX=nX,
    nY=nY,
    nSpecies=Medium.Air_Medium.nX,
    includeVarStream=true) if not simplification annotation (Placement(transformation(extent={{-14,-42},{16,-22}})));

  // Counter rotation between AC and IC
  Components.CrossFlowTopology crossFlowTopology1(
    nX=nY,
    nY=nX,
    nSpecies=Medium.Air_Medium.nX,
    includeVarStream=false) if not simplification annotation (Placement(transformation(extent={{-12,-100},{18,-80}})));

  // Manifolds
  Flow.Manifold                          manifoldFuel(redeclare package Medium =
        OpenTEMPEST.Medium.Fuel_CH4, nPorts_b=nY) if not simplification annotation (Placement(transformation(extent={{-76,48},{-68,68}})));
  Flow.Manifold_out                          manifoldFuel_out(redeclare package Medium =
        OpenTEMPEST.Medium.Fuel_CH4, nPorts_a=nY) if not simplification annotation (Placement(transformation(extent={{54,46},{74,66}})));
  Flow.Manifold                          manifoldAir(redeclare package Medium =
        OpenTEMPEST.Medium.Air_Medium, nPorts_b=nX) if not simplification annotation (Placement(transformation(extent={{-72,-68},{-64,-48}})));
  Flow.Manifold_out                          manifoldAir_out(redeclare package Medium =
        OpenTEMPEST.Medium.Air_Medium, nPorts_a=nX) if not simplification annotation (Placement(transformation(extent={{52,-70},{72,-50}})));

  // Electrical pins
  Modelica.Electrical.Analog.Interfaces.PositivePin pin_p annotation (Placement(transformation(extent={{-110,-16},{-90,4}})));
  Modelica.Electrical.Analog.Interfaces.NegativePin pin_n annotation (Placement(transformation(extent={{-110,16},{-90,36}})));

  // Flanges
  ThermoPower.Gas.FlangeA fuelIn(redeclare package Medium =
        OpenTEMPEST.Medium.Fuel_CH4) annotation (Placement(transformation(extent={{-112,50},{-92,70}})));
  ThermoPower.Gas.FlangeA airIn(redeclare package Medium =
        OpenTEMPEST.Medium.Air_Medium) annotation (Placement(transformation(extent={{-112,-70},{-92,-50}})));
  ThermoPower.Gas.FlangeB fuelOut(redeclare package Medium =
        OpenTEMPEST.Medium.Fuel_CH4) annotation (Placement(transformation(extent={{90,50},{110,70}})));
  ThermoPower.Gas.FlangeB airOut(redeclare package Medium =
        OpenTEMPEST.Medium.Air_Medium) annotation (Placement(transformation(extent={{90,-70},{110,-50}})));

  // Heat ports
  ThermoPower.Thermal.DHTVolumes dHT_x0(N=nY) annotation (Placement(transformation(extent={{-106,-34},{-94,-22}})));
  ThermoPower.Thermal.DHTVolumes dHT_xN(N=nY) annotation (Placement(transformation(extent={{94,0},{106,12}})));
  ThermoPower.Thermal.DHTVolumes dHT_y0(N=nX) annotation (Placement(transformation(extent={{-106,-92},{-94,-80}})));
  ThermoPower.Thermal.DHTVolumes dHT_yN(N=nX) annotation (Placement(transformation(extent={{-106,80},{-94,92}})));
  Heat.DHTVolumes2D dHT2_z1(i=nX, j=nY) annotation (Placement(transformation(extent={{36,94},{46,104}})));
  Heat.DHTVolumes2D dHT2_z0(i=nX, j=nY) annotation (Placement(transformation(extent={{36,-104},{46,-94}})));

  // Radiation objects
  Heat.RadHT2DFV radHT2DFV(
    nX=nX,
    nY=nY,
    specialCase=OpenTEMPEST.Enumerations.RadiationSpecialCases.case1,
    epsilon1=epsilonPEN,
    epsilon2=epsilonIC,
    A1=lX*Bac,
    A2=lX*Bac)
    annotation (Placement(transformation(extent={{46,-48},{66,-28}})));
  Heat.RadHT2DFV radHT2DFV1(
    nX=nX,
    nY=nY,
    specialCase=OpenTEMPEST.Enumerations.RadiationSpecialCases.case1,
    epsilon1=epsilonPEN,
    epsilon2=epsilonIC,
    A1=lX*Bfc,
    A2=lX*Bfc) annotation (Placement(transformation(extent={{46,26},{66,46}})));

equation

  if not simplification then
    // PEN stream variables
    connect(fuelChannel2D.PEN_in, pen.PEN_in) annotation (Line(points={{-16.24,50.36},
          {-14,50.36},{-14,34.4},{0,34.4}}, color={0,0,0}));
    // Manifolds - Channels
    connect(manifoldFuel.ports_b, fuelChannel2D.infl) annotation (Line(points={{-68,58},{-31.78,58.46}}, color={0,127,255}));
    connect(fuelChannel2D.outfl, manifoldFuel_out.ports_a) annotation (Line(points={{37.78,59},{50,59},{50,56},{61.4,56}}, color={159,159,223}));
    connect(manifoldAir.ports_b, airChannel2D.infl) annotation (Line(points={{-64,-58},{-35.54,-57.32}}, color={0,127,255}));
    connect(airChannel2D.outfl, manifoldAir_out.ports_a) annotation (Line(points={{41.54,-58},{41.54,-60},{59.4,-60}}, color={159,159,223}));
    // Manifolds - Flanges
    connect(manifoldFuel.port_a, fuelIn) annotation (Line(points={{-76,58},{-88,58},{-88,60},{-102,60}}, color={0,127,255}));
    connect(manifoldFuel_out.port_b, fuelOut) annotation (Line(points={{69,56},
            {86,56},{86,60},{100,60}},                                                                     color={0,127,255}));
    connect(manifoldAir_out.port_b,airOut)  annotation (Line(points={{67,-60},{100,-60}}, color={0,127,255}));
    connect(manifoldAir.port_a, airIn) annotation (Line(points={{-72,-58},{-88,-58},{-88,-60},{-102,-60}}, color={0,127,255}));
    // Crossflow topology connections
    connect(crossFlowTopology.dHTT_side2, airChannel2D.Q_PEN) annotation (Line(
          points={{1,-37},{1,-41.5},{2.59,-41.5},{2.59,-46.44}}, color={0,0,0}));
    connect(crossFlowTopology1.dHTT_side1, airChannel2D.Q_IC) annotation (Line(points={{3,-85},{3,-78.5},{3,-78.5},{3,-71.6}}, color={0,0,0}));
    // Crossflow topology connections
    connect(crossFlowTopology.dHTT_side1, pen.dhT2_z0) annotation (Line(points={{1,-27},
            {-11.5,-27},{-11.5,-10.4},{-23.8,-10.4}},      color={0,0,0}));
    connect(crossFlowTopology1.dHTT_side2, iCAir.dhT2_z1) annotation (Line(points=
           {{3,-95},{3,-101.5},{-5,-101.5},{-5,-107}}, color={0,0,0}));
    connect(pen.PEN_ina, crossFlowTopology.PEN_side1) annotation (Line(points={{0,-10.4},
            {-6,-10.4},{-6,-27},{-12.5,-27}},        color={0,0,0}));
    // z-direction
    connect(airChannel2D.PEN_in, crossFlowTopology.PEN_side2) annotation (Line(
          points={{-18.32,-47.12},{-18.32,-42.56},{-12.5,-42.56},{-12.5,-37}},
          color={0,0,0}));
    connect(iCFuel.dhT2_z0, fuelChannel2D.Q_IC) annotation (Line(points={{-1,89},{
          -1,79.5},{3,79.5},{3,69.8}}, color={0,0,0}));
    connect(fuelChannel2D.Q_PEN, pen.dhT2_z1) annotation (Line(points={{2.63,49.82},
          {4,49.82},{4,34.4},{-23.8,34.4}}, color={0,0,0}));
  else
    // PEN stream variables
    connect(fuelChannel1D.PEN_in, pen.PEN_in);
    connect(pen.PEN_ina, airChannel1D.PEN_in);
    // Flanges - Channels
    connect(fuelIn, fuelChannel1D.infl);
    connect(fuelChannel1D.outfl, fuelOut);
    connect(airIn, airChannel1D.infl);
    connect(airChannel1D.outfl,airOut);
    // z-direction
    connect(pen.dhT2_z0, airChannel1D.Q_PEN2D);
    connect(airChannel1D.Q_IC2D, iCAir.dhT2_z1);
    connect(iCFuel.dhT2_z0, fuelChannel1D.Q_IC2D);
    connect(fuelChannel1D.Q_PEN2D, pen.dhT2_z1);
  end if;

  // Connect to pins
  connect(pen.pin_n, pin_n) annotation (Line(points={{-34.34,23.2},{-84,23.2},{-84,26},{-100,26}}, color={0,0,255}));
  connect(pen.pin_p, pin_p) annotation (Line(points={{-34.34,-0.48},{-84,-0.48},{-84,-6},{-100,-6}}, color={0,0,255}));

  // x-direction
  connect(dHT_x0, pen.dhT_x0) annotation (Line(points={{-100,-28},{-54,-28},{-54,
          12},{-37.4,12}}, color={255,127,0}));
  connect(dHT_xN, pen.dhT_xN) annotation (Line(points={{100,6},{48,6},{48,12},{37.4,
          12}}, color={255,127,0}));

  // y-direction
  connect(dHT_y0, pen.dhT_y0);
  connect(dHT_yN, pen.dhT_yN);

  // z-direction
  connect(dHT2_z0, iCAir.dhT2_z0);
  connect(dHT2_z1, iCFuel.dhT2_z1);

  connect(iCAir.dhT2_int, radHT2DFV.side2) annotation (Line(points={{2,-114},{0,
          -114},{0,-104},{32,-104},{32,-80},{56,-80},{56,-38.4}}, color={0,0,0}));
  connect(pen.Qrad_AI, radHT2DFV.side1) annotation (Line(points={{17,6.24},{17,-24},
          {56,-24},{56,-33.6}}, color={0,0,0}));
  connect(pen.Qrad_FI, radHT2DFV1.side2)
    annotation (Line(points={{17,19.04},{56,19.04},{56,35.6}}, color={0,0,0}));
  connect(iCFuel.dhT2_int, radHT2DFV1.side1) annotation (Line(points={{6,96},{8,
          96},{8,74},{56,74},{56,40.4}}, color={0,0,0}));

  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -120},{100,100}}),                                  graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={28,108,200},
          fillColor={0,127,127},
          fillPattern=FillPattern.CrossDiag,
          lineThickness=0.5)}),                                  Diagram(coordinateSystem(preserveAspectRatio=false, extent={
            {-100,-120},{100,100}})));
end CellGeneral;

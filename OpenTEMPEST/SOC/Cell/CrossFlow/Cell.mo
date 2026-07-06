within OpenTEMPEST.SOC.Cell.CrossFlow;
model Cell

  import SI = Modelica.SIunits;

  parameter Integer nX=5 "Number of control volumes in the x-direction";
  parameter Integer nY=5 "Number of control volumes in the y-direction";

  // Dimensions from http://dx.doi.org/10.2139/ssrn.3987808
  parameter SI.Length lX = 0.1 "Length of solid" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lY = 0.1 "Width of solid" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lXpen = lX "Length of pen" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lYpen = lY "Width of pen" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lZpen = 3.425e-4 "Thickness of pen" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lXac = lX "Length of air channel" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lYac = lY "Width of air channel" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lZac = 1e-3 "Height of air channel" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lXfc = lX "Length of fuel channel" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lYfc = lY "Width of fuel channel" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lZfc = 1e-3 "Height of fuel channel" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lZsolid = 0.2e-3 "Height of interconnector"  annotation(Dialog(tab="Dimensions"));
  parameter SI.Length Bac = porAC*lYac "Width of air channel without ribs" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length Bfc = porFC*lYfc "Width of fuel channel without ribs" annotation(Dialog(tab="Dimensions"));

  // Initialization
  parameter SI.Temperature TStart=1023.15 "Uniform start temperature" annotation(Dialog(tab="Initialization"));
  parameter SI.AbsolutePressure pStart=101325 "Starting pressure" annotation(Dialog(tab="Initialization"));
  parameter SI.CurrentDensity Jstart = 0 "Starting current density in PEN" annotation(Dialog(tab="Initialization"));
  parameter SI.MassFraction xStartFC[FCMedium.nXi] = FCMedium.reference_X "Starting mass fraction in fuel channel" annotation(Dialog(tab="Initialization"));
  parameter SI.MassFraction xStartAC[ACMedium.nX] = ACMedium.reference_X "Starting mass fraction in air channel" annotation(Dialog(tab="Initialization"));

  // PEN parameters
  replaceable model Electrochem =
      OpenTEMPEST.SOC.Electrochem.Components.Crossflow_Electrochem
      constrainedby OpenTEMPEST.SOC.Electrochem.Components.ElectrochemBase         annotation (
  Dialog(tab="PEN"),
  Placement(transformation(extent={{78,50},{98,70}})),
  choicesAllMatching=true);
  parameter SI.ThermalConductivity kCustom_trans = 2.16 "Thermal Conductivity across layers of PEN in W/mK" annotation(Dialog(tab="PEN"));
  parameter SI.ThermalConductivity kCustom_long = 2.16 "Thermal Conductivity in plane of layers of PEN in W/mK (=k_trans for homogeneous materials)" annotation(Dialog(tab="PEN"));
  parameter SI.Density rhoPEN = 5900 "Density of PEN in kg/m3" annotation(Dialog(tab="PEN"));
  parameter SI.SpecificHeatCapacity cpPEN = 500 "Specific heat capacity of PEN in J/kgK" annotation(Dialog(tab="PEN"));
  parameter SI.SpectralEmissivity epsilonPEN = 0.8 "Emissivity of Anode-Electrolyte-Cathode unit" annotation(Dialog(tab="PEN"));

  // FC parameters
  replaceable package FCMedium = Medium.Fuel_CH4         annotation(Dialog(tab="Fuel Channel"));
  parameter Real porFC = 0.4 "Porosity in fuel channel" annotation(Dialog(tab="Fuel Channel"));
  parameter Real NuFCPEN = 12 "Nusselt number fuel channel on PEN side" annotation(Dialog(tab="Fuel Channel"));
  parameter Real NuFCIC = 10 "Nusselt number fuel channel on IC side" annotation(Dialog(tab="Fuel Channel"));

  // AC parameters
  replaceable package ACMedium = Medium.Air_Medium         annotation(Dialog(tab="Air Channel"));
  parameter Real porAC = 0.4 "Porosity in air channel" annotation(Dialog(tab="Air Channel"));
  parameter Real pDrop(max=0.99) = 0.04 "Pressure drop as a factor of inlet pressure (between 0 and 0.99)" annotation(Dialog(tab="Air Channel"));
  parameter Real NuACPEN = 8 "Nusselt number air channel on PEN side" annotation(Dialog(tab="Air Channel"));
  parameter Real NuACIC = 7.5 "Nusselt number air channel on IC side" annotation(Dialog(tab="Air Channel"));

  // IC parameters
  parameter SI.ThermalConductivity kIC = 40 "Thermal conductivity of CFY (35-45 W/mK for 20-900 °C)" annotation(Dialog(tab="Interconnects"));
  parameter SI.SpecificHeatCapacity cpIC = 451.8 "Interconnect heat capacity" annotation(Dialog(tab="Interconnects"));
  parameter SI.Density rhoIC = 7233   "Interconnect density" annotation(Dialog(tab="Interconnects"));
  parameter SI.SpectralEmissivity epsilonIC = 0.1 "Emissivity of interconnects" annotation(Dialog(tab="Interconnects"));

  PEN2D pen(
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
    redeclare model Electrochem = Electrochem)
    annotation (Placement(transformation(extent={{-34,-20},{34,44}})));
  OpenTEMPEST.SOC.Cell.CrossFlow.AirChannel2D airChannel(
    nX=nY,
    nY=nX,
    TStart=TStart,
    pStart=pStart,
    xStart=xStartAC,
    lX=lYac,
    lY=lXac,
    lZ=lZac,
    por=porAC,
    kRibs=kIC,
    cpRibs=cpIC,
    rhoRibs=rhoIC,
    pDrop=pDrop,
    Nu_PEN=NuACPEN,
    Nu_IC=NuACIC)
    annotation (Placement(transformation(extent={{-38,-24},{44,-92}})));
  CrossFlowTopology crossFlowTopology(
    nX=nX,
    nY=nY,
    nSpecies=Medium.Air_Medium.nX,
    includeVarStream=true)
    annotation (Placement(transformation(extent={{-14,-42},{16,-22}})));
  OpenTEMPEST.SOC.Cell.CrossFlow.FuelChannel2D fuelChannel(
    nX=nX,
    nY=nY,
    TStart=TStart,
    pStart=pStart,
    xStart=xStartFC,
    lX=lXfc,
    lY=lYfc,
    lZ=lZfc,
    por=porFC,
    Nu_PEN=NuFCPEN,
    Nu_IC=NuFCIC,
    kRibs=kIC,
    cpRibs=cpIC,
    rhoRibs=rhoIC)
    annotation (Placement(transformation(extent={{-34,32},{40,86}})));

  Flow.Manifold manifold(redeclare package Medium = Medium.Fuel_CH4, nPorts_b=
        nY) annotation (Placement(transformation(extent={{-76,48},{-68,68}})));
  Flow.Manifold_out manifold_out(redeclare package Medium = Medium.Fuel_CH4,
      nPorts_a=nY)
    annotation (Placement(transformation(extent={{54,46},{74,66}})));
  Flow.Manifold manifold1(redeclare package Medium = Medium.Air_Medium,
      nPorts_b=nX)
    annotation (Placement(transformation(extent={{-72,-68},{-64,-48}})));
  Flow.Manifold_out manifold_out1(redeclare package Medium = Medium.Air_Medium,
      nPorts_a=nX)
    annotation (Placement(transformation(extent={{52,-70},{72,-50}})));
  Modelica.Electrical.Analog.Interfaces.PositivePin pin_p annotation (Placement(transformation(extent={{-106,0},
            {-86,20}}), iconTransformation(extent={{-106,0},{-86,20}})));
  Modelica.Electrical.Analog.Interfaces.NegativePin pin_n annotation (Placement(transformation(extent={{-106,40},
            {-86,60}}), iconTransformation(extent={{-106,40},{-86,60}})));
  ThermoPower.Gas.FlangeA fuelIn(redeclare package Medium = Medium.Fuel_CH4)
    annotation (Placement(transformation(extent={{-112,68},{-92,88}}),
        iconTransformation(extent={{-112,68},{-92,88}})));
  ThermoPower.Gas.FlangeA airIn(redeclare package Medium = Medium.Air_Medium)
    annotation (Placement(transformation(extent={{-110,-110},{-90,-90}}),
        iconTransformation(extent={{-110,-110},{-90,-90}})));
  ThermoPower.Gas.FlangeB fuelOut(redeclare package Medium = Medium.Fuel_CH4)
    annotation (Placement(transformation(extent={{90,50},{110,70}})));
  ThermoPower.Gas.FlangeB airout(redeclare package Medium = Medium.Air_Medium)
    annotation (Placement(transformation(extent={{90,-70},{110,-50}})));
  OpenTEMPEST.Heat.Solid2D iCAir(
    redeclare package SolidMat = OpenTEMPEST.Solid.Material.Crofer22APU,
    nX=nX,
    nY=nY,
    Tstartbar=TStart,
    lX=lX,
    lY=lY,
    lZ=lZsolid)
    annotation (Placement(transformation(extent={{-8,-124},{12,-104}})));
  OpenTEMPEST.Heat.Solid2D iCFuel(
    redeclare package SolidMat = OpenTEMPEST.Solid.Material.Crofer22APU,
    nX=nX,
    nY=nY,
    Tstartbar=TStart,
    lX=lX,
    lY=lY,
    lZ=lZsolid)
    annotation (Placement(transformation(extent={{-4,86},{16,106}})));
  ThermoPower.Thermal.DHTVolumes dHT_x0(N=nY)
    annotation (Placement(transformation(extent={{-100,-60},{-88,-20}}),
        iconTransformation(extent={{-100,-60},{-88,-20}})));
  ThermoPower.Thermal.DHTVolumes dHT_y0(N=nX)
    annotation (Placement(transformation(extent={{-52,-120},{58,-108}}),
        iconTransformation(extent={{-52,-120},{58,-108}})));
  ThermoPower.Thermal.DHTVolumes dHT_yN(N=nX)
    annotation (Placement(transformation(extent={{-54,88},{56,100}}),
        iconTransformation(extent={{-54,88},{56,100}})));
  Heat.DHTVolumes2D dHT2_z1(i=nX, j=nY) annotation (Placement(transformation(
          extent={{-24,114},{36,126}}), iconTransformation(extent={{-54,34},{56,
            58}})));
  Heat.DHTVolumes2D dHT2_z0(i=nX, j=nY) annotation (Placement(transformation(
          extent={{-16,-146},{20,-138}}), iconTransformation(extent={{-52,-68},
            {58,-44}})));
  ThermoPower.Thermal.DHTVolumes dHT_xN(N=nY)
    annotation (Placement(transformation(extent={{86,-20},{100,20}}),
        iconTransformation(extent={{86,-20},{100,20}})));
  CrossFlowTopology crossFlowTopology1(
    nX=nY,
    nY=nX,
    nSpecies=Medium.Air_Medium.nX,
    includeVarStream=false)
    annotation (Placement(transformation(extent={{-12,-100},{18,-80}})));

  OpenTEMPEST.Heat.RadHT2DFV radHT2DFV(
    nX=nX,
    nY=nY,
    specialCase=OpenTEMPEST.Enumerations.RadiationSpecialCases.case1,
    epsilon1=epsilonPEN,
    epsilon2=epsilonIC,
    A1=lX*Bac,
    A2=lX*Bac)
    annotation (Placement(transformation(extent={{46,-48},{66,-28}})));
  OpenTEMPEST.Heat.RadHT2DFV radHT2DFV1(
    nX=nX,
    nY=nY,
    specialCase=OpenTEMPEST.Enumerations.RadiationSpecialCases.case1,
    epsilon1=epsilonPEN,
    epsilon2=epsilonIC,
    A1=lX*Bfc,
    A2=lX*Bfc) annotation (Placement(transformation(extent={{46,26},{66,46}})));
equation

  // Pins connections
  connect(pen.pin_n, pin_n) annotation (Line(points={{-34.34,23.2},{-84,23.2},{-84,
          50},{-96,50}},                                                                           color={0,0,255}));
  connect(pen.pin_p, pin_p) annotation (Line(points={{-34.34,-0.48},{-84,-0.48},
          {-84,10},{-96,10}},                                                                        color={0,0,255}));

  // Flanges connections
  connect(manifold.port_a, fuelIn) annotation (Line(points={{-76,58},{-88,58},{-88,
          78},{-102,78}},                                                                          color={0,127,255}));
  connect(manifold_out.port_b, fuelOut) annotation (Line(points={{69,56},{86,56},{86,60},{100,60}}, color={0,127,255}));
  connect(manifold_out1.port_b, airout) annotation (Line(points={{67,-60},{100,-60}}, color={0,127,255}));
  connect(manifold1.port_a, airIn) annotation (Line(points={{-72,-58},{-88,-58},
          {-88,-100},{-100,-100}},                                                                     color={0,127,255}));

  // x-direction connections
  connect(dHT_x0, pen.dhT_x0);
  connect(dHT_xN, pen.dhT_xN);

  // y-direction connections
  connect(dHT_y0, pen.dhT_y0);
  connect(dHT_yN, pen.dhT_yN);

  // z-direction connections
  connect(dHT2_z0, iCAir.dhT2_z0);
  connect(dHT2_z1, iCFuel.dhT2_z1);
  connect(iCFuel.dhT2_z0, fuelChannel.Q_IC) annotation (Line(points={{-1,89},{-1,
          79.5},{3,79.5},{3,69.8}}, color={0,0,0}));
  connect(fuelChannel.Q_PEN, pen.dhT2_z1) annotation (Line(points={{2.63,49.82},
          {4,49.82},{4,34.4},{-23.8,34.4}}, color={0,0,0}));

  connect(iCAir.dhT2_int, radHT2DFV.side2) annotation (Line(points={{2,-114},{0,
          -114},{0,-104},{32,-104},{32,-80},{56,-80},{56,-38.4}}, color={0,0,0}));
  connect(pen.Qrad_AI, radHT2DFV.side1) annotation (Line(points={{17,6.24},{17,-24},
          {56,-24},{56,-33.6}}, color={0,0,0}));
  connect(pen.Qrad_FI, radHT2DFV1.side2)
    annotation (Line(points={{17,19.04},{56,19.04},{56,35.6}}, color={0,0,0}));
  connect(iCFuel.dhT2_int, radHT2DFV1.side1) annotation (Line(points={{6,96},{8,
          96},{8,74},{56,74},{56,40.4}}, color={0,0,0}));

  connect(fuelChannel.PEN_in, pen.PEN_in) annotation (Line(points={{-16.24,
          50.36},{-14,50.36},{-14,34.4},{0,34.4}}, color={0,0,0}));

  // Connect to Manifolds
  connect(manifold.ports_b, fuelChannel.infl)
    annotation (Line(points={{-68,58},{-31.78,58.46}}, color={0,127,255}));
  connect(fuelChannel.outfl, manifold_out.ports_a) annotation (Line(points={{
          37.78,59},{50,59},{50,56},{61.4,56}}, color={159,159,223}));
  connect(manifold1.ports_b, airChannel.infl)
    annotation (Line(points={{-64,-58},{-35.54,-57.32}}, color={0,127,255}));
  connect(airChannel.outfl, manifold_out1.ports_a) annotation (Line(points={{
          41.54,-58},{41.54,-60},{59.4,-60}}, color={159,159,223}));

  // Crossflow topology connections
  connect(crossFlowTopology.dHTT_side1, pen.dhT2_z0) annotation (Line(points={{1,-27},
          {-11.5,-27},{-11.5,-10.4},{-23.8,-10.4}},      color={0,0,0}));
  connect(crossFlowTopology.dHTT_side2, airChannel.Q_PEN) annotation (Line(
        points={{1,-37},{1,-41.5},{2.59,-41.5},{2.59,-46.44}}, color={0,0,0}));
  connect(crossFlowTopology1.dHTT_side1, airChannel.Q_IC) annotation (Line(
        points={{3,-85},{3,-78.5},{3,-78.5},{3,-71.6}}, color={0,0,0}));
  connect(crossFlowTopology1.dHTT_side2, iCAir.dhT2_z1) annotation (Line(points=
         {{3,-95},{3,-101.5},{-5,-101.5},{-5,-107}}, color={0,0,0}));
  connect(airChannel.PEN_in, crossFlowTopology.PEN_side2) annotation (Line(
        points={{-18.32,-47.12},{-18.32,-42.56},{-12.5,-42.56},{-12.5,-37}},
        color={0,0,0}));
  connect(pen.PEN_ina, crossFlowTopology.PEN_side1) annotation (Line(points={{0,-10.4},
          {-6,-10.4},{-6,-27},{-12.5,-27}},        color={0,0,0}));


  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -120},{100,100}}),                                  graphics={Rectangle(
          extent={{-100,100},{100,-120}},
          lineColor={28,108,200},
          fillColor={0,127,127},
          fillPattern=FillPattern.CrossDiag,
          lineThickness=0.5)}),                                  Diagram(coordinateSystem(preserveAspectRatio=false, extent={
            {-100,-120},{100,100}})),
    Documentation(info="<html>
<h2>2D Crossflow Cell</h2>
<p>
The model represents a two-dimensional crossflow cell discretized 
into control volumes (CVs) in the x- and y-directions.
</p>

<p>
The cell consists of a discretized PEN, perpendicular air and fuel channels, 
and metallic interconnects. The electrochemistry is provided by a replaceable 
submodel, and gas manifolds distribute the reactants to the channels. 
Geometrical dimensions, material properties, discretization, and transport 
parameters are fully parameterized.
</p>
</html>"));
end Cell;

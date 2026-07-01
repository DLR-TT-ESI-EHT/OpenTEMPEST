within OpenTEMPEST.SOC.Cell.CrossFlow;
model CellCFSimp
  "Simplified crossflow cell with 1D fuel and air channels, 2D heat ports and variable stream connectors directly in channel models."

  import SI = Modelica.SIunits;

  // General parameters
  parameter Integer nX=3 "Number of control volumes in first direction";
  parameter Integer nY=3 "Number of control volumes in second direction";
  parameter Boolean LUDS=false  "Set true if Linear upwind difference wanted (more accuracy), false for Upwind difference scheme (more speed and stability)";
  parameter Boolean heatTransferCorrelationFormDuct=true "true for Nusselt correlation duct geometry with characteristic length=2*lZ (default), false for plate geometry with characteristic length=lX";

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
  parameter SI.Temperature TStart=1023.15 "Uniform starting temperature" annotation(Dialog(tab="Initialization"));
  parameter SI.AbsolutePressure pStart=101325 "Starting pressure" annotation(Dialog(tab="Initialization"));
  parameter SI.CurrentDensity Jstart = 0 "Starting current density in PEN" annotation(Dialog(tab="Initialization"));
  parameter SI.MassFraction xStartFC[FCMedium.nXi] = FCMedium.reference_X "Starting mass fraction in fuel channel" annotation(Dialog(tab="Initialization"));
  parameter SI.MassFraction xStartAC[ACMedium.nX] = ACMedium.reference_X "Starting mass fraction in air channel" annotation(Dialog(tab="Initialization"));

  // PEN parameters
  replaceable model Electrochem =
      OpenTEMPEST.SOC.Electrochem.Components.Crossflow_Electrochem
      constrainedby OpenTEMPEST.SOC.Electrochem.Components.ElectrochemBase
                                                                   "Electrochemical model used" annotation (
  Dialog(tab="PEN"),
  Placement(transformation(extent={{78,50},{98,70}})),
  choicesAllMatching=true);
  parameter SI.ThermalConductivity kCustom_trans = 2.16 "Thermal Conductivity across layers of PEN in W/mK" annotation(Dialog(tab="PEN"));
  parameter SI.ThermalConductivity kCustom_long = 2.16 "Thermal Conductivity in plane of layers of PEN in W/mK (=k_trans for homogeneous materials)" annotation(Dialog(tab="PEN"));
  parameter SI.Density rhoPEN=5900   "Density of PEN in kg/m3" annotation(Dialog(tab="PEN"));
  parameter SI.SpecificHeatCapacity cpPEN = 500 "Specific heat capacity of PEN in J/kgK" annotation(Dialog(tab="PEN"));
  parameter SI.SpectralEmissivity epsilonPEN = 0.8 "Emissivity of Anode-Electrolyte-Cathode unit" annotation(Dialog(tab="PEN"));

  // FC parameters
  replaceable package FCMedium = Medium.Fuel_CH4         annotation(Dialog(tab="Fuel Channel"));
  parameter Real porFC = 0.4 "Porosity in fuel channel" annotation(Dialog(tab="Fuel Channel"));
  parameter Real NuFCPEN = 12 "Nusselt number fuel channel on PEN side" annotation(Dialog(tab="Fuel Channel"));
  parameter Real NuFCIC = 9.86 "Nusselt number fuel channel on IC side" annotation(Dialog(tab="Fuel Channel"));
  parameter Real alfafc=0.015  "Weight for 2D temperature in z-direction convection for fuel channel" annotation(Dialog(tab="Fuel Channel"));

  // AC parameters
  replaceable package ACMedium = Medium.Air_Medium         annotation(Dialog(tab="Air Channel"));
  parameter Real porAC = 0.4 "Porosity in air channel" annotation(Dialog(tab="Air Channel"));
  parameter Real pDrop(max=0.99) = 0.04 "Pressure drop as a factor of inlet pressure (between 0 and 0.99)" annotation(Dialog(tab="Air Channel"));
  parameter Real NuACPEN = 8 "Nusselt number air channel on PEN side" annotation(Dialog(tab="Air Channel"));
  parameter Real NuACIC = 7.5 "Nusselt number air channel on IC side" annotation(Dialog(tab="Air Channel"));
  parameter Real alfaac=0.12   "Weight for 2D temperature in z-direction convection for air channel" annotation(Dialog(tab="Air Channel"));

  // IC parameters
  parameter SI.ThermalConductivity kIC = 40 "Thermal conductivity of CFY (35-45 W/mK for 20-900 °C)" annotation(Dialog(tab="Interconnects"));
  parameter SI.SpecificHeatCapacity cpIC = 451.8 "Interconnect heat capacity" annotation(Dialog(tab="Interconnects"));
  parameter SI.Density rhoIC = 7233   "Interconnect density" annotation(Dialog(tab="Interconnects"));
  parameter SI.SpectralEmissivity epsilonIC = 0.1 "Emissivity of interconnects" annotation(Dialog(tab="Interconnects"));

  OpenTEMPEST.SOC.Cell.CrossFlow.FuelChannel2DSimp fuelChannel(
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
    alfa=alfafc,
    heatTransferCorrelationFormDuct=heatTransferCorrelationFormDuct)
    annotation (Placement(transformation(extent={{-34,68},{34,116}})));

  OpenTEMPEST.SOC.Cell.CrossFlow.AirChannel2DSimp airChannel(
    nX=nY,
    nY=nX,
    TStart=TStart,
    pStart=pStart,
    xStart=xStartAC,
    lX=lYac,
    lY=lXac,
    lZ=lZac,
    por=porAC,
    Nu_PEN=NuACPEN,
    Nu_IC=NuACIC,
    alfa=alfaac,
    pDrop=0.01,
    heatTransferCorrelationFormDuct=heatTransferCorrelationFormDuct)
    annotation (Placement(transformation(extent={{-36,-48},{34,-96}})));

  OpenTEMPEST.Heat.Solid2D iCFuel(
    redeclare package SolidMat = OpenTEMPEST.Solid.Material.Crofer22APU,
    nX=nX,
    nY=nY,
    Tstartbar=TStart,
    lX=lX,
    lY=lY,
    lZ=lZsolid)
    annotation (Placement(transformation(extent={{-16,124},{16,156}})));
  OpenTEMPEST.Heat.Solid2D iCAir(
    redeclare package SolidMat = OpenTEMPEST.Solid.Material.Crofer22APU,
    nX=nX,
    nY=nY,
    Tstartbar=TStart,
    lX=lX,
    lY=lY,
    lZ=lZsolid)
    annotation (Placement(transformation(extent={{-16,-148},{16,-114}})));

  Modelica.Electrical.Analog.Interfaces.PositivePin pin_p annotation (Placement(transformation(extent={{-148,18},
            {-126,40}}),       iconTransformation(extent={{-148,18},{-126,40}})));
  Modelica.Electrical.Analog.Interfaces.NegativePin pin_n annotation (Placement(transformation(extent={{-148,58},
            {-126,80}}), iconTransformation(extent={{-148,58},{-126,80}})));
  ThermoPower.Gas.FlangeA fuelIn(redeclare package Medium = Medium.Fuel_CH4)
    annotation (Placement(transformation(extent={{-148,100},{-128,120}}),
        iconTransformation(extent={{-148,100},{-128,120}})));
  ThermoPower.Gas.FlangeA airIn(redeclare package Medium = Medium.Air_Medium)
    annotation (Placement(transformation(extent={{-150,-118},{-130,-98}}),
        iconTransformation(extent={{-150,-118},{-130,-98}})));
  ThermoPower.Gas.FlangeB fuelOut(redeclare package Medium = Medium.Fuel_CH4)
    annotation (Placement(transformation(extent={{128,74},{148,94}})));
  ThermoPower.Gas.FlangeB airout(redeclare package Medium = Medium.Air_Medium)
    annotation (Placement(transformation(extent={{130,-86},{150,-66}}),
        iconTransformation(extent={{130,-86},{150,-66}})));
  ThermoPower.Thermal.DHTVolumes dHT_x0(N=nY)
    annotation (Placement(transformation(extent={{-140,-60},{-122,-20}}),
        iconTransformation(extent={{-140,-60},{-122,-20}})));
  ThermoPower.Thermal.DHTVolumes dHT_y0(N=nX)
    annotation (Placement(transformation(extent={{-6,-6},{6,6}}, origin={140,-26}),
        iconTransformation(extent={{-70,-140},{72,-120}})));
  ThermoPower.Thermal.DHTVolumes dHT_yN(N=nX)
    annotation (Placement(transformation(extent={{-6,-6},{6,6}}, origin={140,26}),
        iconTransformation(extent={{-70,122},{72,140}})));
  ThermoPower.Thermal.DHTVolumes dHT_xN(N=nY)
    annotation (Placement(transformation(extent={{-5,-16},{5,16}},
                                                                 origin={135,4}),
        iconTransformation(extent={{122,-20},{140,20}})));
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
    annotation (Placement(transformation(extent={{-36,-6},{24,44}})));

  Heat.DHTVolumes2D dHT2_z0(i=nX, j=nY) annotation (Placement(transformation(
          extent={{-90,-118},{-40,-96}}), iconTransformation(extent={{-82,-60},
            {86,-38}})));

  Heat.DHTVolumes2D dHT2_z1(i=nX, j=nY) annotation (Placement(transformation(
          extent={{48,118},{174,138}}), iconTransformation(extent={{-82,52},{86,
            74}})));

  CrossFlowTopology crossFlowTopologyPENAC(
    nX=nX,
    nY=nY,
    nSpecies=Medium.Air_Medium.nX,
    includeVarStream=true)
    annotation (Placement(transformation(extent={{-14,-44},{22,-18}})));
  CrossFlowTopology crossFlowTopologyACIC(
    nX=nY,
    nY=nX,
    nSpecies=Medium.Air_Medium.nX,
    includeVarStream=false)
    annotation (Placement(transformation(extent={{-6,-110},{32,-84}})));
equation

  // Pins connections
  connect(pin_n, pen.pin_n) annotation (Line(points={{-137,69},{-48,69},{-48,27.75},
          {-36.3,27.75}},        color={0,0,255}));
  connect(pin_p, pen.pin_p) annotation (Line(points={{-137,29},{-42,29},{-42,9.25},
          {-36.3,9.25}},       color={0,0,255}));

  // Flanges connections
  connect(fuelChannel.infl, fuelIn) annotation (Line(points={{-34,93.44},{-78,
          93.44},{-78,110},{-138,110}}, color={159,159,223}));
  connect(fuelChannel.outfl, fuelOut) annotation (Line(points={{34,93.44},{62,
          93.44},{62,90},{138,90},{138,84}}, color={159,159,223}));
  connect(airIn, airChannel.infl) annotation (Line(points={{-140,-108},{-106,-108},
          {-106,-76},{-36,-76},{-36,-73.44}}, color={159,159,223}));
  connect(airChannel.outfl, airout) annotation (Line(points={{34,-73.44},{140,-73.44},
          {140,-76}}, color={159,159,223}));

  // x-direction connections
  connect(dHT_x0, pen.dhT_x0) annotation (Line(points={{-131,-40},{-48,-40},{-48,
          19},{-39,19}},
                     color={255,127,0}));
  connect(dHT_xN, pen.dhT_xN) annotation (Line(points={{135,4},{38,4},{38,19},{27,
          19}},    color={255,127,0}));

  // y-direction connections
  connect(pen.dhT_y0, dHT_y0)
    annotation (Line(points={{15,1.5},{15,-6},{68,-6},{68,-18},{132,-18},{132,-26},
          {140,-26}},                                      color={255,127,0}));
  connect(dHT_yN, pen.dhT_yN) annotation (Line(points={{140,26},{36,26},{36,
          36.5},{15,36.5}}, color={255,127,0}));

  // z-direction connections
  connect(iCFuel.dhT2_z1, dHT2_z1) annotation (Line(points={{-11.2,151.2},{-11.2,
          160},{42,160},{42,128},{111,128}},
                                           color={0,0,0}));
  connect(iCAir.dhT2_z0, dHT2_z0) annotation (Line(points={{-11.2,-142.9},{-65,-142.9},
          {-65,-107}},                      color={0,0,0}));
  connect(crossFlowTopologyACIC.dHTT_side2, iCAir.dhT2_z1) annotation (Line(
        points={{13,-103.5},{13,-110},{-11.2,-110},{-11.2,-119.1}}, color={0,0,0}));

  connect(pen.PEN_in, fuelChannel.PEN_in) annotation (Line(points={{-6,36.5},{-6,
          64},{-17,64},{-17,84.8}}, color={0,0,0}));

  connect(crossFlowTopologyPENAC.dHTT_side2, airChannel.Q_PEN) annotation (
      Line(points={{4,-37.5},{32,-37.5},{32,-63.84},{-1,-63.84}}, color={0,0,0}));
  connect(crossFlowTopologyPENAC.PEN_side2, airChannel.PEN_in) annotation (Line(
        points={{-12.2,-37.5},{-18.5,-37.5},{-18.5,-64.8}}, color={0,0,0}));
  connect(pen.PEN_ina, crossFlowTopologyPENAC.PEN_side1) annotation (Line(
        points={{-6,1.5},{-6,-14},{-20,-14},{-20,-24.5},{-12.2,-24.5}}, color={
          0,0,0}));
  connect(iCFuel.dhT2_z0, fuelChannel.Q_IC) annotation (Line(points={{-11.2,
          128.8},{-11.2,120},{0,120},{0,101.6}}, color={0,0,0}));
  connect(fuelChannel.Q_PEN, pen.dhT2_z1) annotation (Line(points={{0,83.84},
          {0,64},{-8,64},{-8,48},{-30,48},{-30,36.5},{-27,36.5}}, color={0,0,0}));
  connect(pen.dhT2_z0, crossFlowTopologyPENAC.dHTT_side1) annotation (Line(
        points={{-27,1.5},{-28,1.5},{-28,-10},{4,-10},{4,-24.5}}, color={0,0,0}));
  connect(crossFlowTopologyACIC.dHTT_side1, airChannel.Q_IC) annotation (Line(
        points={{13,-90.5},{56,-90.5},{56,-84},{-2,-84},{-2,-81.6},{-1,-81.6}},
        color={0,0,0}));

  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-140,-140},
            {140,140}}),                                        graphics={Rectangle(
          extent={{-140,140},{140,-140}},
          lineColor={28,108,200},
          fillColor={0,127,127},
          fillPattern=FillPattern.CrossDiag,
          lineThickness=0.5)}),                                  Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-140,
            -140},{140,140}})),
    Documentation(info="<html>
<p>
This model represents a simplified crossflow cell model based on the detailed cell model already present in the library
<a href=\"TEMPEST.ECReactorModels.Cell.CrossFlow.Cell\">Detailed CF Cell Model</a>.
<p>
This simplification approach consists in combining 1D discretised air and fuel channel in their respective flow direction with 2D discretised PEN and interconnect models.
The adaptations needed for connecting and interfacing the 2D models with the 1D models are managed internally, specifically inside the channel models.
</p>
</html>", revisions="<html>
</html>"));
end CellCFSimp;

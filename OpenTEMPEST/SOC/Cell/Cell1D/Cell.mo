within OpenTEMPEST.SOC.Cell.Cell1D;
model Cell
  import SI = Modelica.SIunits;
  replaceable package FCMedium = OpenTEMPEST.Medium.Fuel_CH4
                                                         annotation (Dialog(tab="Fuel Channel"));
  replaceable package ACMedium = OpenTEMPEST.Medium.Air_Medium
                                                           annotation (Dialog(tab="Air Channel"));

  // General parameters
  parameter Integer N(min=3) "number of CVs in each layer" annotation (Dialog(tab="General"));
  parameter SI.Length lX=0.1 annotation (Dialog(tab="General"));
  parameter SI.Length lY=0.15 annotation (Dialog(tab="General"));
  constant SI.AbsolutePressure p0=1e5 annotation (Dialog(tab="General"));
  replaceable function fluxInterp =
      OpenTEMPEST.Flow.FluxInterpolators.UDSinterp  constrainedby OpenTEMPEST.Flow.FluxInterpolators.DifferencingSchemeInterpBase
                                                                                         annotation(choicesAllMatching = true, Dialog(tab="General"));

  //initialization
  parameter SI.Temperature TStart=1073.15 annotation (Dialog(tab="Initialization"));
  parameter SI.Temperature TStart_pen=TStart annotation (Dialog(tab="Initialization"));
  parameter SI.Temperature TStart_FCin=TStart annotation (Dialog(tab="Initialization"));
  parameter SI.Temperature TStart_FCout=TStart annotation (Dialog(tab="Initialization"));
  parameter SI.Temperature TStart_ACin=TStart annotation (Dialog(tab="Initialization"));
  parameter SI.Temperature TStart_ACout=TStart annotation (Dialog(tab="Initialization"));
  parameter SI.Temperature TStart_IC=TStart annotation (Dialog(tab="Initialization"));
  parameter SI.AbsolutePressure pStartAC=101325 annotation (Dialog(tab="Initialization"));
  parameter SI.AbsolutePressure pStartInFC=101325 annotation (Dialog(tab="Initialization"));
  parameter SI.AbsolutePressure pStartOutFC=pStartInFC annotation (Dialog(tab="Initialization"));
  parameter SI.CurrentDensity Jstart annotation (Dialog(tab="Initialization"));
  parameter SI.MassFraction xStartAC[ACMedium.nXi]=ACMedium.reference_X annotation (Dialog(tab="Initialization"));
  parameter SI.MassFraction xStartFC[FCMedium.nX]=FCMedium.reference_X annotation (Dialog(tab="Initialization"));
  parameter SI.PressureDifference dpNomFC=0.01 "nominal pressure loss for initialization" annotation (Dialog(tab="Initialization"));

  // PEN parameters
  replaceable model Electrochem =
      OpenTEMPEST.SOC.Electrochem.Components.BV_Steam
  constrainedby OpenTEMPEST.SOC.Electrochem.Components.ElectrochemBase
     annotation (Placement(transformation(extent={{78,50},{98,70}})),choicesAllMatching=true);
  parameter SI.Length lXpen=lX annotation (Dialog(tab="PEN"));
  parameter SI.Length lYpen=lY annotation (Dialog(tab="PEN"));
  parameter SI.Length lZpen=1.51e-4 annotation (Dialog(tab="PEN"));
  parameter SI.ThermalConductivity kPen=2 annotation (Dialog(tab="PEN"));
  parameter SI.Density rhoPen=5900 annotation (Dialog(tab="PEN"));
  parameter SI.SpecificHeatCapacity cpPen=500.0 annotation (Dialog(tab="PEN"));
  parameter SI.SpectralEmissivity epsilon_pen=0.8 "emissivity of Anode-Electrolyte-Cathode unit" annotation (Dialog(tab="PEN"));
  // FC parameters
  parameter SI.Length lXfc=lX annotation (Dialog(tab="Fuel Channel"));
  parameter SI.Length lYfc=lY annotation (Dialog(tab="Fuel Channel"));
  parameter SI.Length lZfc=1e-3 annotation (Dialog(tab="Fuel Channel"));
  parameter Real porFC=0.9 annotation (Dialog(tab="Fuel Channel"));
  parameter Real Nu_PENfc=12 "Nusselt number fuel channel on PEN side" annotation (Dialog(tab="Fuel Channel"));
  parameter Real Nu_ICfc=10 "Nusselt Number fuel channel on IC side" annotation (Dialog(tab="Fuel Channel"));
  parameter SI.ThermalConductivity kFoam=4 annotation (Dialog(tab="Fuel Channel"));
  parameter SI.SpecificHeatCapacity cpFoam=400 annotation (Dialog(tab="Fuel Channel"));
  parameter SI.Density rhoFoam=9000 annotation (Dialog(tab="Fuel Channel"));
  parameter ThermoPower.Units.HydraulicResistance HyR_Fuel=50000000 "Hydraulic resistance of fuel channel of each single cell in stack"
    annotation (Dialog(tab="Fuel Channel"));

  // AC parameters
  parameter SI.Length lXac=lX annotation (Dialog(tab="Air Channel"));
  parameter SI.Length lYac=0.08 annotation (Dialog(tab="Air Channel"));
  parameter SI.Length lZac=1e-3 annotation (Dialog(tab="Air Channel"));
  parameter Real porAC=1 annotation (Dialog(tab="Air Channel"));
  parameter Real pDropAC=0.04 annotation (Dialog(tab="Air Channel"));
  parameter Real Nu_PENac=8 "Nusselt number air channel on PEN side" annotation (Dialog(tab="Air Channel"));
  parameter Real Nu_ICac=8 "Nusselt Number air channel on IC side" annotation (Dialog(tab="Air Channel"));
  parameter ThermoPower.Units.HydraulicResistance HyR_Air=6000000 "Hydraulic resistance of air channel of each single cell in stack" annotation (Dialog(tab="Air Channel"));

  //IC parameters
  parameter SI.SpectralEmissivity epsilon_ic=0.1 "emissivity of surface of IC" annotation (Dialog(tab="Interconnect"));
  parameter SI.Length lZIC=3e-4 "half of IC height" annotation (Dialog(tab="Interconnect"));
  parameter SI.Density rhoIC=8000 "Density of IC unit" annotation (Dialog(tab="Interconnect"));
  parameter SI.SpecificHeatCapacity cpIC=500 "specific heat capcity of IC" annotation (Dialog(tab="Interconnect"));
  parameter SI.ThermalConductivity LambdaIC=25 "thermal conductivity of IC" annotation (Dialog(tab="Interconnect"));

  ThermoPower.Gas.FlangeA fuelFlangeIn(redeclare package Medium = FCMedium)
    annotation (Placement(transformation(extent={{-110,44},{-90,64}}), iconTransformation(extent={{-110,44},{-90,64}})));
  ThermoPower.Gas.FlangeB fuelFlangeOut(redeclare package Medium = FCMedium)
    annotation (Placement(transformation(extent={{90,44},{110,64}}), iconTransformation(extent={{90,44},{110,64}})));
  ThermoPower.Gas.FlangeA airFlangeIn(redeclare package Medium = ACMedium)
    annotation (Placement(transformation(extent={{-110,-76},{-90,-56}}), iconTransformation(extent={{-110,-76},{-90,-56}})));
  ThermoPower.Gas.FlangeB airFlangeOut(redeclare package Medium = ACMedium)
    annotation (Placement(transformation(extent={{90,-58},{110,-38}}), iconTransformation(extent={{90,-58},{110,-38}})));
  Modelica.Electrical.Analog.Interfaces.PositivePin pin_p
    annotation (Placement(transformation(extent={{-106,-54},{-86,-34}}), iconTransformation(extent={{-106,-54},{-86,-34}})));
  Modelica.Electrical.Analog.Interfaces.NegativePin pin_n
    annotation (Placement(transformation(extent={{-106,-32},{-86,-12}}), iconTransformation(extent={{-106,-32},{-86,-12}})));
  ThermoPower.Thermal.DHTVolumes dHT_z1(N=N) annotation (Placement(transformation(extent={{-60,-80},{-40,-60}}), iconTransformation(extent={{-62,64},{-42,84}})));
  ThermoPower.Thermal.DHTVolumes dHT_z0(N=N) annotation (Placement(transformation(extent={{-60,72},{-40,92}}), iconTransformation(extent={{-60,-80},{-40,-60}})));
  ThermoPower.Thermal.DHTVolumes dHT_y0(N=N) annotation (Placement(transformation(extent={{24,72},{44,92}}), iconTransformation(extent={{24,84},{44,92}})));
  ThermoPower.Thermal.DHTVolumes dHT_y1(N=N) annotation (Placement(transformation(extent={{72,-80},{92,-60}}), iconTransformation(extent={{16,-94},{36,-86}})));
  ThermoPower.Thermal.HT hT_x0 annotation (Placement(transformation(extent={{-106,-10},{-86,10}}), iconTransformation(extent={{-106,-10},{-86,10}})));
  ThermoPower.Thermal.HT hT_xN annotation (Placement(transformation(extent={{88,-10},{108,10}}), iconTransformation(extent={{88,-10},{108,10}})));

  OpenTEMPEST.SOC.Cell.Cell1D.AirChannel airChannel(
    N=N,
    redeclare function fluxInterp = fluxInterp,
    TStartIn=TStart_ACin,
    TStartOut=TStart_ACout,
    pStartIn=pStartAC,
    lX=lX,
    lY=lYac,
    lZ=lZac,
    xStart=xStartAC,
    por=porAC,
    Nu_PEN=Nu_PENac,
    Nu_IC=Nu_ICac,
    pDrop=pDropAC)
    annotation (Placement(transformation(extent={{-8,-32},{12,-52}})));

  OpenTEMPEST.SOC.Cell.Cell1D.FuelChannel fuelChannel(
    N=N,
    redeclare function fluxInterp = fluxInterp,
    TStartIn=TStart_FCin,
    TStartOut=TStart_FCout,
    pStartIn=pStartInFC,
    pStartOut=pStartOutFC,
    xStart=xStartFC,
    lX=lX,
    lY=lYfc,
    lZ=lZfc,
    por=porFC,
    Nu_PEN=Nu_PENfc,
    Nu_IC=Nu_ICfc,
    k_eff_Ni=kFoam,
    cp_Ni=cpFoam,
    rho_Ni=rhoFoam,
    dpNom=dpNomFC)
    annotation (Placement(transformation(extent={{-10,28},{10,48}})));

  OpenTEMPEST.SOC.Cell.Cell1D.PEN1D pen(
    redeclare package SolidMat = OpenTEMPEST.Solid.Material.Custom,
    N=N,
    Tstartbar=TStart_pen,
    kCustom_trans=kPen,
    kCustom_long=kPen,
    rhoCustom=rhoPen,
    cpCustom=cpPen,
    Jstart=Jstart,
    lX=lXpen,
    lY=lYpen,
    lZ=lZpen,
    redeclare model Electrochem = Electrochem)
    annotation (Placement(transformation(extent={{-28,-26},{22,20}})));

  OpenTEMPEST.Heat.RadHTFV radHTFV(
    N=N,
    specialCase=OpenTEMPEST.Enumerations.RadiationSpecialCases.case1,
    epsilon1=epsilon_pen,
    epsilon2=epsilon_ic,
    A1=lX*lYac,
    A2=lX*lYac)
    annotation (Placement(transformation(extent={{14,-44},{34,-24}})));

  ThermoPower.Thermal.ConvHTFV convHTFV(Nv=N, G=LambdaIC/lZac*(lY - lYac)*lX) "why is this here?"
                                                                              annotation (Placement(transformation(extent={{-52,-42},
            {-32,-22}})));

  OpenTEMPEST.Heat.Solid1D ICAir(
    redeclare package SolidMat = OpenTEMPEST.Solid.Material.Custom,
    N=N,
    Tstartbar=TStart_IC,
    lX=lX,
    lY=lY,
    lZ=lZIC,
    kCustom_trans=LambdaIC,
    rhoCustom=rhoIC,
    cpCustom=cpIC)
    annotation (Placement(transformation(extent={{8,-50},{28,-70}})));

  OpenTEMPEST.Heat.Solid1D ICFuel(
    redeclare package SolidMat = OpenTEMPEST.Solid.Material.Custom,
    N=N,
    Tstartbar=TStart_IC,
    lX=lX,
    lY=lY,
    lZ=lZIC,
    kCustom_trans=LambdaIC,
    rhoCustom=rhoIC,
    cpCustom=cpIC)
    annotation (Placement(transformation(extent={{-24,78},{-4,58}})));

  ThermoPower.Gas.PressDropLin pressDropLin1(
    redeclare package Medium = FCMedium,
    final R=HyR_Fuel,
    final allowFlowReversal=false) annotation (Placement(transformation(extent={{40,28},{60,48}})));
  ThermoPower.Gas.PressDropLin pressDropLin3(
    redeclare package Medium = ACMedium,
    final R=HyR_Air,
    final allowFlowReversal=false) annotation (Placement(transformation(extent={{58,-58},{78,-38}})));
equation
  connect(pen.pin_n, pin_n) annotation (Line(points={{-28.25,5.97},{-82,5.97},{-82,-22},{-96,-22}}, color={0,0,255}));
  connect(pen.pin_p, pin_p) annotation (Line(points={{-28.25,-12.43},{-78,-12.43},{-78,-44},{-96,-44}},
                                                                                                      color={0,0,255}));
  connect(pen.PEN_in, fuelChannel.PEN_in) annotation (Line(points={{-3,13.1},{-3,28},{-5.2,28},{-5.2,34.8}}, color={0,0,0}));
  connect(airChannel.PEN_in, pen.PEN_ina) annotation (Line(points={{-3.2,-38.8},{-3.2,-34},{-3,-34},{-3,-19.1}}, color={0,0,0}));
  connect(airChannel.Q_IC, ICAir.dhT_z0) annotation (Line(points={{2,-46},{2,-50},{11,-50},{11,-53}}, color={255,127,0}));
  connect(ICAir.dhT_z1, dHT_z1) annotation (Line(points={{11,-67},{11,-70},{-50,-70}}, color={255,127,0}));
  connect(dHT_z0, ICFuel.dhT_z0) annotation (Line(points={{-50,82},{-21,82},{-21,75}}, color={255,127,0}));
  connect(fuelChannel.Q_IC, ICFuel.dhT_z1) annotation (Line(points={{0,42},{2,42},{2,52},{-21,52},{-21,61}}, color={255,127,0}));
  connect(convHTFV.side2, ICAir.dhT_z0) annotation (Line(points={{-42,-35.1},{
          -42,-54},{11,-54},{11,-53}},                                                                                     color={255,127,0}));
  connect(airChannel.outfl, pressDropLin3.inlet) annotation (Line(points={{11.4,-42},{36,-42},{36,-48},{58,-48}}, color={159,159,223}));
  connect(pressDropLin3.outlet, airFlangeOut) annotation (Line(points={{78,-48},{100,-48}}, color={159,159,223}));
  connect(fuelChannel.outfl, pressDropLin1.inlet) annotation (Line(points={{9.4,38},{40,38}}, color={159,159,223}));
  connect(pressDropLin1.outlet, fuelFlangeOut) annotation (Line(points={{60,38},{66,38},{66,54},{100,54}}, color={159,159,223}));
  connect(fuelFlangeIn, fuelChannel.infl) annotation (Line(points={{-100,54},{-80,54},{-80,44},{-46,44},{-46,37.8},{-9.4,37.8}}, color={159,159,223}));
  connect(airFlangeIn, airChannel.infl) annotation (Line(points={{-100,-66},{-84,-66},{-84,-50},{-22,-50},{-22,-41.8},{-7.4,-41.8}}, color={159,159,223}));
  connect(pen.dhT_z1, fuelChannel.Q_PEN) annotation (Line(points={{-20.5,13.1},{-20.5,26},{-14,26},{-14,24},{-0.1,24},{-0.1,34.6}}, color={255,127,0}));
  connect(pen.dhT_z0, airChannel.Q_PEN) annotation (Line(points={{-20.5,-19.1},{-20.5,-38.6},{1.9,-38.6}}, color={255,127,0}));
  connect(pen.dhT_y1, dHT_y1) annotation (Line(points={{14.5,13.1},{116,13.1},{116,-70},{82,-70}}, color={255,127,0}));
  connect(pen.dhT_y0, dHT_y0) annotation (Line(points={{14.5,-19.1},{34,-19.1},{34,82}}, color={255,127,0}));
  connect(radHTFV.side1, pen.Qrad_AI) annotation (Line(points={{24,-29.6},{22,-29.6},{22,-7.6},{11,-7.6}}, color={255,127,0}));
  connect(ICAir.dhT_int, radHTFV.side2) annotation (Line(points={{18,-60.4},{18,-34.4},{24,-34.4}}, color={255,127,0}));
  connect(pen.hT_xN, hT_xN) annotation (Line(points={{24.5,-3},{84,-3},{84,0},{98,0}}, color={191,0,0}));
  connect(hT_x0, pen.hT_x0) annotation (Line(points={{-96,0},{-36,0},{-36,-3},{-30.5,-3}}, color={191,0,0}));
  connect(convHTFV.side1, pen.dhT_z0) annotation (Line(points={{-42,-29},{-42,
          -19.1},{-20.5,-19.1}}, color={255,127,0}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,86},{100,-92}},
          lineColor={0,0,0},
          fillColor={0,127,127},
          fillPattern=FillPattern.CrossDiag)}),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<h2>1D Cell</h2>

<p>
The model represents a one-dimensional discretized 
into N control volumes along the channel direction.
</p>

<p>
The cell consists of a discretized PEN with electrochemistry, axial fuel and air channels, 
and metallic interconnects. Heat transfer between PEN, channels, and interconnects 
is modeled via convective and radiative elements. Gas manifolds supply the 
reactants, while Darcy-based pressure drop and electrochemical reactions are 
computed along the channels. Geometrical dimensions, material properties, 
discretization, and transport parameters are fully parameterized.
</p>
</html>"));
end Cell;

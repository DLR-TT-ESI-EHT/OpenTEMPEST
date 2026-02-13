within OpenTEMPEST.SOC.Stack;
model Stack1D
  "detailed stack based on Sunfire cell. Each cell is fully calculated"
  import SI = Modelica.SIunits;

  parameter Integer N(min=3)=3;
  parameter Integer Ncell(min=2)=5;
  replaceable function fluxInterp =
      Flow.FluxInterpolators.UDSinterp              constrainedby Flow.FluxInterpolators.DifferencingSchemeInterpBase
                                                                                         annotation(choicesAllMatching = true);

  // Initialisation
  parameter SI.PressureDifference dpNomFC=0.01 "nominal pressure loss for initialization" annotation (Dialog(tab="Initialization"));
  parameter SI.Temperature TStart=1023.15
                                   "Starting temperature of cell" annotation (Dialog(tab = "Initialization"));
  parameter SI.CurrentDensity JStart=0
                                     "Starting current density in PEN" annotation (Dialog(tab="Initialization"));
  parameter SI.AbsolutePressure pStartInFC=pStartOutFC+dpNomFC "Starting pressure in FC" annotation (Dialog(tab = "Initialization"));
  parameter SI.AbsolutePressure pStartOutFC=140000 "Starting pressure in FC" annotation (Dialog(tab = "Initialization"));
  parameter SI.MassFraction xStartFC[FCMedium.nXi]= FCMedium.reference_X "Starting Mass Fractions in FC" annotation (Dialog(tab = "Initialization"));
  parameter SI.AbsolutePressure pStartAC=140000 "Starting Pressure in AC" annotation (Dialog(tab = "Initialization"));
  parameter SI.MassFraction xStartAC[ACMedium.nXi] = ACMedium.reference_X "Starting composition in AC" annotation (Dialog(tab = "Initialization"));

  parameter SI.Temperature TStart_pen=TStart "Starting Temperature of PEN" annotation (Dialog(tab="Initialization"));
  parameter SI.Temperature TStart_FCin=TStart "Starting temperature at FC inlet" annotation (Dialog(tab = "Initialization"));
  parameter SI.Temperature TStart_FCout=TStart "Starting temperature at FC outlet" annotation (Dialog(tab = "Initialization"));
  parameter SI.Temperature TStart_ACin=TStart "Starting temperature at AC inlet" annotation (Dialog(tab = "Initialization"));
  parameter SI.Temperature TStart_ACout=TStart "Starting Temperature at AC outlet " annotation (Dialog(tab = "Initialization"));
  parameter SI.Temperature TStart_IC = TStart "Starting temperature of Interconnect" annotation (Dialog(tab = "Initialization"));
  parameter SI.Temperature TStart_TopPlate = TStart "starting Temperature in Top Plate" annotation (Dialog(tab = "Initialization"));
  parameter SI.Temperature TStart_BottomPlate= TStart "starting Temperature in Bottom Plate" annotation (Dialog(tab = "Initialization"));
  parameter SI.Temperature TStart_IntermediatePlate= TStart "starting Temperature in Bottom Plate" annotation (Dialog(tab = "Initialization"));

  //Dimensions
  parameter SI.Length lY=0.142 "Width of FC"
    annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length lYfc=lY "Width of FC"
    annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length lZfc=1.23e-3 "Height of FC"
    annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length lX=0.09 "Length of PEN"
    annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length lXac=lX "Length of AC"
    annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length lYac=0.08448 "Width of AC"
    annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length lZac=1.2e-3 "Height of AC"
    annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length lYpen=lY "Width of PEN"
    annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length lZpen=150e-6 "Thickness of PEN including GDC layer"
    annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length lZIC=250e-6 "Thickness of IC over the active part of cell"
    annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length lZTopPlate=11e-3 "Top Plate thickness" annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length lZBottomPlate=11e-3 "Bottom Plate thickness" annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length lZIntermediatePlate=11e-3 "Bottom Plate thickness" annotation (Dialog(tab="Cell Dimensions"));

  // PEN parameters
  replaceable model Electrochem =
      OpenTEMPEST.SOC.Electrochem.Components.BV_Steam                                     constrainedby
    OpenTEMPEST.SOC.Electrochem.Components.ElectrochemBase
    annotation (Placement(transformation(extent={{78,50},{98,70}})), choicesAllMatching=true);
  parameter SI.ThermalConductivity kPen=2 "Thermal conductivity of PEN in W/mK" annotation (Dialog(tab="PEN"));
  parameter SI.Density rhoPen=5900 "Density of PEN in kg/m3" annotation (Dialog(tab="PEN"));
  parameter SI.SpecificHeatCapacity cpPen=500 "Specific heat capacity of PEN in J/kgK" annotation (Dialog(tab="PEN"));
  parameter SI.SpectralEmissivity epsilon_pen=0.8 "emissivity of Anode-Electrolyte-Cathode unit" annotation (Dialog(tab="PEN"));
  // FC parameters
  replaceable package FCMedium = Medium.Fuel_CH4          annotation (Dialog(tab="Fuel Channel"));
  parameter Real porFC=0.867 "Porosity in FC" annotation (Dialog(tab="Fuel Channel"));
  parameter Real Nu_PENfc=12 "Nusselt number fuel channel on PEN side" annotation (Dialog(tab="Fuel Channel"));
  parameter Real Nu_ICfc=9.86 "Nusselt Number fuel channel on IC side" annotation (Dialog(tab="Fuel Channel"));
  parameter SI.ThermalConductivity kFoam=3.576 "Thermal Conductivity of Ni Foam" annotation (Dialog(tab="Fuel Channel"));
  parameter SI.SpecificHeatCapacity cpFoam=440       "Specific heat capacity of Ni in FC" annotation (Dialog(tab="Fuel Channel"));
  parameter SI.Density rhoFoam(displayUnit="kg/m3")=8908 "Density of Ni in FC" annotation (Dialog(tab="Fuel Channel"));
  parameter ThermoPower.Units.HydraulicResistance HyR_Fuel=48991852 "Hydraulic resistance of fuel channel of each single cell in stack"
    annotation (Dialog(tab="Fuel Channel"));

  // AC parameters
  replaceable package ACMedium = Medium.Air_Medium         annotation (Dialog(tab="Air Channel"));
  parameter Real porAC=1 "Porosity in AC" annotation (Dialog(tab="Air Channel"));
  parameter Real pDropAC=0 "Pressure drop factor, pOut=(1-pDropAC)*pIn" annotation (Dialog(tab="Air Channel"));
  parameter Real Nu_PENac=8.235 "Nusselt number air channel on PEN side" annotation (Dialog(tab="Air Channel"));
  parameter Real Nu_ICac=7.54 "Nusselt Number air channel on IC side" annotation (Dialog(tab="Air Channel"));
  parameter ThermoPower.Units.HydraulicResistance HyR_Air=5593470 "Hydraulic resistance of air channel of each single cell in stack" annotation (Dialog(tab="Air Channel"));

  // Interconnect parameters
  parameter SI.ThermalConductivity k_solidCustom=0.2812 "Thermal Conductivity of solid parallel to windows" annotation (Dialog(tab="Interconnect"));
  parameter SI.SpecificHeatCapacity cp_solidCustom=463.8 "Specific heat capacity of solid parallel to windows" annotation (Dialog(tab="Interconnect"));
  parameter SI.Density rho_solidCustom(displayUnit="kg/m3")=1330 "Density of solid parallel to windows" annotation (Dialog(tab="Interconnect"));
  parameter SI.SpectralEmissivity epsilon_ic=0.1 "emissivity of surface of IC" annotation (Dialog(tab = "Interconnect"));

  //End / Intermediate plates parameters
  parameter Integer nIMPs = integer((floor((Ncell - 1)/10))) "Number of intermediate plates";
  parameter SI.ThermalConductivity kTopPlate=k_solidCustom "Thermal Conductivity of solid parallel to windows" annotation (Dialog(tab="End / Intermediate plates"));
  parameter SI.SpecificHeatCapacity cpTopPlate=cp_solidCustom "Specific heat capacity of solid parallel to windows" annotation (Dialog(tab="End / Intermediate plates"));
  parameter SI.Density rhoTopPlate(displayUnit="kg/m3")=rho_solidCustom "Density of solid parallel to windows" annotation (Dialog(tab="End / Intermediate plates"));
  parameter SI.ThermalConductivity kBottomPlate=k_solidCustom "Thermal Conductivity of solid parallel to windows" annotation (Dialog(tab="End / Intermediate plates"));
  parameter SI.SpecificHeatCapacity cpBottomPlate=cp_solidCustom "Specific heat capacity of solid parallel to windows" annotation (Dialog(tab="End / Intermediate plates"));
  parameter SI.Density rhoBottomPlate(displayUnit="kg/m3")=rho_solidCustom "Density of solid parallel to windows" annotation (Dialog(tab="End / Intermediate plates"));
  parameter SI.ThermalConductivity kIntermediatePlate=k_solidCustom "Thermal Conductivity of solid parallel to windows" annotation (Dialog(tab="End / Intermediate plates"));
  parameter SI.SpecificHeatCapacity cpIntermediatePlate=cp_solidCustom "Specific heat capacity of solid parallel to windows" annotation (Dialog(tab="End / Intermediate plates"));
  parameter SI.Density rhoIntermediatePlate(displayUnit="kg/m3")=rho_solidCustom "Density of solid parallel to windows" annotation (Dialog(tab="End / Intermediate plates"));

  Flow.Manifold Fuel_Manifold(redeclare package Medium = Medium.Fuel_CH4,
      nPorts_b=Ncell)
    annotation (Placement(transformation(extent={{-82,56},{-70,88}})));
  Flow.Manifold Air_Manifold(redeclare package Medium = Medium.Air_Medium,
      nPorts_b=Ncell)
    annotation (Placement(transformation(extent={{-90,-70},{-78,-38}})));
  ThermoPower.Gas.FlangeB fuelOut(redeclare package Medium = Medium.Fuel_CH4)
    annotation (Placement(transformation(extent={{100,50},{120,70}}),
        iconTransformation(extent={{100,50},{120,70}})));
  ThermoPower.Gas.FlangeB airOut(redeclare package Medium = Medium.Air_Medium)
    annotation (Placement(transformation(extent={{110,-72},{130,-52}}),
        iconTransformation(extent={{110,-72},{130,-52}})));
  Flow.Manifold_out manifold_outFuel(redeclare package Medium = Medium.Fuel_CH4,
      final nPorts_a=Ncell)
    annotation (Placement(transformation(extent={{70,34},{90,54}})));
  Flow.Manifold_out manifold_outAir(redeclare package Medium =
        Medium.Air_Medium, final nPorts_a=Ncell)
    annotation (Placement(transformation(extent={{66,-20},{86,0}})));
  ThermoPower.Gas.FlangeA fuelIn(redeclare package Medium = Medium.Fuel_CH4)
    annotation (Placement(transformation(extent={{-124,50},{-104,70}}),
        iconTransformation(extent={{-124,50},{-104,70}})));
  ThermoPower.Gas.FlangeA airIn(redeclare package Medium = Medium.Air_Medium)
    annotation (Placement(transformation(extent={{-124,-92},{-104,-72}}),
        iconTransformation(extent={{-124,-92},{-104,-72}})));

  Modelica.Electrical.Analog.Interfaces.PositivePin pin_p
    annotation (Placement(transformation(extent={{-132,-58},{-112,-38}}),
        iconTransformation(extent={{-132,-58},{-112,-38}})));
  Modelica.Electrical.Analog.Interfaces.NegativePin pin_n
    annotation (Placement(transformation(extent={{-132,-32},{-112,-12}}),
        iconTransformation(extent={{-132,-32},{-112,-12}})));

  ThermoPower.Thermal.DHTVolumes dHT_top(N=N)       annotation (Placement(
        transformation(extent={{-16,68},{4,88}}), iconTransformation(extent={{-50,
            30},{44,48}})));
  ThermoPower.Thermal.DHTVolumes dHT_bottom(N=N)       annotation (Placement(
        transformation(extent={{-16,-76},{4,-56}}), iconTransformation(extent={{
            -50,-52},{44,-34}})));

  ThermoPower.Thermal.DHTVolumes dHT_y0(N=N)       annotation (Placement(
        transformation(extent={{-32,100},{26,110}}), iconTransformation(extent={{-60,90},
            {60,100}})));

  ThermoPower.Thermal.DHTVolumes dHT_y1(N=N)       annotation (Placement(
        transformation(extent={{-30,-110},{28,-100}}), iconTransformation(
          extent={{-60,-100},{60,-90}})));
  ThermoPower.Thermal.DHTVolumes dHT_x0(N=Ncell)
    annotation (Placement(transformation(extent={{-136,8},{-116,28}}),
        iconTransformation(extent={{-124,0},{-112,40}})));
  ThermoPower.Thermal.DHTVolumes dHT_xN(N=Ncell)
    annotation (Placement(transformation(extent={{110,2},{130,22}}),
        iconTransformation(extent={{108,0},{120,40}})));
  ThermoPower.Thermal.HT hTPlates_x0[2 + nIMPs] annotation (Placement(transformation(extent={{-124,80},
            {-104,100}}),                                                                                           iconTransformation(extent={{-124,80},
            {-104,100}})));
  ThermoPower.Thermal.HT hTPlates_xN[2 + nIMPs] annotation (Placement(transformation(extent={{108,24},{128,44}}), iconTransformation(extent={{100,80},
            {120,100}})));
  SI.EnergyFlowRate EfSum "check energy balance through sum of energy flows";
  SI.MassFlowRate massSum "check mass balance through sum of mass flows";
  Cell.Cell1D.Cell cell[Ncell](
    each N=N,
    redeclare function fluxInterp = fluxInterp,
    each TStart=TStart,
    each Jstart=JStart,
    each pStartInFC=pStartInFC,
    each pStartOutFC=pStartOutFC,
    each pStartAC=pStartAC,
    each xStartAC=xStartAC,
    each TStart_pen=TStart_pen,
    each TStart_FCin=TStart_FCin,
    each TStart_FCout=TStart_FCout,
    each TStart_ACin=TStart_ACin,
    each TStart_ACout=TStart_ACout,
    each TStart_IC=TStart_IC,
    redeclare model Electrochem = Electrochem,
    each epsilon_pen=epsilon_pen,
    each lYfc=lYfc,
    each lZfc=lZfc,
    each lX=lX,
    each lY=lY,
    each Nu_PENfc=Nu_PENfc,
    each Nu_ICfc=Nu_ICfc,
    each HyR_Fuel=HyR_Fuel,
    each lXac=lXac,
    each lYac=lYac,
    each lZac=lZac,
    each lYpen=lYpen,
    each lZpen=lZpen,
    each Nu_PENac=Nu_PENac,
    each Nu_ICac=Nu_ICac,
    each HyR_Air=HyR_Air,
    each lZIC=lZIC,
    each kPen=kPen,
    each xStartFC=xStartFC,
    each rhoPen=rhoPen,
    each cpPen=cpPen,
    redeclare package FCMedium = FCMedium,
    redeclare package ACMedium = ACMedium,
    each porFC=porFC,
    each kFoam=kFoam,
    each cpFoam=cpFoam,
    each rhoFoam=rhoFoam,
    each porAC=porAC,
    each pDropAC=pDropAC,
    each LambdaIC=k_solidCustom,
    each cpIC=cp_solidCustom,
    each rhoIC=rho_solidCustom,
    each epsilon_ic=epsilon_ic,
    each dpNomFC=dpNomFC)
    annotation (Placement(transformation(extent={{-20,-14},{16,22}})));

  OpenTEMPEST.Heat.Solid1D plates[2 + nIMPs](
    redeclare package SolidMat = TEMPEST.Solid.Material.Custom,
    each N=N,
    each TstartX0=TStart_IntermediatePlate,
    each TstartXN=TStart_IntermediatePlate,
    each lX=lX,
    each lY=lY,
    lZ=cat(
        1,
        {lZBottomPlate},
        fill(lZIntermediatePlate, nIMPs),
        {lZTopPlate}),
    each kCustom_trans=kIntermediatePlate,
    each rhoCustom=rhoIntermediatePlate,
    each cpCustom=cpIntermediatePlate)
    annotation (Placement(transformation(extent={{24,36},{44,56}})));

protected
  Heat.HTs_DHT hTs_DHT_x0(each N=Ncell)
    annotation (Placement(transformation(extent={{-60,10},{-72,22}})));
  Heat.HTs_DHT hTs_DHT_xN(each N=Ncell)
    annotation (Placement(transformation(extent={{74,2},{86,14}})));

initial equation
  //Fuel_Manifold.ports_b.m_flow = fill(fuelIn.m_flow/Ncell,Ncell);
equation
  EfSum = fuelIn.m_flow*inStream(fuelIn.h_outflow) + fuelOut.m_flow*fuelOut.h_outflow
              + airIn.m_flow*inStream(airIn.h_outflow) + airOut.m_flow*airOut.h_outflow
              + (pin_p.v-pin_n.v)*pin_n.i
              + sum(dHT_x0.Q[:]) + sum(dHT_xN.Q[:])
              + sum(dHT_top.Q[:]) + sum(dHT_bottom.Q[:])
              + sum(dHT_y0.Q[:]) + sum(dHT_y1.Q[:])
              + sum(hTPlates_x0[:].Q_flow)
              + sum(hTPlates_xN[:].Q_flow); // useful for steady state, but not subtracting all heat capacitors, as in ECReactorModels.Stack.SimplifiedStack.SimplifiedStackVertMult

  massSum = fuelIn.m_flow + fuelOut.m_flow + airIn.m_flow + airOut.m_flow;

  // Electrical connections - External
  connect(pin_n, cell[1].pin_n) annotation (Line(points={{-122,-22},{-46,-22},{-46,
          0.04},{-19.28,0.04}},                                                                        color={0,0,255}));
  connect(pin_p, cell[Ncell].pin_p) annotation (Line(points={{-122,-48},{-38,-48},
          {-38,-3.92},{-19.28,-3.92}},                                                                         color={0,0,255}));
  //Z-direction Heat ports
  // z-direction HT cell<->endplates
  connect(cell[1].dHT_z0, plates[1].dhT_z1) annotation (Line(points={{-11,-8.6},{-11,-22},{-36,-22},{-36,62},{27,62},{27,53}}, color={255,127,0}));
  connect(cell[end].dHT_z1, plates[end].dhT_z0) annotation (Line(points={{-11.36,17.32},{-11.36,26},{27,26},{27,39}}, color={255,127,0}));
  //Z-direction HT external<->endplates
  connect(plates[end].dhT_z1, dHT_top) annotation (Line(points={{27,53},{27,58},{-6,58},{-6,78}}, color={255,127,0}));
  connect(plates[1].dhT_z0, dHT_bottom) annotation (Line(points={{27,39},{27,-38},{-8,-38},{-8,-66},{-6,-66}}, color={255,127,0}));
  //Inter cell connections
  for i in 1:Ncell - 1 loop
    if mod(i, 10) == 0 and i > 10 then
      //z-direction heat transfer intermediate plate<->cell
      connect(cell[i + 1].dHT_z0, plates[integer(i/10)].dhT_z1) annotation (Line(points={{-11,-8.6},{-11,-22},{-36,-22},{-36,62},{27,62},{27,53}}, color={255,127,0}));
      connect(cell[i].dHT_z1, plates[integer(i/10)].dhT_z0) annotation (Line(points={{-11.36,17.32},{-11.36,26},{27,26},{27,39}}, color={255,127,0}));
    else
      // z-direction HT cell<->cell
      connect(cell[i].dHT_z1, cell[i + 1].dHT_z0);
    end if;
    // Series/Internal Electrical Connections
    connect(cell[i].pin_p, cell[i + 1].pin_n);
  end for;

  // Gas flow connections
  //Fuel In
  connect(fuelIn, Fuel_Manifold.port_a) annotation (Line(points={{-114,60},{-87,
          60},{-87,72},{-82,72}},                                                                       color={159,159,223}));
  connect(Fuel_Manifold.ports_b[1:Ncell], cell.fuelFlangeIn) annotation (Line(points={{-70,72},{-32,72},{-32,13.72},{-20,13.72}}, color={0,127,255}));
  //AirIn
  connect(airIn, Air_Manifold.port_a) annotation (Line(points={{-114,-82},{-92,-82},
          {-92,-54},{-90,-54}},                                                                           color={159,159,223}));
  connect(Air_Manifold.ports_b, cell.airFlangeIn) annotation (Line(points={{-78,-54},{-60,-54},{-60,-58},{-32,-58},{-32,-7.88},{-20,-7.88}}, color={0,127,255}));
  //Fuel out
  connect(manifold_outFuel.ports_a, cell.fuelFlangeOut) annotation (Line(points={{77.4,44},{66,44},{66,13.72},{16,13.72}}, color={0,127,255}));
  connect(manifold_outFuel.port_b, fuelOut) annotation (Line(points={{85,44},{94,
          44},{94,60},{110,60}},                                                                        color={0,127,255}));
  //Air out
  connect(manifold_outAir.port_b, airOut) annotation (Line(points={{81,-10},{92,
          -10},{92,-62},{120,-62}},                                                                       color={0,127,255}));
  connect(manifold_outAir.ports_a, cell.airFlangeOut) annotation (Line(points={{73.4,-10},{70,-10},{70,-4.64},{16,-4.64}}, color={0,127,255}));

  // x-direction Heat ports
  // x-direction dHT<>hTs_DHT<>cells
  connect(hTs_DHT_x0.DHT_port, dHT_x0) annotation (Line(points={{-72.6,16},{-126,16},{-126,18}}, color={255,127,0}));
  connect(hTs_DHT_xN.DHT_port, dHT_xN) annotation (Line(points={{86.6,8},{98,8},{98,12},{120,12}}, color={255,127,0}));
  connect(hTs_DHT_x0.HT_ports, cell.hT_x0) annotation (Line(points={{-58.8,16.12},{-44,16.12},{-44,4},{-19.28,4}}, color={191,0,0}));
  connect(cell.hT_xN, hTs_DHT_xN.HT_ports) annotation (Line(points={{15.64,4},{60,4},{60,8.12},{72.8,8.12}}, color={191,0,0}));

  //X-direction intermediate + end plates to external
  for i in 1:(nIMPs+2) loop
    connect(plates.hT_xN, hTPlates_xN) annotation (Line(points={{45,46},{72,46},{72,34},{118,34}}, color={191,0,0}));
    connect(plates.hT_x0, hTPlates_x0) annotation (Line(points={{23,46},{-92,46},
            {-92,90},{-114,90}},                                                                      color={191,0,0}));
  end for;

  // y-direction ports
  //Y-direction cells to external
   for i in 1:Ncell loop
     connect(cell[i].dHT_y0, dHT_y0) annotation (Line(points={{4.12,19.84},{4.12,30},{16,30},{16,92},{-3,92},{-3,105}}, color={255,127,0}));
     connect(cell[i].dHT_y1, dHT_y1) annotation (Line(points={{2.68,-12.2},{2.68,-20},{24,-20},{24,-90},{0,-90},{0,-105},{-1,-105}}, color={255,127,0}));
   end for;

  //Y-direction intermediate and end plates to external
  for i in 1:(nIMPs+2) loop
    connect(plates[i].dhT_y0, dHT_y0);
    connect(plates[i].dhT_y1, dHT_y1);
  end for;

  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Rectangle(
          extent={{-124,100},{120,-100}},
          lineColor={28,108,200},
          fillColor={28,108,200},
          fillPattern=FillPattern.CrossDiag),
        Rectangle(
          extent={{-102,80},{102,-82}},
          lineColor={28,108,200},
          fillColor={0,0,99},
          fillPattern=FillPattern.Forward),
        Text(
          extent={{-98,24},{98,-28}},
          lineColor={255,255,255},
          fillColor={28,108,200},
          fillPattern=FillPattern.CrossDiag,
          textStyle={TextStyle.Bold},
          textString="Detailed Stack")}),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(StopTime=50000, __Dymola_Algorithm="Cvode"),
    __Dymola_experimentSetupOutput(equidistant=false),
    __Dymola_experimentFlags(
      Advanced(GenerateVariableDependencies=false, OutputModelicaCode=true),
      Evaluate=false,
      OutputCPUtime=true,
      OutputFlatModelica=false),
    Documentation(revisions="<html>
<ul>
<li><i>23 Dez 2021</i> by Faisal Sedeqi</a>:<br>First release. </li>
</ul>
</html>"));
end Stack1D;

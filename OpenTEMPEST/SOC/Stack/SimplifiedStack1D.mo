within OpenTEMPEST.SOC.Stack;
model SimplifiedStack1D
  "Stack based on Sunfire cells with simplification approach. From a total of Ncells, Ndetailed cells are calculated fully and the rest approximated through simplified cell heat transfer blocks, with flow and voltage interpolation."
  import SI = Modelica.SIunits;
  parameter Integer N(min=3)=3 annotation (Dialog(tab = "General"));
  parameter Integer Ncell(min=2)=30 annotation (Dialog(tab = "General"));
  parameter Integer[:] verticalBlockSize={3, 3, 2, 2, 2, 1, 1, 1, 1, 2, 2, 3, 3} "Vertical size of the simplified blocks (e. g. 3 represents three cells merged in one simplfied block). Sum needs to be Ncell-NdetailedCell" annotation (Dialog(tab = "General"));
  parameter Integer NdetailedCell(min=1, max=Ncell) = 4 "number of cells that are calculated in detail" annotation (Dialog(tab = "General"));
  replaceable function fluxInterp =
      Flow.FluxInterpolators.UDSinterp              constrainedby Flow.FluxInterpolators.DifferencingSchemeInterpBase
                                                                                         annotation(choicesAllMatching = true, Dialog(tab="General"));

  // Initialisation
  parameter SI.Temperature TStart=1000  "Starting Temperature of Cell" annotation (Dialog(tab = "Initialization"));
  parameter SI.CurrentDensity JStart=0 "Starting Current Density in PEN" annotation (Dialog(tab="Initialization"));
  parameter SI.AbsolutePressure pStart=110000  "default cell Starting Pressure" annotation (Dialog(tab="Initialization"));
  parameter SI.AbsolutePressure pStartInFC=pStart "Starting pressure in FC" annotation (Dialog(tab = "Initialization"));
  parameter SI.AbsolutePressure pStartOutFC=pStart "Starting pressure in FC" annotation (Dialog(tab = "Initialization"));
  parameter SI.MassFraction xStartFC[FCMedium.nXi]= FCMedium.reference_X "Starting Mass Fractions in FC" annotation (Dialog(tab = "Initialization"));
  parameter SI.AbsolutePressure pStartAC=pStart "Starting Pressure in AC" annotation (Dialog(tab = "Initialization"));
  parameter SI.MassFraction xStartAC[ACMedium.nXi] = ACMedium.reference_X "Starting composition in AC" annotation (Dialog(tab = "Initialization"));

  parameter SI.Temperature TStart_pen=TStart "Starting Temperature of PEN" annotation (Dialog(tab="Initialization"));
  parameter SI.Temperature TStart_FCin=TStart "Starting temperature at FC inlet" annotation (Dialog(tab = "Initialization"));
  parameter SI.Temperature TStart_FCout=TStart "Starting temperature at FC outlet" annotation (Dialog(tab = "Initialization"));
  parameter SI.Temperature TStart_ACin=TStart "Starting temperature at AC inlet" annotation (Dialog(tab = "Initialization"));
  parameter SI.Temperature TStart_ACout=TStart "Starting Temperature at AC outlet " annotation (Dialog(tab = "Initialization"));
  parameter SI.Temperature TStart_IC = TStart "Starting temperature of Interconnect" annotation (Dialog(tab = "Initialization"));
  parameter SI.Temperature TStart_FCvol = TStart "Starting temperature of FC in passive volume" annotation (Dialog(tab = "Initialization"));
  parameter SI.Temperature TStart_ICvol = TStart "Starting temperature of IC in passive volume" annotation (Dialog(tab = "Initialization"));
  parameter SI.Temperature TStart_ACvol = TStart "Starting temperature of AC in passive volume" annotation (Dialog(tab = "Initialization"));
  parameter SI.Temperature TStart_plates = TStart "starting Temperature in Top Plate" annotation (Dialog(tab = "Initialization"));
  parameter SI.PressureDifference dpNomFC=0.01 "nominal pressure loss for initialization" annotation (Dialog(tab="Initialization"));

  //Dimensions
  parameter SI.Length lY=0.142 "Width of FC" annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length lYfc=lY "Width of FC" annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length lZfc=1.23e-3 "Height of FC" annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length lX=0.09 "Length of PEN" annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length lXac=lX "Length of AC" annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length lYac=0.08448 "Width of AC" annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length lZac=1.2e-3 "Height of AC" annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length lYpen=lY "Width of PEN" annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length lZpen=150e-6 "Thickness of PEN including GDC layer" annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length lZIC=2*250e-6 "Thickness of IC over the active part of cell" annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length lZTopPlate=9.3e-3
                                       "Top Plate thickness" annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length lZBottomPlate=9.3e-3
                                          "Bottom Plate thickness" annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length lZIntermediatePlate=2.5e-3
                                                "Bottom Plate thickness" annotation (Dialog(tab="Cell Dimensions"));
  parameter SI.Length tau_ae=30e-6 " thickness of fuel electrode" annotation (Dialog(tab="Electrode-Electrolyte parameters"));
  parameter SI.Length tau_ce=55e-6 "thickness of air electrode" annotation (Dialog(tab="Electrode-Electrolyte parameters"));
  parameter SI.Length tau_se=90e-6 "thickness of electrolyte" annotation (Dialog(tab="Electrode-Electrolyte parameters"));

  // PEN parameters
  replaceable model Electrochem =
      OpenTEMPEST.SOC.Electrochem.Components.BV_Steam                                  constrainedby
    OpenTEMPEST.SOC.Electrochem.Components.ElectrochemBase                                                                                                               annotation (Placement(transformation(extent={{78,50},{98,70}})),choicesAllMatching=true);
  parameter SI.ThermalConductivity kPen=2 "Thermal conductivity of PEN in W/mK" annotation (Dialog(tab="PEN"));
  parameter SI.Density rhoPen=5900 "Density of PEN in kg/m3" annotation (Dialog(tab="PEN"));
  parameter SI.SpecificHeatCapacity cpPen=500 "Specific heat capacity of PEN in J/kgK" annotation (Dialog(tab="PEN"));
  parameter SI.SpectralEmissivity epsilon_pen=0.8 "emissivity of Anode-Electrolyte-Cathode unit" annotation (Dialog(tab = "PEN"));

    // FC parameters
  replaceable package FCMedium = Medium.Fuel_CH4          annotation (Dialog(tab="Fuel Channel"));
  parameter Real porFC=0.867 "Porosity in FC" annotation (Dialog(tab="Fuel Channel"));
  parameter Real Nu_PENfc=12 "Nusselt number fuel channel on PEN side" annotation (Dialog(tab="Fuel Channel"));
  parameter Real Nu_ICfc=9.86 "Nusselt Number fuel channel on IC side" annotation (Dialog(tab="Fuel Channel"));
  parameter SI.ThermalConductivity kFoam=3.576 "Thermal Conductivity of Ni Foam" annotation (Dialog(tab="Fuel Channel"));
  parameter SI.SpecificHeatCapacity cpFoam=440       "Specific heat capacity of Ni in FC" annotation (Dialog(tab="Fuel Channel"));
  parameter SI.Density rhoFoam(displayUnit="kg/m3")=8908 "Density of Ni in FC" annotation (Dialog(tab="Fuel Channel"));
  parameter ThermoPower.Units.HydraulicResistance HyR_Fuel=2000/2     "Hydraulic resistance of fuel channel of each single cell in stack" annotation (Dialog(tab="Fuel Channel"));

  // AC parameters
  replaceable package ACMedium = Medium.Air_Medium         annotation (Dialog(tab="Air Channel"));
  parameter Real porAC=1 "Porosity in AC" annotation (Dialog(tab="Air Channel"));
  parameter Real pDropAC=0 "Pressure drop factor, pOut=(1-pDropAC)*pIn" annotation (Dialog(tab="Air Channel"));
  parameter Real Nu_PENac=8.235 "Nusselt number air channel on PEN side" annotation (Dialog(tab="Air Channel"));
  parameter Real Nu_ICac=7.54 "Nusselt Number air channel on IC side" annotation (Dialog(tab="Air Channel"));
  parameter SI.ThermalConductivity kAC=0.0262 "Thermal Conductivity of air channel check with Faisal" annotation (Dialog(tab="End / Intermediate plates"));
  parameter ThermoPower.Units.HydraulicResistance HyR_Air=55000 "Hydraulic resistance of air channel of each single cell in stack" annotation (Dialog(tab="Air Channel"));

  // Interconnect parameters
  parameter SI.ThermalConductivity kIC=0.2812 "Thermal Conductivity of solid parallel to windows" annotation (Dialog(tab="Interconnect"));
  parameter SI.SpecificHeatCapacity cpIC=463.8 "Specific heat capacity of solid parallel to windows" annotation (Dialog(tab="Interconnect"));
  parameter SI.Density rhoIC(displayUnit="kg/m3")=1330 "Density of solid parallel to windows" annotation (Dialog(tab="Interconnect"));
  parameter SI.SpectralEmissivity epsilon_ic=0.1 "emissivity of surface of IC" annotation (Dialog(tab="Interconnect"));

  //End / Intermediate plates parameters
  parameter SI.ThermalConductivity kTopPlate=kIC "Thermal Conductivity of solid parallel to windows" annotation (Dialog(tab="End / Intermediate plates"));
  parameter SI.SpecificHeatCapacity cpTopPlate=cpIC "Specific heat capacity of solid parallel to windows" annotation (Dialog(tab="End / Intermediate plates"));
  parameter SI.Density rhoTopPlate(displayUnit="kg/m3")=rhoIC "Density of solid parallel to windows" annotation (Dialog(tab="End / Intermediate plates"));
  parameter SI.ThermalConductivity kBottomPlate=kIC "Thermal Conductivity of solid parallel to windows" annotation (Dialog(tab="End / Intermediate plates"));
  parameter SI.SpecificHeatCapacity cpBottomPlate=cpIC "Specific heat capacity of solid parallel to windows" annotation (Dialog(tab="End / Intermediate plates"));
  parameter SI.Density rhoBottomPlate(displayUnit="kg/m3")=rhoIC "Density of solid parallel to windows" annotation (Dialog(tab="End / Intermediate plates"));
  parameter SI.ThermalConductivity kIntermediatePlate=kIC "Thermal Conductivity of solid parallel to windows" annotation (Dialog(tab="Intermediate Plates"));
  parameter SI.SpecificHeatCapacity cpIntermediatePlate=cpIC "Specific heat capacity of solid parallel to windows" annotation (Dialog(tab="Intermediate Plates"));
  parameter SI.Density rhoIntermediatePlate(displayUnit="kg/m3")=rhoIC "Density of solid parallel to windows" annotation (Dialog(tab="Intermediate Plates"));

  // Simplification
  parameter SI.Length lzCell=lZpen + lZfc + lZIC + lZac "total height of a cell" annotation (Dialog(tab="Simplification"));
  parameter SI.Density rhoSimple=(lZpen*rhoPen + lZfc*rhoFoam*porFC + lZIC*rhoIC)/(lZpen + lZfc + lZIC + lZac) "dl" annotation (Dialog(tab="Simplification"));
  parameter SI.SpecificHeatCapacity cpSimple=(lZpen*rhoPen*cpPen + lZfc*rhoFoam*porFC*cpFoam + lZIC*rhoIC*cpIC)/(rhoSimple*lzCell) "height and density weighted average specific heat capacity of simple cell" annotation (Dialog(tab="End / Intermediate plates"));

  parameter SI.ThermalConductivity LambdaAir=72.98e-3 "Thermal conductivity of Gas in Air Channel" annotation (Dialog(tab="Simplification"));
  parameter SI.ThermalConductivity LambdaACConv=2*LambdaAir/(1/Nu_ICac + 1/Nu_PENac) " Thermal Conductivity of simple cell transversal" annotation (Dialog(tab="Simplification"));

  parameter SI.Temperature TIcRad = 815+273.15 "IC temperature assumption for radiative heat transfer" annotation (Dialog(tab="Simplification"));
  parameter SI.Temperature TPenRad=845 + 273.15 "PEN temperature assumption for radiative heat transfer" annotation (Dialog(tab="Simplification"));
  parameter SI.ThermalConductivity LambdaACRad = Modelica.Constants.sigma *epsilon_pen*epsilon_ic*(TPenRad^4-TIcRad^4)/(epsilon_pen+epsilon_ic-epsilon_ic*epsilon_pen)/(TPenRad-TIcRad)*lZac " Thermal Conductivity of simple cell transversal" annotation (Dialog(tab="Simplification"));
  parameter SI.ThermalConductivity LambdaAC = ((lY - lYac)*kIC+lYac*LambdaACRad+lYac*LambdaACConv)/lY "total thermal conductivity of air channel"  annotation (Dialog(tab="Simplification"));

  parameter SI.ThermalConductivity LambdaFuel=0.0935 "Thermal conductivity of Gas in Air Channel" annotation (Dialog(tab="Simplification"));
  parameter SI.ThermalConductivity LambdaFCConv=2*LambdaFuel/(1/Nu_ICfc + 1/Nu_PENfc) " Thermal Conductivity of simple cell transversal" annotation (Dialog(tab="Simplification"));

  parameter SI.ThermalConductivity LambdaFC = kFoam+LambdaFCConv "total thermal conductivity of air channel"  annotation (Dialog(tab="Simplification"));

  parameter Integer nCellVertMult[nSimpMult,N]=transpose(fill(verticalBlockSize, N)) "number of blocks of simplified cells";
  parameter SI.ThermalConductivity kSimpleLong=(lZpen*lY*kPen + lZfc*lY*LambdaFC + (lZIC*lY + lZac*(lY - lYac))*kIC + lZac*lYac*LambdaAC)/(lY*(lZpen + lZfc + lZIC + lZac)) " Thermal Conductivity of simple cell transversal" annotation (Dialog(tab="Simplification"));
  parameter SI.ThermalConductivity kSimpleTrans=(lZpen + lZfc + lZIC + lZac)/(lZpen/kPen + lZfc/LambdaFC + lZIC/kIC + lZac*lY/(LambdaAC*lYac + kIC*(lY - lYac))) "          0.2812 Thermal Conductivity of simple cell longitudinal" annotation (Dialog(tab="Simplification"));

  //parameter SI.ThermalConductivity kSimpleTrans=(lZpen*lY*kPen + lZfc*lY*kFoam + (lZIC*lY + lZac*(lY - lYac))*kIC + lZac*lYac*kAC)/(lY*(lZpen + lZfc + lZIC + lZac)) " Thermal Conductivity of simple cell transversal" annotation (Dialog(tab="Simplification"));
  //parameter SI.ThermalConductivity kSimpleLong=(lZpen + lZfc + lZIC + lZac)/(lZpen/kPen + lZfc/kFoam + lZIC/kIC + lZac/(kAC*lYac + kIC*(lY - lYac))) "          0.2812 Thermal Conductivity of simple cell longitudinal" annotation (Dialog(tab="Simplification"));
  parameter Integer nSimpMult=size(verticalBlockSize,1) "number of blocks of simplified cells";
  parameter Integer nIMPs = integer((floor((Ncell - 1)/10))) "Number of intermediate plates";
protected
  parameter Integer stackingHelper[Ncell] = {integer(i/(Ncell)*(Ncell-NdetailedCell) + 0.5) - integer((i - 1)/(Ncell)*(Ncell-NdetailedCell) + 0.5) for i in 1:Ncell} "define which cell is detailed (0) or simplified (1)";
  parameter Integer j_simp[Ncell] = {sum(stackingHelper[1:i]) for i in 1:Ncell} "counter for simplified cells";
  parameter Integer j_detl[Ncell] = {i - j_simp[i] for i in 1:Ncell} "counter for detailed cells";
  parameter Boolean isSimpl[Ncell] = {stackingHelper[i] > 0.5 for i in 1:Ncell} "turn integer in stackingHelper to boolean";
  parameter Integer distancesCellsToCells[Ncell,Ncell]=transpose(fill(1:Ncell,Ncell))-fill(1:Ncell,Ncell) "distances from each cell to each other cell";
  parameter Integer j_simp_positions[Ncell] = {Modelica.Math.Vectors.find(i,j_simp) for i in 1:Ncell} "list of simplified cells and trailing zeroes";
  parameter Integer j_detl_positions[Ncell] = {Modelica.Math.Vectors.find(i,j_detl) for i in 1:Ncell} "list of detailed cells and trailing zeroes";
  parameter Integer verticalBlockRanges[nSimpMult,2] = {{j_simp_positions[sum(verticalBlockSize[1:i])]-verticalBlockSize[i]+1,j_simp_positions[sum(verticalBlockSize[1:i])]} for i in 1:size(verticalBlockSize,1)} "list with cells indices of limiting upper and lower cell in each simplified vertical block";
  parameter Integer j_next_detailed_cell[nSimpMult-1] = {if j_detl[verticalBlockRanges[i+1,2]]==0 then j_detl_positions[1] else j_detl_positions[max(1,j_detl[verticalBlockRanges[i+1,2]])] for i in 1:nSimpMult-1} "list the next detailed cell corresponding to each cell. Avoid translation error by using max(1,index), where index can be 0.";
  parameter Integer verticalBlockRangesHighActual[nSimpMult] = cat(1,{if j_detl[verticalBlockRanges[i+1,2]]==1+j_detl[verticalBlockRanges[i,2]] then j_next_detailed_cell[i]-1 else verticalBlockRanges[i,2] for i in 1:nSimpMult-1},{verticalBlockRanges[nSimpMult,2]}) "corrected limiting upper cell in each simplified vertical block. If verticalBlockSize[:] is such that every detailed cell is placed exactly between two simplified cell vertical blocks, this is equal to verticalBlockRanges[:]. Otherwise it is adapted correspondingly.";
  parameter Integer verticalBlockRangesLowActual[nSimpMult] = cat(1,{verticalBlockRanges[1,1]},{if j_detl[verticalBlockRangesHighActual[i]]==1+j_detl[verticalBlockRangesHighActual[i-1]] then verticalBlockRangesHighActual[i-1]+2 else verticalBlockRangesHighActual[i-1]+(j_detl[verticalBlockRangesHighActual[i]]-j_detl[verticalBlockRangesHighActual[i-1]])+1 for i in 2:nSimpMult}) "corrected limiting lower cell in each simplified vertical block. If verticalBlockSize[:] is such that every detailed cell is placed exactly between two simplified cell vertical blocks, this is equal to verticalBlockRanges[:,1]. Otherwise it is adapted correspondingly.";
  parameter Integer verticalBlockRangesActual[nSimpMult,2] = {{verticalBlockRangesLowActual[i],verticalBlockRangesHighActual[i]} for i in 1:nSimpMult} "corrected limiting lower cell in each simplified vertical block. If verticalBlockSize[:] is such that every detailed cell is placed exactly between two simplified cell vertical blocks, this is equal to verticalBlockRanges. Otherwise it is adapted correspondingly.";
  parameter Integer distancesCellsToVerticalBlocks[nSimpMult,Ncell] = {{abs(sum(distancesCellsToCells[i,l_ind] for i in verticalBlockRangesActual[j,1]:verticalBlockRangesActual[j,2])) for l_ind in 1:Ncell} for j in 1:nSimpMult} "helper with distances of each cell (columns) to each simplified cells vertical block (rows)";
  parameter Integer distancesHelperVerticalBlocks[nSimpMult,Ncell] = {{if isSimpl[i] then max(distancesCellsToVerticalBlocks) else distancesCellsToVerticalBlocks[j,i] for i in 1:Ncell} for j in 1:nSimpMult} "helper with total distance of each simplified cell block to each detailed cell, as the sum of the distances of each simplified cell contained in the block";
  parameter Integer j_corresponding[nSimpMult] = {Modelica.Math.Vectors.find(min(distancesHelperVerticalBlocks[i,:]),distancesHelperVerticalBlocks[i,:]) for i in 1:nSimpMult} "list detailed cells corresponding to each simplified cell vertical block";
  parameter Integer j_verticalBlock[Ncell] = {Modelica.Math.Vectors.find(1,{if verticalBlockRangesActual[j,1]-i==0 or verticalBlockRangesActual[j,2]-i==0 or (verticalBlockRangesActual[j,1]-i<0 and verticalBlockRangesActual[j,2]-i>0) then 1 else 0 for j in 1:size(verticalBlockRangesActual,1)}) for i in 1:Ncell} "index of the simplified block corresponding to each cell (0 for none)";
  parameter Integer j_SimpDet_z0[Ncell] = cat(1,{0},{sum({if -(stackingHelper[i+1]-stackingHelper[i])>0 then 1 else 0 for i in 1:j}) for j in 1:Ncell-1}) "indexes of interfaces between bottom of detailed cell and top of previous simplified cell";
  parameter Integer j_SimpDet_z1[Ncell] = cat(1,{sum({if stackingHelper[i]-stackingHelper[i-1]>0 then 1 else 0 for i in 2:j}) for j in 2:Ncell},{sum({if stackingHelper[i]-stackingHelper[i-1]>0 then 1 else 0 for i in 2:Ncell})}) "indexes of interfaces between top of detailed cell and bottom of next simplified cell";
  parameter Integer Ncontacts_SimpDet_z0 = max(j_SimpDet_z0) "total contact points between simplified and detailed cells at the detailed-cell vertical bottom position";
  parameter Integer Ncontacts_SimpDet_z1 = max(j_SimpDet_z1) "total contact points between simplified and detailed cells at the detailed-cell vertical top position";
  parameter Integer simpleCellHeatPortIndex[Ncell] = {j-verticalBlockRangesActual[max(1,j_verticalBlock[j]),1]+1 for j in 1:Ncell};

public
  DummyManifold Fuel_Manifold(
    redeclare package Medium = Medium.Fuel_CH4,
    final nPorts_b=Ncell,
    nDummyPorts_b=Ncell - NdetailedCell)
    annotation (Placement(transformation(extent={{-74,50},{-62,82}})));
  DummyManifold Air_Manifold(
    redeclare package Medium = Medium.Air_Medium,
    final nPorts_b=Ncell,
    nDummyPorts_b=Ncell - NdetailedCell)
    annotation (Placement(transformation(extent={{-76,-80},{-64,-48}})));
  ThermoPower.Gas.FlangeB fuelOut(redeclare package Medium = Medium.Fuel_CH4)
    annotation (Placement(transformation(extent={{110,56},{130,76}}),
        iconTransformation(extent={{110,56},{130,76}})));
  ThermoPower.Gas.FlangeB airOut(redeclare package Medium = Medium.Air_Medium)
    annotation (Placement(transformation(extent={{110,-72},{130,-52}}),
        iconTransformation(extent={{110,-72},{130,-52}})));
  Flow.Manifold_out manifold_outFuel(redeclare package Medium = Medium.Fuel_CH4,
      nPorts_a=NdetailedCell)
    annotation (Placement(transformation(extent={{36,10},{56,30}})));
  Flow.Manifold_out manifold_outAir(redeclare package Medium =
        Medium.Air_Medium, nPorts_a=NdetailedCell)
    annotation (Placement(transformation(extent={{36,-18},{56,2}})));
  ThermoPower.Gas.FlangeA fuelIn(redeclare package Medium = Medium.Fuel_CH4)
    annotation (Placement(transformation(extent={{-138,56},{-118,76}}),
        iconTransformation(extent={{-138,56},{-118,76}})));
  ThermoPower.Gas.FlangeA airIn(redeclare package Medium = Medium.Air_Medium)
    annotation (Placement(transformation(extent={{-140,-80},{-120,-60}}),
        iconTransformation(extent={{-140,-80},{-120,-60}})));

  Modelica.Electrical.Analog.Interfaces.PositivePin pin_p
    annotation (Placement(transformation(extent={{-136,-60},{-116,-40}}),
        iconTransformation(extent={{-136,-60},{-116,-40}})));
  Modelica.Electrical.Analog.Interfaces.NegativePin pin_n
    annotation (Placement(transformation(extent={{-136,-18},{-116,2}}),
        iconTransformation(extent={{-136,-18},{-116,2}})));

  ThermoPower.Thermal.DHTVolumes dHT_top(N=N) annotation (Placement(
        transformation(extent={{-22,78},{-2,98}}),iconTransformation(extent={{-50,
            30},{44,48}})));
  ThermoPower.Thermal.DHTVolumes dHT_bottom(N=N) annotation (Placement(
        transformation(extent={{-32,-96},{-12,-76}}),
                                                    iconTransformation(extent={{
            -50,-52},{44,-34}})));

  ThermoPower.Thermal.DHTVolumes dHT_y0(N=N) annotation (Placement(
        transformation(extent={{6,100},{64,110}}),   iconTransformation(extent={{-56,100},{64,110}})));

  ThermoPower.Thermal.DHTVolumes dHT_y1(N=N) annotation (Placement(
        transformation(extent={{12,-110},{70,-100}}),  iconTransformation(
          extent={{-60,-110},{60,-100}})));
  ThermoPower.Thermal.DHTVolumes dHT_x0(N=Ncell)
    annotation (Placement(transformation(extent={{-132,-38},{-112,-18}}),
        iconTransformation(extent={{-132,-36},{-120,4}})));
  ThermoPower.Thermal.DHTVolumes dHT_xN(N=Ncell)
    annotation (Placement(transformation(extent={{108,-10},{128,10}}),
        iconTransformation(extent={{116,-20},{128,20}})));
  FlowInterpolator flowInterpolator(
    nCell=Ncell,
    nSimplified=Ncell - NdetailedCell,
    nNonUnitOrSimpCells=0,
    isRedu=isSimpl,
    TStart=TStart,
    pStart=pStart,
    XFuelStart=xStartFC,
    redeclare package Fuel = FCMedium,
    redeclare package Air = ACMedium)
    annotation (Placement(transformation(extent={{60,-4},{80,16}})));
  OpenTEMPEST.Flow.SensGasProperty sensCpFuel(
    mfOutput=false,
    pOutput=false,
    hOutput=false,
    XOutput=false,
    YOutput=false,
    HfOutput=false,
    redeclare package Medium = Medium.Fuel_CH4,
    cpOutput=true,
    lambdaOutput=false,
    etaOutput=false,
    TOutput=false,
    rhoOutput=false)
    annotation (Placement(transformation(extent={{-108,60},{-88,80}})));
  OpenTEMPEST.Flow.SensGasProperty sensCpAir(
    mfOutput=false,
    pOutput=false,
    hOutput=false,
    XOutput=false,
    YOutput=false,
    HfOutput=false,
    redeclare package Medium = Medium.Air_Medium,
    cpOutput=true,
    lambdaOutput=false,
    etaOutput=false,
    TOutput=false,
    rhoOutput=false)
    annotation (Placement(transformation(extent={{-108,-70},{-88,-50}})));
    SI.EnergyFlowRate EfSum "check energy balance through sum of energy flows";
    SI.MassFlowRate massSum "check mass balance through sum of mass flows";
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature TdetailedToFVI[
    NdetailedCell]
    annotation (Placement(transformation(extent={{70,-40},{60,-30}})));
  Modelica.Blocks.Sources.RealExpression realTdetailedOut[NdetailedCell](
    y=(detailedCells.fuelChannel.hTNi_xN.T))
    annotation (Placement(transformation(extent={{144,-38},{108,-20}})));

  Cell.Cell1D.SimplifiedCellBlock simplifiedCellBlock[nSimpMult](
    each N=N,
    each TStart=TStart_IC,
    nCellVertMult=verticalBlockSize,
    each lX=lX,
    each lY=lY,
    each lZ=lzCell,
    each cp_solidCustom=cpIC,
    each rho_solidCustom=rhoIC,
    each rho_custom=rhoSimple,
    each cp_custom=cpSimple,
    each k_solidCustomLong=kSimpleLong,
    each k_solidCustomTrans=kSimpleTrans)
    annotation (Placement(transformation(extent={{-12,14},{8,34}})));
    //each matOpt=TEMPEST.ECReactorModels.Cell.baseModelsFV_SOC.solidMaterialOptions.Crofer22APU,

  Cell.Cell1D.Cell detailedCells[NdetailedCell](
    redeclare model Electrochem = Electrochem,
    each N=N,
    each TStart=TStart,
    each Jstart=JStart,
    each pStartInFC=pStartInFC,
    each pStartOutFC=pStartOutFC,
    each xStartFC=xStartFC,
    each pStartAC=pStartAC,
    each xStartAC=xStartAC,
    each TStart_pen=TStart_pen,
    each TStart_FCin=TStart_FCin,
    each TStart_FCout=TStart_FCout,
    each TStart_ACin=TStart_ACin,
    each TStart_ACout=TStart_ACout,
    each TStart_IC=TStart_IC,
    each epsilon_pen=epsilon_pen,
    each lYfc=lYfc,
    each lZfc=lZfc,
    each lX=lX,
    each lY=lY,
    each Nu_PENfc=Nu_PENfc,
    each Nu_ICfc=Nu_ICfc,
    each lXac=lXac,
    each lYac=lYac,
    each lZac=lZac,
    each lYpen=lYpen,
    each lZpen=lZpen,
    each Nu_PENac=Nu_PENac,
    each Nu_ICac=Nu_ICac,
    each lZIC=lZIC/2,
    each kPen=kPen,
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
    each LambdaIC=kIC,
    each cpIC=cpIC,
    each rhoIC=rhoIC,
    each epsilon_ic=epsilon_ic,
    each HyR_Fuel=HyR_Fuel,
    each HyR_Air=HyR_Air,
    each dpNomFC=dpNomFC,
    redeclare function fluxInterp = fluxInterp)
    annotation (Placement(transformation(extent={{-22,-34},{16,2}})));

  OpenTEMPEST.Heat.Solid1D plates[nIMPs + 2](
    redeclare package SolidMat = TEMPEST.Solid.Material.Custom,
    each N=N,
    each TstartX0=TStart,
    each TstartXN=TStart,
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
    annotation (Placement(transformation(extent={{48,36},{68,56}})));

  ThermoPower.Thermal.HT hTPlates_x1[nIMPs+2]
    annotation (Placement(transformation(extent={{110,30},{130,50}}), iconTransformation(extent={{110,34},{130,54}})));
  ThermoPower.Thermal.HT hTPlates_x0[nIMPs+2]
    annotation (Placement(transformation(extent={{-126,34},{-106,54}}), iconTransformation(extent={{-134,34},{-114,54}})));
//protected
  Heat.HTs_DHT hTs_DHT_x0(each N=Ncell)
    annotation (Placement(transformation(extent={{-68,-22},{-80,-10}})));
  Heat.HTs_DHT hTs_DHT_xN(each N=Ncell)
    annotation (Placement(transformation(extent={{72,-22},{84,-10}})));
  Modelica.Electrical.Analog.Sources.SignalVoltage signalVoltage annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={-96,-4})));
  Modelica.Blocks.Sources.RealExpression USimp(y=-sum(simplifiedCellBlock[:].UCalc[:] .* nCellVertMult[:, :])/N)
annotation (Placement(transformation(extent={{-172,12},{-80,34}})));

equation
   assert(not Modelica.Math.BooleanVectors.anyTrue(Modelica.Math.BooleanVectors.anyTrue({{j_detl_positions[i]-verticalBlockRangesActual[j,1]>=0 and verticalBlockRangesActual[j,2]-j_detl_positions[i]>=0 for i in 1:size(j_detl_positions,1)} for j in 1:size(verticalBlockRangesActual,1)})),"error: some vertical blocks are overlapping with detailed cell positions.");
   assert(sum(verticalBlockSize)==Ncell-NdetailedCell,"Vertical block sizes must total up to Ncell-NdetailedCell");
   assert(Modelica.Math.BooleanVectors.anyTrue(Modelica.Math.BooleanVectors.anyTrue({{verticalBlockRangesLowActual[i]==j*10 or verticalBlockRangesHighActual[i]==j*10 for i in 1:size(verticalBlockRangesLowActual,1)} for j in 1:nIMPs})) or
                                                                                                                                                                                                        Ncell <=10,"error: some vertical blocks are overlapping with intermediate plate positions");

      EfSum = fuelIn.m_flow*inStream(fuelIn.h_outflow) + fuelOut.m_flow*fuelOut.h_outflow
                + airIn.m_flow*inStream(airIn.h_outflow) + airOut.m_flow*airOut.h_outflow
                + (pin_p.v-pin_n.v)*pin_n.i
                + sum(dHT_x0.Q[:]) + sum(dHT_xN.Q[:])
                + sum(dHT_top.Q[:]) + sum(dHT_bottom.Q[:])
                + sum(dHT_y0.Q[:]) + sum(dHT_y1.Q[:])
                + sum(hTPlates_x1[:].Q_flow)
                + sum(hTPlates_x0[:].Q_flow); // useful for steady state, but not subtracting all heat capacitors, as in ECReactorModels.Stack.SimplifiedStack.SimplifiedStackVertMult

      massSum = fuelIn.m_flow + fuelOut.m_flow + airIn.m_flow + airOut.m_flow;

  // simplification voltage and temperature
  flowInterpolator.TOutAirSimp = sum(sum(simplifiedCellBlock[j].simpFlowInfoPortOut[N].TAir.*nCellVertMult[j,N] for j in 1:nSimpMult))/Ncell+sum(detailedCells.airChannel.Gas[N].T)/Ncell;
  flowInterpolator.TOutFuelSimp = sum(sum(simplifiedCellBlock[j].simpFlowInfoPortOut[N].TFuel.*nCellVertMult[j,N] for j in 1:nSimpMult))/Ncell+sum(detailedCells.fuelChannel.hTNi_xN.T)/Ncell;

  // pass information to activeBlocks
  for j in 1:nSimpMult loop
    for m_idx in 1:N loop
     // here: j_detl[j_corresponding[j]] is the number of detailed cell corresponding to detailed cell index, corresponding to current simplified cells block
     simplifiedCellBlock[j].reactionHeatInfoPort[m_idx].H_r = 0;  // not used //detailedCells[k_idx].activeAreaLayer.sPCellActiveBlockRS[i,j].pen.q_electroChem[l_idx];
     simplifiedCellBlock[j].reactionHeatInfoPort[m_idx].LHVFlow = 0; // not used
     simplifiedCellBlock[j].reactionHeatInfoPort[m_idx].Vop = detailedCells[j_detl[j_corresponding[j]]].pen.electrochem[m_idx].Uop;
     simplifiedCellBlock[j].reactionHeatInfoPort[m_idx].Vid = detailedCells[j_detl[j_corresponding[j]]].pen.electrochem[m_idx].Uideal;
     simplifiedCellBlock[j].reactionHeatInfoPort[m_idx].I = detailedCells[j_detl[j_corresponding[j]]].pen.Iv[m_idx];
     simplifiedCellBlock[j].reactionHeatInfoPort[m_idx].TPEN = detailedCells[j_detl[j_corresponding[j]]].pen.electrochem[m_idx].Tpen;

     simplifiedCellBlock[j].simpFlowInfoPortIn[m_idx].mfAir= +detailedCells[j_detl[j_corresponding[j]]].airChannel.mfv[m_idx];
     simplifiedCellBlock[j].simpFlowInfoPortIn[m_idx].mfFuel= +detailedCells[j_detl[j_corresponding[j]]].fuelChannel.mfv[m_idx];
     simplifiedCellBlock[j].simpFlowInfoPortIn[m_idx].HfAirRef = +detailedCells[j_detl[j_corresponding[j]]].airChannel.mfv[m_idx] * detailedCells[j_detl[j_corresponding[j]]].airChannel.hv[m_idx];
     simplifiedCellBlock[j].simpFlowInfoPortIn[m_idx].HfFuelRef =+detailedCells[j_detl[j_corresponding[j]]].fuelChannel.mfv[m_idx] * detailedCells[j_detl[j_corresponding[j]]].fuelChannel.hv[m_idx];

     if m_idx==1 then
       simplifiedCellBlock[j].simpFlowInfoPortIn[m_idx].cpAir= sensCpAir.cp;
       simplifiedCellBlock[j].simpFlowInfoPortIn[m_idx].cpFuel= sensCpFuel.cp;
       simplifiedCellBlock[j].simpFlowInfoPortIn[m_idx].TAirRef = detailedCells[j_detl[j_corresponding[j]]].airChannel.Gas[1].T;
       simplifiedCellBlock[j].simpFlowInfoPortIn[m_idx].TFuelRef = detailedCells[j_detl[j_corresponding[j]]].fuelChannel.hTNi_x0.T;
       simplifiedCellBlock[j].simpFlowInfoPortIn[m_idx].TAir = detailedCells[j_detl[j_corresponding[j]]].airChannel.Gas[1].T;
       simplifiedCellBlock[j].simpFlowInfoPortIn[m_idx].TFuel = detailedCells[j_detl[j_corresponding[j]]].fuelChannel.hTNi_x0.T;
     elseif m_idx>1 then
       simplifiedCellBlock[j].simpFlowInfoPortIn[m_idx].cpAir= ACMedium.specificHeatCapacityCp(detailedCells[j_detl[j_corresponding[j]]].airChannel.Gas[m_idx-1].state);
       simplifiedCellBlock[j].simpFlowInfoPortIn[m_idx].cpFuel= FCMedium.specificHeatCapacityCp(detailedCells[j_detl[j_corresponding[j]]].fuelChannel.Gas[m_idx-1].state);
       simplifiedCellBlock[j].simpFlowInfoPortIn[m_idx].TAirRef = detailedCells[j_detl[j_corresponding[j]]].airChannel.Gas[m_idx-1].T;
       simplifiedCellBlock[j].simpFlowInfoPortIn[m_idx].TFuelRef = detailedCells[j_detl[j_corresponding[j]]].fuelChannel.Gas[m_idx-1].T;
       simplifiedCellBlock[j].simpFlowInfoPortIn[m_idx].TAir = simplifiedCellBlock[j].simpFlowInfoPortOut[m_idx-1].TAir;
       simplifiedCellBlock[j].simpFlowInfoPortIn[m_idx].TFuel = simplifiedCellBlock[j].simpFlowInfoPortOut[m_idx-1].TFuel;
     end if;

     simplifiedCellBlock[j].simpFlowInfoPortOut[m_idx].cpAir= ACMedium.specificHeatCapacityCp(detailedCells[j_detl[j_corresponding[j]]].airChannel.Gas[m_idx].state);
     simplifiedCellBlock[j].simpFlowInfoPortOut[m_idx].cpFuel= FCMedium.specificHeatCapacityCp(detailedCells[j_detl[j_corresponding[j]]].fuelChannel.Gas[m_idx].state);

     simplifiedCellBlock[j].simpFlowInfoPortOut[m_idx].mfAir= -detailedCells[j_detl[j_corresponding[j]]].airChannel.mfv[m_idx+1];
     simplifiedCellBlock[j].simpFlowInfoPortOut[m_idx].mfFuel= -detailedCells[j_detl[j_corresponding[j]]].fuelChannel.mfv[m_idx+1];

     simplifiedCellBlock[j].simpFlowInfoPortOut[m_idx].HfAirRef = -detailedCells[j_detl[j_corresponding[j]]].airChannel.mfv[m_idx+1] * detailedCells[j_detl[j_corresponding[j]]].airChannel.hv[m_idx+1];
     simplifiedCellBlock[j].simpFlowInfoPortOut[m_idx].HfFuelRef = -detailedCells[j_detl[j_corresponding[j]]].fuelChannel.mfv[m_idx+1] * detailedCells[j_detl[j_corresponding[j]]].fuelChannel.hv[m_idx+1];

     if m_idx<N then
       simplifiedCellBlock[j].simpFlowInfoPortOut[m_idx].TAirRef = detailedCells[j_detl[j_corresponding[j]]].airChannel.Gas[m_idx].T;
       simplifiedCellBlock[j].simpFlowInfoPortOut[m_idx].TFuelRef = detailedCells[j_detl[j_corresponding[j]]].fuelChannel.Gas[m_idx].T;
       simplifiedCellBlock[j].simpFlowInfoPortOut[m_idx].TAir = simplifiedCellBlock[j].solid1D.T[m_idx];
       simplifiedCellBlock[j].simpFlowInfoPortOut[m_idx].TFuel = simplifiedCellBlock[j].solid1D.T[m_idx];
     elseif m_idx==N then
       simplifiedCellBlock[j].simpFlowInfoPortOut[m_idx].TAirRef = detailedCells[j_detl[j_corresponding[j]]].airChannel.Gas[N].T;
       simplifiedCellBlock[j].simpFlowInfoPortOut[m_idx].TFuelRef = detailedCells[j_detl[j_corresponding[j]]].fuelChannel.hTNi_xN.T;
       simplifiedCellBlock[j].simpFlowInfoPortOut[m_idx].TAir = simplifiedCellBlock[j].solid1D.hT_xN.T;
       simplifiedCellBlock[j].simpFlowInfoPortOut[m_idx].TFuel = simplifiedCellBlock[j].solid1D.hT_xN.T;
     end if;
    end for;
   end for;

  // Electrical connections - External
  connect(USimp.y, signalVoltage.v) annotation (Line(points={{-75.4,23},{-61.7,23},{-61.7,-4},{-84,-4}}, color={0,0,127}));
  connect(signalVoltage.n, detailedCells[1].pin_n) annotation (Line(points={{-96,-14},{-96,-19.96},{-21.24,-19.96}}, color={0,0,255}));
  connect(signalVoltage.p, pin_n) annotation (Line(points={{-96,6},{-106,6},{-106,-8},{-126,-8}}, color={0,0,255}));
  connect(detailedCells[NdetailedCell].pin_p, pin_p) annotation (Line(points={{-21.24,-23.92},{-52,-23.92},{-52,-44},{-118,-44},{-118,-50},{-126,-50}},
                                                                                                                                        color={0,0,255}));

  // Electrical connections - Inter cell (detailed cells)

  for j in 1:Ncell loop
    if not isSimpl[j] then
      // here: j_detl[j] is the detailed cell index corresponding to the complete cells index.
      if NdetailedCell > 1 and j_detl[j] < NdetailedCell then // more than one detailed cell in between
      // Electrical Connections (Series/Internal)
        connect(detailedCells[j_detl[j]].pin_p, detailedCells[j_detl[j] + 1].pin_n);
      end if;
     end if;
  end for;

  // connect z-direction HT ports
  connect(plates[1].dhT_z0, dHT_bottom) annotation (Line(points={{51,39},{51,-72},{-22,-72},{-22,-86}}, color={255,127,0}));
  connect(plates[nIMPs + 2].dhT_z1, dHT_top) annotation (Line(points={{51,53},{51,78.5},{-12,78.5},{-12,88}}, color={255,127,0}));

  //connect endplates
  //topplate
  if isSimpl[Ncell] then
    // connect simplified cell to top plate
    connect(plates[nIMPs + 2].dhT_z0, simplifiedCellBlock[j_verticalBlock[Ncell]].dHTactive_z1);
  else
    // connect detailed cell to top plate
    connect(plates[nIMPs + 2].dhT_z0, detailedCells[j_detl[Ncell]].dHT_z1) annotation (Line(points={{51,39},
            {51,36},{-18,36},{-18,6},{-12.88,6},{-12.88,-2.68}},                                                                                                 color={255,127,0}));
  end if;
  //bottomplate
  if isSimpl[1] then
    // connect simplified cell to bottom plate below
    connect(plates[1].dhT_z1, simplifiedCellBlock[j_verticalBlock[1]].dHTactive_z0);
  else
    // connect detailed cell to bottom plate below
    connect(plates[1].dhT_z1, detailedCells[1].dHT_z0) annotation (Line(points={{51,53},
            {51,60},{30,60},{30,-46},{-10,-46},{-10,-28.6},{-12.5,-28.6}},                                         color={255,127,0}));
  end if;

  // consider possible cases: detl[j] to simp[j+1], simp[j] to detl[j+1], simp[1] to bottom plate, etc.
  for j in 1:Ncell-1 loop
    if isSimpl[j] then
      if mod(j, 10) == 0 and j > 10 then
        //connect simplified block to intermediate plate above
        connect(simplifiedCellBlock[j_verticalBlock[j]].dHTactive_z1, plates[integer(1+j/10)].dhT_z0);
        if isSimpl[j + 1] then
          // connect intermediate plate to simplified block above
          connect(plates[integer(1+j/10)].dhT_z1, simplifiedCellBlock[j_verticalBlock[j + 1]].dHTactive_z0);
        else
          // connect intermediate plate to detailed cell above
          connect(plates[integer(1+j/10)].dhT_z1, detailedCells[j_detl[j] + 1].dHT_z0);
        end if;
      elseif not isSimpl[j + 1] and verticalBlockRangesActual[max(1, j_verticalBlock[j]), 2] - j == 0 then
        // only connect upwards if simplified cell j is at the top of its vertical block (simplified cells in between are not separate objects)
        // connect simplified cell to detailed cell above (through interface hts_DHT)
        connect(simplifiedCellBlock[j_verticalBlock[j]].dHTactive_z1, detailedCells[j_detl[j + 1]].dHT_z0);
      elseif verticalBlockRangesActual[max(1, j_verticalBlock[j]), 2] - j == 0 then
        // only connect upwards if simplified cell j is at the top of its vertical block (simplified cells in between are not separate objects)
        // connect simplified cell to simplified cell above
        connect(simplifiedCellBlock[j_verticalBlock[j]].dHTactive_z1, simplifiedCellBlock[j_verticalBlock[j + 1]].dHTactive_z0);
      end if;
    elseif not isSimpl[j] then
      if mod(j, 10) == 0 and j > 10 then
         //connect detailed cell to intermediate plate above
           connect(detailedCells[j_detl[j]].dHT_z1, plates[integer(1+j/10)].dhT_z0);
        if isSimpl[j + 1] then
          // connect intermediate plate to simplified block above
          connect(plates[integer(1+j/10)].dhT_z1, simplifiedCellBlock[j_verticalBlock[j + 1]].dHTactive_z0);
        else
          // connect intermediate plate to detailed cell above
          connect(plates[integer(1+j/10)].dhT_z1, detailedCells[j_detl[j] + 1].dHT_z0);
        end if;
      elseif isSimpl[j + 1] then
        // connect detailed cell to simplified block above
        connect(detailedCells[j_detl[j]].dHT_z1, simplifiedCellBlock[j_verticalBlock[j + 1]].dHTactive_z0);
      else
        // connect detailed cell to detailed cell above
        connect(detailedCells[j_detl[j]].dHT_z1, detailedCells[j_detl[j] + 1].dHT_z0);
      end if;
    end if;
  end for;

  // Gas flow connections
  connect(Fuel_Manifold.ports_b[1:NdetailedCell], detailedCells.fuelFlangeIn) annotation (
      Line(points={{-62,66},{-50,66},{-50,-6},{-22,-6},{-22,-6.28}},
                                                      color={0,127,255}));
  connect(Air_Manifold.ports_b[1:NdetailedCell], detailedCells.airFlangeIn) annotation (Line(
        points={{-64,-64},{-48,-64},{-48,-32},{-22,-32},{-22,-27.88}},
                                               color={0,127,255}));
  connect(manifold_outAir.port_b, flowInterpolator.airInlet)
    annotation (Line(points={{51,-8},{56,-8},{56,2},{60,2}}, color={0,127,255}));
  connect(detailedCells.airFlangeOut, manifold_outAir.ports_a[1:NdetailedCell]) annotation (Line(
        points={{16,-24.64},{34,-24.64},{34,-8},{43.4,-8}},         color={159,159,
          223}));
  connect(manifold_outFuel.port_b, flowInterpolator.fuelInlet)
    annotation (Line(points={{51,20},{56,20},{56,10},{60,10}}, color={0,127,255}));
  connect(flowInterpolator.fuelExit, fuelOut) annotation (Line(points={{78,10},{78,66},{120,66}},
                                 color={159,159,223}));
  connect(flowInterpolator.airExit, airOut) annotation (Line(points={{78,2},{78,-62},{120,-62}}, color={159,159,223}));

  connect(fuelIn, sensCpFuel.inlet)
    annotation (Line(points={{-128,66},{-104,66}}, color={159,159,223}));
  connect(detailedCells.fuelFlangeOut, manifold_outFuel.ports_a[1:NdetailedCell]) annotation (Line(
        points={{16,-6.28},{28.66,-6.28},{28.66,20},{43.4,20}},            color=
         {159,159,223}));
  connect(sensCpFuel.outlet, Fuel_Manifold.port_a)
    annotation (Line(points={{-92,66},{-74,66}}, color={159,159,223}));
  connect(airIn, sensCpAir.inlet) annotation (Line(points={{-130,-70},{-114,-70},{-114,-64},{-104,-64}},
                                  color={159,159,223}));
  connect(sensCpAir.outlet, Air_Manifold.port_a)
    annotation (Line(points={{-92,-64},{-76,-64}}, color={159,159,223}));

  // connect x-direction HT ports (flow direction)
  // Cells
  connect(hTs_DHT_xN.DHT_port, dHT_xN) annotation (Line(points={{84.6,-16},{106,-16},{106,0},{118,0}}, color={255,127,0}));
  connect(hTs_DHT_x0.DHT_port, dHT_x0) annotation (Line(points={{-80.6,-16},{-84,-16},{-84,-28},{-122,-28}}, color={255,127,0}));

  for j in 1:Ncell loop
    if not isSimpl[j] then
      // connect detailed cell inlet and outlet port in flow direction
      connect(hTs_DHT_x0.HT_ports[j], detailedCells[j_detl[j]].hT_x0) annotation (Line(points={{-66.8,-15.88},{-44,-15.88},{-44,-16},{-21.24,-16}}, color={191,0,0}));
      connect(hTs_DHT_xN.HT_ports[j], detailedCells[j_detl[j]].hT_xN) annotation (Line(points={{70.8,-15.88},{15.62,-15.88},{15.62,-16}}, color={191,0,0}));
      connect(flowInterpolator.hT[j], TdetailedToFVI[j_detl[j]].port) annotation (Line(points={{61,6},{42,6},{42,-35},{60,-35}}, color={191,0,0}));
    elseif isSimpl[j] then
      // connect simplified cell port. Blocks with simplified cells get connected to all corresponding hTports (simplified cells lumped)
      connect(hTs_DHT_x0.HT_ports[j], simplifiedCellBlock[j_verticalBlock[j]].hTactive_x0[simpleCellHeatPortIndex[j]]) annotation (Line(points={{-66.8,-15.88},{-38,-15.88},{-38,24.3},{-12.3,24.3}}, color={191,0,0}));
      connect(hTs_DHT_xN.HT_ports[j], simplifiedCellBlock[j_verticalBlock[j]].hTactive_x1[simpleCellHeatPortIndex[j]]) annotation (Line(points={{70.8,-15.88},{32,-15.88},{32,24.3},{8.1,24.3}}, color={191,0,0}));
      connect(flowInterpolator.hT[j], simplifiedCellBlock[j_verticalBlock[j]].hTactive_x1[simpleCellHeatPortIndex[j]]);
    end if;
  end for;

  /*
  connect(hTs_DHT_x0.HT_ports[2], detailedCells[1].hT_x0);
  connect(hTs_DHT_xN.HT_ports[2], detailedCells[1].hT_xN);
  connect(flowInterpolator.hT[2], TdetailedToFVI[1].port);

  connect(hTs_DHT_x0.HT_ports[1], simplifiedCellBlock[1].hTactive_x0[1]);
  connect(hTs_DHT_xN.HT_ports[1], simplifiedCellBlock[1].hTactive_x1[1]);
  connect(flowInterpolator.hT[1], simplifiedCellBlock[1].hTactive_x1[1]);

  connect(hTs_DHT_x0.HT_ports[3], simplifiedCellBlock[2].hTactive_x0[1]);
  connect(hTs_DHT_xN.HT_ports[3], simplifiedCellBlock[2].hTactive_x1[1]);
  connect(flowInterpolator.hT[3], simplifiedCellBlock[2].hTactive_x1[1]);
*/

  // Top plate
  connect(plates[nIMPs + 2].dhT_y0, dHT_y0) annotation (Line(points={{65,39},{65,
          32},{12,32},{12,105},{35,105}},                                                                        color={255,127,0}));
  connect(plates[nIMPs + 2].dhT_y1, dHT_y1) annotation (Line(points={{65,53},{65,70},{41,70},{41,-105}}, color={255,127,0}));

  // Bottom plate
  connect(plates[1].dhT_y0, dHT_y0) annotation (Line(points={{65,39},{52,39},{52,
          34},{26,34},{26,66},{35,66},{35,105}},                                         color={255,127,0}));
  connect(plates[1].dhT_y1, dHT_y1) annotation (Line(points={{65,53},{65,64},{92,
          64},{92,24},{38,24},{38,-34},{41,-34},{41,-105}},                                                         color={255,127,0}));

  //intermediate plates to external
  //y-direction
  for i in 1:nIMPs loop
    connect(plates[i+1].dhT_y0, dHT_y0) annotation (Line(points={{65,39},{65,30},{86,30},{86,82},{35,82},{35,105}}, color={255,127,0}));
    connect(plates[i+1].dhT_y1, dHT_y1) annotation (Line(points={{65,53},{65,68},{41,68},{41,-105}}, color={255,127,0}));
  end for;
  //x-direction
  for i in 1:(nIMPs+2) loop
    connect(plates[i].hT_xN, hTPlates_x1[i]) annotation (Line(points={{69,46},{98,46},{98,40},{120,40}}, color={191,0,0}));
    connect(plates[i].hT_x0, hTPlates_x0[i]) annotation (Line(points={{47,46},{-34,46},{-34,44},{-116,44}}, color={191,0,0}));
  end for;
 // Cells
  for j in 1:Ncell loop
    if not isSimpl[j] then // connect detailed cell ports
      connect(detailedCells[j_detl[j]].dHT_y1, dHT_y1) annotation (Line(points={{1.94,-32.2},{1.94,-40},{41,-40},{41,-105}}, color={255,127,0}));
      connect(detailedCells[j_detl[j]].dHT_y0, dHT_y0) annotation (Line(points={{3.46,-0.16},{3.46,10},{26,10},{26,80},{35,80},{35,105}}, color={255,127,0}));
    elseif isSimpl[j] then // connect simplified cell ports
      connect(simplifiedCellBlock[j_verticalBlock[j]].dHTactive_y1, dHT_y1) annotation (Line(points={{-1.5,
              14.4},{-1.5,6},{22,6},{22,-90},{42,-90},{42,-105},{41,-105}},                                                                              color={255,127,0}));
      connect(simplifiedCellBlock[j_verticalBlock[j]].dHTactive_y0, dHT_y0) annotation (Line(points={{-1.7,33.8},{-1.7,40},{35,40},{35,105}}, color={255,127,0}));
    end if;
  end for;

  connect(realTdetailedOut.y, TdetailedToFVI.T) annotation (Line(points={{106.2,-29},{80,-29},{80,-35},{71,-35}}, color={0,0,127}));

  annotation (
    Line(points={{-110,26},{-124,26}}, color={0,0,255}),
    Placement(transformation(extent={{160,-50},{180,-30}})),
    Line(points={{152.4,-40},{158,-40}}, color={0,0,127}),
    Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Rectangle(
          extent={{-124,100},{120,-100}},
          lineColor={28,108,200},
          fillColor={28,108,200},
          fillPattern=FillPattern.CrossDiag),
        Rectangle(
          extent={{-100,80},{98,-80}},
          lineColor={28,108,200},
          fillColor={0,0,99},
          fillPattern=FillPattern.Forward),
        Text(
          extent={{-94,24},{94,-26}},
          lineColor={255,255,255},
          fillColor={28,108,200},
          fillPattern=FillPattern.CrossDiag,
          textString="Simple Stack",
          textStyle={TextStyle.Bold})}),
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
<li><i>23 Dez 2021</i> by Hans Wiggenhauser</a>:<br>First release. </li>
<li><i>02 Oct 2024</i> by Anis Taissir</a>:<br> Bugfix: replaced LambdaAir by LambdaFuel in definition of LambdaFCConv in line 120. </li>
</ul>
</html>"));
end SimplifiedStack1D;

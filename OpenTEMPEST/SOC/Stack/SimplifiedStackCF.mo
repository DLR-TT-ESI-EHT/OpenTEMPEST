within OpenTEMPEST.SOC.Stack;
model SimplifiedStackCF
  import SI = Modelica.SIunits;

  parameter Integer nX(min=3)=5   "Number of control volumes in the x-direction" annotation (Dialog(tab = "General"));
  parameter Integer nY(min=3)=5   "Number of control volumes in the y-direction" annotation (Dialog(tab = "General"));
  parameter Integer Ncell(min=3)=30 "Total number of cells in the stack" annotation (Dialog(tab = "General"));
  parameter Boolean LUDS=true "Set true if Linear upwind difference wanted (more accuracy), false for Upwind difference scheme (more speed and stability)" annotation (Dialog(tab = "General"));

  // Dimensions from http://dx.doi.org/10.2139/ssrn.3987808
  parameter SI.Length lX = 0.1 "total length of solid" annotation(Dialog(tab="Dimensions"));
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
  parameter SI.Length lZsolid = 0.2e-3 "Height of interconnector" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length Bac = porAC*lY "Width of air channel without ribs" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length Bfc = porFC*lY "WIdth of fuel channel without ribs" annotation(Dialog(tab="Dimensions"));

  // Initialization
  parameter SI.Temperature TStart=973.15   "Uniform start temperature" annotation(Dialog(tab="Initialization"));
  parameter SI.AbsolutePressure pStart=101325 "Starting pressure" annotation(Dialog(tab="Initialization"));
  parameter SI.CurrentDensity JStart = 0 "Starting current density" annotation(Dialog(tab="Initialization"));
  parameter SI.MassFraction xStartFC[FCMedium.nXi] = FCMedium.reference_X "Starting mass fraction in fuel channel" annotation(Dialog(tab="Initialization"));
  parameter SI.MassFraction xStartAC[ACMedium.nXi] = ACMedium.reference_X "Starting mass fraction in air channel" annotation(Dialog(tab="Initialization"));

  // PEN parameters
  replaceable model Electrochem =
      OpenTEMPEST.SOC.Electrochem.Components.Crossflow_Electrochem
  constrainedby OpenTEMPEST.SOC.Electrochem.Components.ElectrochemBase         annotation (
  Dialog(tab="PEN"),
  Placement(transformation(extent={{78,50},{98,70}})),
  choicesAllMatching=true);
  parameter SI.ThermalConductivity kCustom_trans = 2.16 "Thermal conductivity across layers of PEN in W/mK" annotation(Dialog(tab="PEN"));
  parameter SI.ThermalConductivity kCustom_long = 2.16 "Thermal conductivity in plane of layers of PEN in W/mK" annotation(Dialog(tab="PEN"));
  parameter SI.Density rhoPEN=5900 "Density of PEN in kg/m3" annotation(Dialog(tab="PEN"));
  parameter SI.SpecificHeatCapacity cpPEN = 500 "Specific heat capacity of PEN in J/kgK" annotation(Dialog(tab="PEN"));
  parameter SI.SpectralEmissivity epsilonPEN = 0.8 "Emissivity of PEN" annotation(Dialog(tab="PEN"));
  parameter SI.ThermalConductivity lambdaPEN = 2.16 "Thermal conductivity of simplified PEN" annotation(Dialog(tab="PEN"));
  parameter SI.Temperature TPENRad=1118.15   "IC temperature assumption for radiative heat transfer" annotation(Dialog(tab="PEN"));

  // FC parameters
  replaceable package FCMedium = Medium.Fuel_CH4         annotation(Dialog(tab="Fuel Channel"));
  parameter Real porFC = 0.4 "Porosity in fuel channel" annotation(Dialog(tab="Fuel Channel"));
  parameter Real NuFCPEN = 12 "Nusselt number fuel channel on PEN side" annotation(Dialog(tab="Fuel Channel"));
  parameter Real NuFCIC = 10 "Nusselt number fuel channel on IC side" annotation(Dialog(tab="Fuel Channel"));
  parameter SI.ThermalConductivity lambdaFuel = 0.0935 "Themal conductivity of the gas in the fuel channel" annotation(Dialog(tab="Fuel Channel"));

  // AC parameters
  replaceable package ACMedium = Medium.Air_Medium         annotation(Dialog(tab="Air Channel"));
  parameter Real pDrop(max=0.99) = 0.04 "pressure drop as a factor of inlet pressure (between 0 and 0.99)" annotation(Dialog(tab="Air Channel"));
  parameter Real NuACPEN = 8 "Nusselt number air channel on PEN side" annotation(Dialog(tab="Air Channel"));
  parameter Real NuACIC = 7.5 "Nusselt number air channel on IC side" annotation(Dialog(tab="Air Channel"));
  parameter SI.ThermalConductivity lambdaAir = 72.98e-3 "Thermal conductivity of gas in air channel" annotation(Dialog(tab="Air Channel"));
  parameter Real porAC = 0.4 "Porosity in air channel" annotation(Dialog(tab="Air Channel"));

  // IC parameters
  parameter SI.Density rhoIC = 8000 "Density of IC" annotation(Dialog(tab="Interconnects"));
  parameter SI.SpecificHeatCapacity cpIC = 500 "Specific heat capacity of IC" annotation(Dialog(tab="Interconnects"));
  parameter SI.SpectralEmissivity epsilonIC = 0.1 "Emissivity of interconnects" annotation(Dialog(tab="Interconnects"));
  parameter SI.ThermalConductivity lambdaIC = 27 "Thermal conductivity of interconnects" annotation(Dialog(tab="Interconnects"));
  parameter SI.Temperature TICRad=1088.15 "IC temperature assumption for radiative heat transfer" annotation(Dialog(tab="Interconnects"));

  // Thermal conductivities for the simplified block

  // Resistances in z-direction
  parameter SI.ThermalResistance RIC = lZsolid/(lambdaIC*lX*lY) "Resistance for interconnect thermal conduction" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RPEN = lZpen/(lambdaPEN*lXpen*lYpen) "Resistance for PEN thermal conduction" annotation(Dialog(tab="Simplifications"));

  parameter SI.ThermalResistance RCondRibsAC = lZac/(lambdaIC*(lYac - Bac)*lXac) "Resistance for ribs conduction air channel side" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RConvACPEN = 2*lZac/(NuACPEN*lambdaAir*lXac*Bac) "Air channel PEN side convection resistance" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RConvACIC = 2*lZac/(NuACIC*lambdaAir*lXac*Bac) "Air channel IC side convection resistance" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RConvAC = RConvACPEN + RConvACIC "Air channel total convection resistance" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RRadAC = (epsilonPEN + epsilonIC - epsilonPEN*epsilonIC)*(TPENRad-TICRad)/(Modelica.Constants.sigma*epsilonPEN*epsilonIC*Bac*lXac*(TPENRad^4 - TICRad^4)) "Resistance for radiation air channel side" annotation(Dialog(tab="Simplifications"));

//   parameter SI.ThermalResistance RAC = 1/(1/RCondRibsAC + 1/RConvAC + 1/RRadAC) "Air channel total resistance" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RAC = 1/(1/RConvAC + 1/RRadAC) "Air channel total resistance without ribs conduction" annotation(Dialog(tab="Simplifications"));

  parameter SI.ThermalResistance RCondRibsFC = lZfc/(lambdaIC*(lYfc - Bfc)*lXfc) "Resistance for ribs conduction fuel channel side" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RConvFCPEN = 2*lZfc/(NuFCPEN*lambdaFuel*lXfc*Bfc) "Fuel channel PEN side convection" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RConvFCIC = 2*lZfc/(NuFCIC*lambdaFuel*lXfc*Bfc) "Fuel channel IC side convection" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RConvFC = RConvFCPEN + RConvFCIC "Fuel channel total convection resistance" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RRadFC = (epsilonPEN + epsilonIC - epsilonPEN*epsilonIC)*(TPENRad-TICRad)/(Modelica.Constants.sigma*epsilonPEN*epsilonIC*Bfc*lXfc*(TPENRad^4 - TICRad^4)) "Resistance for radiation fuel channel side" annotation(Dialog(tab="Simplifications"));

//   parameter SI.ThermalResistance RFC = 1/(1/RCondRibsFC + 1/RConvFC + 1/RRadFC) "Fuel channel total resistance" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RFC = 1/(1/RConvFC + 1/RRadFC) "Fuel channel total resistance without ribs conduction" annotation(Dialog(tab="Simplifications"));

  parameter SI.ThermalResistance RTrans = 2*RIC + RAC + RPEN + RFC "Cell equivalent transversal resistance" annotation(Dialog(tab="Simplifications"));

  // Resistances in longitudinal direction
  parameter SI.ThermalResistance RcondIC = lX/(lambdaIC*lZsolid*lY) "Longitudinal conductivity resistance for interconnect" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RcondPEN = lX/(lambdaPEN*lZpen*lY) "Longitudinal conductivity resistance for PEN" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RcondFC = lX/(lambdaIC*lZsolid*(lYfc - Bfc)) "Longitudinal conductivity resistance for fuel channel ribs" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RcondAC = lX/(lambdaIC*lZsolid*(lYac - Bac)) "Longitudinal conductivity resistance for air channel ribs" annotation(Dialog(tab="Simplifications"));

//   parameter SI.ThermalResistance RLong = 1/(2*1/RcondIC + 1/RcondPEN + 1/RcondFC + 1/RcondAC) "Cell equivalent longitudinal resistance" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RLong = 1/(2*1/RcondIC + 1/RcondPEN) "Cell equivalent longitudinal resistance without fuel and air channel ribs longitudinal conductivity" annotation(Dialog(tab="Simplifications"));

  // Equivalent thermal conductivities for the simplified block
  parameter SI.ThermalConductivity kCellLong = lX/(RLong*lZCell*lY) "Cell equivalent longitudinal conductivity" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalConductivity kCellTrans = lZCell/(RTrans*lX*lY) "Cell equivalent transversal conductivity" annotation(Dialog(tab="Simplifications"));

  // Simplification Parameters
  parameter SI.Length lZCell = (2*lZsolid + lZfc + lZac + lZpen) "Total height of a cell" annotation(Dialog(tab="Simplifications"));
//   parameter SI.Density rhoBlock = (rhoPEN*lZpen + 2*rhoIC*lZsolid + rhoIC*lZfc*(1-porFC) + rhoIC*lZac*(1-porAC))/lZCell "Total density of a cell with ribs in AC/FC" annotation(Dialog(tab="Simplifications"));
//   parameter SI.SpecificHeatCapacity cpBlock = (cpPEN*rhoPEN*lZpen + 2*cpIC*rhoIC*lZsolid + cpIC*rhoIC*lZfc*(1-porFC) + cpIC*rhoIC*lZac*(1-porAC))/(rhoPEN*lZpen + 2*rhoIC*lZsolid + rhoIC*lZfc*(1-porFC) + rhoIC*lZac*(1-porAC)) "Specific heat capacity of a cell" annotation(Dialog(tab="Simplifications"));

  parameter SI.Density rhoBlock = (rhoPEN*lZpen + 2*rhoIC*lZsolid)/lZCell "Total density of a cell without ribs in AC/FC" annotation(Dialog(tab="Simplifications"));
  parameter SI.SpecificHeatCapacity cpBlock = (cpPEN*rhoPEN*lZpen + 2*cpIC*rhoIC*lZsolid)/(rhoPEN*lZpen + 2*rhoIC*lZsolid) "Specific heat capacity of a cell without ribs in AC/FC" annotation(Dialog(tab="Simplifications"));

  // Simplification parameters
  parameter Integer NdetailedCell(min=1, max=Ncell)=4 "number of cells that are calculated in detail" annotation(Dialog(tab = "General"));
  parameter Integer[:] verticalBlockSize={3, 3, 2, 2, 2, 1, 1, 1, 1, 2, 2, 3, 3} "Vertical size of the simplified blocks (e.g. 3 represents three cells merged in one simplfied block). Sum needs to be Ncell-NdetailedCell" annotation (Dialog(tab = "General"));

  parameter Integer nCellVertMultHelp[nY, nX, nSimpMult] = fill(verticalBlockSize, nY, nX) "Helping parameter to build nCellVertMult";
  parameter Integer nCellVertMult[nSimpMult, nX, nY] = {{{nCellVertMultHelp[i,j,k] for i in 1:nY} for j in 1:nX} for k in 1:nSimpMult} "Vertical size of the simplified blocks information passed to all control volumes";

  parameter Integer nSimpMult = size(verticalBlockSize, 1) "number of blocks of simplified cells";
  parameter Integer nIMPs = integer((floor((Ncell-1)/10))) "number of intermediate plates";
  parameter Integer stackingHelper[Ncell] = {integer(i/(Ncell)*(Ncell-NdetailedCell) + 0.5) - integer((i - 1)/(Ncell)*(Ncell-NdetailedCell) + 0.5) for i in 1:Ncell} "define which cell is detailed (0) or simplified (1)";

protected
  parameter Integer j_simp[Ncell] = {sum(stackingHelper[1:i]) for i in 1:Ncell} "counter for simplified cells";
  parameter Integer j_detl[Ncell] = {i - j_simp[i] for i in 1:Ncell} "counter for detailed cells";
  parameter Boolean isSimpl[Ncell] = {stackingHelper[i] > 0.5 for i in 1:Ncell} "turn integer in stackingHelper to boolean";
  parameter Integer distancesCellsToCells[Ncell, Ncell] = transpose(fill(1:Ncell, Ncell)) - fill(1:Ncell, Ncell) "Distances from each cell to each other cell";
  parameter Integer j_simp_positions[Ncell] = {Modelica.Math.Vectors.find(i,j_simp) for i in 1:Ncell} "list of simplified cells and trailing zeroes";
  parameter Integer j_detl_positions[Ncell] = {Modelica.Math.Vectors.find(i,j_detl) for i in 1:Ncell} "list of detailed cells and trailing zeroes";
  parameter Integer verticalBlockRanges[nSimpMult,2] = {{j_simp_positions[sum(verticalBlockSize[1:i])]-verticalBlockSize[i]+1,j_simp_positions[sum(verticalBlockSize[1:i])]} for i in 1:size(verticalBlockSize,1)}
  "list with cells indices of limiting upper and lower cell in each simplified vertical block";
  parameter Integer j_next_detailed_cell[nSimpMult-1] = {if j_detl[verticalBlockRanges[i+1,2]]==0 then j_detl_positions[1] else j_detl_positions[max(1,j_detl[verticalBlockRanges[i+1,2]])] for i in 1:nSimpMult-1}
  "list the next detailed cell corresponding to each cell. Avoid translation error by using max(1,index), where index can be 0.";
  parameter Integer verticalBlockRangesHighActual[nSimpMult] = cat(1,{if j_detl[verticalBlockRanges[i+1,2]]==1+j_detl[verticalBlockRanges[i,2]] then j_next_detailed_cell[i]-1
                    else verticalBlockRanges[i,2] for i in 1:nSimpMult-1},{verticalBlockRanges[nSimpMult,2]})
  "corrected limiting upper cell in each simplified vertical block. If verticalBlockSize[:] is such that every detailed cell is placed exactly
   between two simplified cell vertical blocks, this is equal to verticalBlockRanges[:]. Otherwise it is adapted correspondingly.";
  parameter Integer verticalBlockRangesLowActual[nSimpMult] = cat(1,{verticalBlockRanges[1,1]},{if j_detl[verticalBlockRangesHighActual[i]]==1+j_detl[verticalBlockRangesHighActual[i-1]] then verticalBlockRangesHighActual[i-1]+2
                    else verticalBlockRangesHighActual[i-1]+(j_detl[verticalBlockRangesHighActual[i]]-j_detl[verticalBlockRangesHighActual[i-1]])+1 for i in 2:nSimpMult})
  "corrected limiting lower cell in each simplified vertical block. If verticalBlockSize[:] is such that every detailed cell is placed exactly between two simplified cell vertical blocks,
   this is equal to verticalBlockRanges[:,1]. Otherwise it is adapted correspondingly.";
  parameter Integer verticalBlockRangesActual[nSimpMult,2] = {{verticalBlockRangesLowActual[i],verticalBlockRangesHighActual[i]} for i in 1:nSimpMult}
  "corrected limiting lower cell in each simplified vertical block. If verticalBlockSize[:] is such that every detailed cell is placed exactly between two simplified cell vertical blocks,
   this is equal to verticalBlockRanges. Otherwise it is adapted correspondingly.";
  parameter Integer distancesCellsToVerticalBlocks[nSimpMult,Ncell] = {{abs(sum(distancesCellsToCells[i,l_ind] for i in verticalBlockRangesActual[j,1]:verticalBlockRangesActual[j,2])) for l_ind in 1:Ncell} for j in 1:nSimpMult}
  "helper with distances of each cell (columns) to each simplified cells vertical block (rows)";
  parameter Integer distancesHelperVerticalBlocks[nSimpMult,Ncell] = {{if isSimpl[i] then max(distancesCellsToVerticalBlocks) else distancesCellsToVerticalBlocks[j,i] for i in 1:Ncell} for j in 1:nSimpMult}
  "helper with total distance of each simplified cell block to each detailed cell, as the sum of the distances of each simplified cell contained in the block";
  parameter Integer j_corresponding[nSimpMult] = {Modelica.Math.Vectors.find(min(distancesHelperVerticalBlocks[i,:]),distancesHelperVerticalBlocks[i,:]) for i in 1:nSimpMult}
  "list detailed cells corresponding to each simplified cell vertical block";
  parameter Integer j_verticalBlock[Ncell] = {Modelica.Math.Vectors.find(1,{if verticalBlockRangesActual[j,1]-i==0 or verticalBlockRangesActual[j,2]-i==0 or (verticalBlockRangesActual[j,1]-i<0 and verticalBlockRangesActual[j,2]-i>0) then 1
                    else 0 for j in 1:size(verticalBlockRangesActual,1)}) for i in 1:Ncell}
  "index of the simplified block corresponding to each cell (0 for none)";
  parameter Integer j_SimpDet_z0[Ncell] = cat(1,{0},{sum({if -(stackingHelper[i+1]-stackingHelper[i])>0 then 1 else 0 for i in 1:j}) for j in 1:Ncell-1})
  "indexes of interfaces between bottom of detailed cell and top of previous simplified cell";
  parameter Integer j_SimpDet_z1[Ncell] = cat(1,{sum({if stackingHelper[i]-stackingHelper[i-1]>0 then 1 else 0 for i in 2:j}) for j in 2:Ncell},{sum({if stackingHelper[i]-stackingHelper[i-1]>0 then 1 else 0 for i in 2:Ncell})})
  "indexes of interfaces between top of detailed cell and bottom of next simplified cell";
  parameter Integer Ncontacts_SimpDet_z0 = max(j_SimpDet_z0) "total contact points between simplified and detailed cells at the detailed-cell vertical bottom position";
  parameter Integer Ncontacts_SimpDet_z1 = max(j_SimpDet_z1) "total contact points between simplified and detailed cells at the detailed-cell vertical top position";
  parameter Integer simpleCellHeatPortIndex[Ncell] = {j-verticalBlockRangesActual[max(1,j_verticalBlock[j]),1]+1 for j in 1:Ncell};

public
  OpenTEMPEST.SOC.Cell.CrossFlow.Cell detailedCells[NdetailedCell](
    each nX=nX,
    each nY=nY,
    each lX=lX,
    each lY=lY,
    each lXpen=lXpen,
    each lYpen=lYpen,
    each lZpen=lZpen,
    each lXac=lXac,
    each lYac=lYac,
    each lZac=lZac,
    each lXfc=lXfc,
    each lYfc=lYfc,
    each lZfc=lZfc,
    each lZsolid=lZsolid,
    each TStart=TStart,
    each pStart=pStart,
    each Jstart=JStart,
    each xStartFC=xStartFC,
    each xStartAC=xStartAC,
    redeclare model Electrochem = Electrochem,
    each kCustom_trans=kCustom_trans,
    each kCustom_long=kCustom_long,
    each rhoPEN=rhoPEN,
    each cpPEN=cpPEN,
    each epsilonPEN=epsilonPEN,
    each porFC=porFC,
    each NuFCPEN=NuFCPEN,
    each NuFCIC=NuFCIC,
    each kIC=lambdaIC,
    each cpIC=cpIC,
    each rhoIC=rhoIC,
    each porAC=porAC,
    each pDrop=pDrop,
    each NuACPEN=NuACPEN,
    each NuACIC=NuACIC,
    each epsilonIC=epsilonIC)
    annotation (Placement(transformation(extent={{-20,-66},{20,-22}})));

  SimplifiedCellBlockCF simplifiedBlock[nSimpMult](
    each nX=nX,
    each nY=nY,
    each TStart=TStart,
    nCellVertMult=verticalBlockSize,
    each lX=lX,
    each lY=lY,
    each lZpen=lZpen,
    each lZac=lZac,
    each lZfc=lZfc,
    each lZsolid=lZsolid,
    each rhoPEN=rhoPEN,
    each cpPEN=cpPEN,
    each epsilonPEN=epsilonPEN,
    each lambdaPEN=lambdaPEN,
    each TPENRad=TPENRad,
    each rhoIC=rhoIC,
    each cpIC=cpIC,
    each epsilonIC=epsilonIC,
    each lambdaIC=lambdaIC,
    each TICRad=TICRad,
    each NuACPEN=NuACPEN,
    each NuACIC=NuACIC,
    each lambdaAir=lambdaAir,
    each porAC=porAC,
    each NuFCPEN=NuFCPEN,
    each NuFCIC=NuFCIC,
    each lambdaFuel=lambdaFuel,
    each porFC=porFC,
    each RIC=RIC,
    each RPEN=RPEN,
    each RCondRibsAC=RCondRibsAC,
    each RConvACPEN=RConvACPEN,
    each RConvACIC=RConvACIC,
    each RConvAC=RConvAC,
    each RRadAC=RRadAC,
    each RAC=RAC,
    each RCondRibsFC=RCondRibsFC,
    each RConvFCPEN=RConvFCPEN,
    each RConvFCIC=RConvFCIC,
    each RConvFC=RConvFC,
    each RRadFC=RRadFC,
    each RFC=RFC,
    each RTrans=RTrans,
    each RcondIC=RcondIC,
    each RcondPEN=RcondPEN,
    each RcondFC=RcondFC,
    each RcondAC=RcondAC,
    each RLong=RLong,
    each kCellLong=kCellLong,
    each kCellTrans=kCellTrans,
    each rhoBlock=rhoBlock,
    each cpBlock=cpBlock,
    each lZCell=lZCell)
    annotation (Placement(transformation(extent={{-10,-2},{10,18}})));

  OpenTEMPEST.Heat.Solid2D plates[2 + nIMPs](
    redeclare package SolidMat = OpenTEMPEST.Solid.Material.Crofer22APU,
    each nX=nX,
    each nY=nY,
    each Tstartbar=TStart,
    each lX=lX,
    each lY=lY,
    each lZ=lZsolid)
    annotation (Placement(transformation(extent={{-10,52},{10,72}})));

  Heat.DHTVolumes2D dHT2_top(i=nX, j=nY) annotation (Placement(transformation(
          extent={{-76,92},{-62,106}}), iconTransformation(extent={{-76,92},{-62,
            106}})));
  Heat.DHTVolumes2D dHT2_bottom(i=nX, j=nY) annotation (Placement(
        transformation(extent={{-56,-146},{-42,-132}}), iconTransformation(
          extent={{-76,-146},{-62,-132}})));
  Heat.DHTVolumes2D dHT2_x0(i=Ncell, j=nY) annotation (Placement(transformation(
          extent={{-126,-38},{-112,-24}}), iconTransformation(extent={{-126,-38},
            {-112,-24}})));
  Heat.DHTVolumes2D dHT2_xN(i=Ncell, j=nY) annotation (Placement(transformation(
          extent={{114,-38},{128,-24}}), iconTransformation(extent={{114,-38},{
            128,-24}})));

  DummyManifold Fuel_Manifold(
    redeclare package Medium = Medium.Fuel_CH4,
    final nPorts_b=Ncell,
    nDummyPorts_b=Ncell - NdetailedCell)
    annotation (Placement(transformation(extent={{-66,4},{-56,30}})));

  DummyManifold Air_Manifold(
    redeclare package Medium = Medium.Air_Medium,
    final nPorts_b=Ncell,
    nDummyPorts_b=Ncell - NdetailedCell)
    annotation (Placement(transformation(extent={{-68,-96},{-56,-66}})));

  ThermoPower.Gas.FlangeA fuelIn(redeclare package Medium = Medium.Fuel_CH4)
    annotation (Placement(transformation(extent={{-130,38},{-110,58}})));
  ThermoPower.Gas.FlangeA airIn(redeclare package Medium = Medium.Air_Medium)
    annotation (Placement(transformation(extent={{-130,-86},{-110,-66}})));
  ThermoPower.Gas.FlangeB fuelOut(redeclare package Medium = Medium.Fuel_CH4)
    annotation (Placement(transformation(extent={{110,38},{130,58}})));
  ThermoPower.Gas.FlangeB airOut(redeclare package Medium = Medium.Air_Medium)
    annotation (Placement(transformation(extent={{108,-86},{128,-66}})));

  Flow.Manifold_out manifold_outFuel(redeclare package Medium = Medium.Fuel_CH4,
      nPorts_a=NdetailedCell)
    annotation (Placement(transformation(extent={{34,22},{62,50}})));

  Flow.Manifold_out manifold_outAir(redeclare package Medium =
        Medium.Air_Medium, nPorts_a=NdetailedCell)
    annotation (Placement(transformation(extent={{34,-38},{62,-8}})));

  Modelica.Electrical.Analog.Interfaces.NegativePin pin_n annotation (Placement(transformation(extent={{-130,14},
            {-110,34}}), iconTransformation(extent={{-130,14},{-110,34}})));
  Modelica.Electrical.Analog.Interfaces.PositivePin pin_p annotation (Placement(transformation(extent={{-130,
            -64},{-110,-44}})));
  Heat.DHT_DHTVolumes2D dHT_DHT2_x0(i=Ncell, j=nY)
    annotation (Placement(transformation(extent={{-62,-44},{-70,-36}})));
  Heat.DHT_DHTVolumes2D dHT_DHT2_xN(i=Ncell, j=nY)
    annotation (Placement(transformation(extent={{70,-44},{78,-36}})));

  SI.EnergyFlowRate EfSum "Check energy balance through sum of energy flows";
  SI.MassFlowRate massSum "Check mass balance through sum of mass flows";

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
    annotation (Placement(transformation(extent={{-92,34},{-72,54}})));
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
    annotation (Placement(transformation(extent={{-100,-84},{-80,-64}})));
  Modelica.Blocks.Sources.RealExpression USimp(y=-sum(simplifiedBlock[:].UCalc[:,:] .* nCellVertMult[:,:,:])/(nX*nY))
    annotation (Placement(transformation(extent={{-74,-10},{-94,10}})));
  Modelica.Electrical.Analog.Sources.SignalVoltage signalVoltage annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={-120,0})));
  Modelica.Blocks.Sources.RealExpression realTdetailedOut[NdetailedCell,nY](each y=
        sum(detailedCells.fuelChannel.T[:, nY])/(nX*NdetailedCell))
    annotation (Placement(transformation(extent={{130,-114},{110,-94}})));

  Heat.DHT_DHTVolumes2D dHT_DHT2_yN(i=Ncell, j=nX) annotation (Placement(
        transformation(
        extent={{4,-4},{-4,4}},
        rotation=270,
        origin={56,74})));
  Heat.DHT_DHTVolumes2D dHT_DHT2_y0(i=Ncell, j=nX) annotation (Placement(
        transformation(
        extent={{4,-4},{-4,4}},
        rotation=90,
        origin={56,-96})));
  Heat.DHTVolumes2D dHT2_yN(i=Ncell, j=nX) annotation (Placement(transformation(
          extent={{52,92},{66,106}}), iconTransformation(extent={{52,92},{66,
            106}})));
  Heat.DHTVolumes2D dHT2_y0(i=Ncell, j=nX) annotation (Placement(transformation(
          extent={{52,-146},{66,-132}}), iconTransformation(extent={{52,-146},{
            66,-132}})));
  FlowInterpolator1D flowInterpolator1D(
    nCell=Ncell,
    N=nY,
    nSimplified=Ncell - NdetailedCell,
    nNonUnitOrSimpCells=0,
    isRedu=isSimpl,
    TStart=TStart,
    pStart=pStart)
    annotation (Placement(transformation(extent={{74,-6},{98,18}})));
  OpenTEMPEST.Heat.PrescribedTemperature1D TdetailedToFVI1D[NdetailedCell](
      each N=nY) annotation (Placement(transformation(
        extent={{-8,-8},{8,8}},
        rotation=180,
        origin={84,-96})));

  Heat.DHT_DHTVolumes2D dHT_DHT2Plates_y0(i=2 + nIMPs, j=nX) annotation (
      Placement(transformation(
        extent={{4,-4},{-4,4}},
        rotation=90,
        origin={18,42})));
  Heat.DHT_DHTVolumes2D dHT_DHT2Plates_x0(i=2 + nIMPs, j=nY) annotation (
      Placement(transformation(
        extent={{4,-4},{-4,4}},
        rotation=0,
        origin={-54,62})));
  Heat.DHT_DHTVolumes2D dHT_DHT2Plates_yN(i=2 + nIMPs, j=nX) annotation (
      Placement(transformation(
        extent={{4,-4},{-4,4}},
        rotation=270,
        origin={18,80})));
  Heat.DHT_DHTVolumes2D dHT_DHT2Plates_xN(i=2 + nIMPs, j=nY) annotation (
      Placement(transformation(
        extent={{4,-4},{-4,4}},
        rotation=180,
        origin={36,62})));
  Heat.DHTVolumes2D dHT2Plates_yN(i=2 + nIMPs, j=nX) annotation (Placement(
        transformation(extent={{26,92},{40,106}}), iconTransformation(extent={{
            30,92},{44,106}})));
  Heat.DHTVolumes2D dHT2Plates_y0(i=2 + nIMPs, j=nX) annotation (Placement(
        transformation(extent={{26,-146},{40,-132}}), iconTransformation(extent
          ={{30,-146},{44,-132}})));
  Heat.DHTVolumes2D dHT2Plates_x0(i=2 + nIMPs, j=nX) annotation (Placement(
        transformation(extent={{-126,64},{-112,78}}), iconTransformation(extent
          ={{-126,72},{-112,86}})));
  Heat.DHTVolumes2D dHT2Plates_xN(i=2 + nIMPs, j=nX) annotation (Placement(
        transformation(extent={{114,64},{128,78}}), iconTransformation(extent={
            {112,72},{126,86}})));
equation
  assert(not Modelica.Math.BooleanVectors.anyTrue(Modelica.Math.BooleanVectors.anyTrue({{j_detl_positions[i]-verticalBlockRangesActual[j,1]>=0 and verticalBlockRangesActual[j,2]-j_detl_positions[i]>=0 for i in 1:size(j_detl_positions,1)} for j in 1:size(verticalBlockRangesActual,1)})),
  "error: some vertical blocks are overlapping with detailed cell positions.");
  assert(sum(verticalBlockSize)==Ncell-NdetailedCell,"Vertical block sizes must total up to Ncell-NdetailedCell");
  assert(Modelica.Math.BooleanVectors.anyTrue(Modelica.Math.BooleanVectors.anyTrue({{verticalBlockRangesLowActual[i]==j*10 or verticalBlockRangesHighActual[i]==j*10 for i in 1:size(verticalBlockRangesLowActual,1)} for j in 1:nIMPs})) or Ncell <=10,"error: some vertical blocks are overlapping with intermediate plate positions");

  // Checking energy balance
  EfSum = fuelIn.m_flow * inStream(fuelIn.h_outflow) + fuelOut.m_flow * fuelOut.h_outflow
        + airIn.m_flow * inStream(airIn.h_outflow) + airOut.m_flow * airOut.h_outflow
        + (pin_p.v - pin_n.v) * pin_n.i
        + sum(dHT2_x0.Q[:,:]) + sum(dHT2_xN.Q[:,:])
        + sum(dHT2_top.Q[:,:]) + sum(dHT2_bottom.Q[:,:])
        + sum(dHT2_y0.Q[:,:]) + sum(dHT2_yN.Q[:,:])
        + sum(dHT2Plates_x0.Q[:,:]) + sum(dHT2Plates_xN.Q[:,:])
        + sum(dHT2Plates_y0.Q[:,:]) + sum(dHT2Plates_yN.Q[:,:]);

  // Checking mass balance
  massSum = fuelIn.m_flow + fuelOut.m_flow + airIn.m_flow + airOut.m_flow;

  // Voltage and temperature simplification
  flowInterpolator1D.TOutAirSimp = sum(sum(simplifiedBlock[j].simpFlowInfoPortOut[:,nY].TAir)*verticalBlockSize[j] for j in 1:nSimpMult)/(Ncell*nX) + sum(detailedCells.airChannel.Gas[nY,:].T)/(Ncell*nX);
  flowInterpolator1D.TOutFuelSimp = sum(sum(simplifiedBlock[j].simpFlowInfoPortOut[nX,:].TFuel)*verticalBlockSize[j] for j in 1:nSimpMult)/(Ncell*nY) + sum(detailedCells.fuelChannel.Gas[nX,:].T)/(Ncell*nY);

  // Pass information to simplifiedBlock
  for i in 1:nSimpMult loop
    for m_idx in 1:nX loop
      for n_idx in 1:nY loop
        // here: j_detl[l_corresponding[j]] is the number of detailed cells corresponding to the detailed cell index, corresponding to the current simplified cells block
        // reactionHeatInfoPort
        simplifiedBlock[i].reactionHeatInfoPort[m_idx, n_idx].H_r = 0;
        simplifiedBlock[i].reactionHeatInfoPort[m_idx, n_idx].LHVFlow = 0;
        simplifiedBlock[i].reactionHeatInfoPort[m_idx, n_idx].Vop = detailedCells[j_detl[j_corresponding[i]]].pen.electrochem[m_idx, n_idx].Uop;
        simplifiedBlock[i].reactionHeatInfoPort[m_idx, n_idx].Vid = detailedCells[j_detl[j_corresponding[i]]].pen.electrochem[m_idx, n_idx].Uideal;
        simplifiedBlock[i].reactionHeatInfoPort[m_idx, n_idx].I = detailedCells[j_detl[j_corresponding[i]]].pen.Iv[m_idx, n_idx];
        simplifiedBlock[i].reactionHeatInfoPort[m_idx, n_idx].TPEN = detailedCells[j_detl[j_corresponding[i]]].pen.electrochem[m_idx, n_idx].Tpen;

        // simpFlowInfoPortIn
        // Fuel channel
        simplifiedBlock[i].simpFlowInfoPortIn[m_idx, n_idx].mfFuel = +detailedCells[j_detl[j_corresponding[i]]].fuelChannel.mfv[m_idx, n_idx];
        simplifiedBlock[i].simpFlowInfoPortIn[m_idx, n_idx].HfFuelRef = +detailedCells[j_detl[j_corresponding[i]]].fuelChannel.mfv[m_idx, n_idx] * detailedCells[j_detl[j_corresponding[i]]].fuelChannel.hv[m_idx, n_idx];

        if m_idx==1 then
          simplifiedBlock[i].simpFlowInfoPortIn[m_idx, n_idx].cpFuel= sensCpFuel.cp;
          simplifiedBlock[i].simpFlowInfoPortIn[m_idx, n_idx].TFuelRef = detailedCells[j_detl[j_corresponding[i]]].fuelChannel.T[1, n_idx];
          simplifiedBlock[i].simpFlowInfoPortIn[m_idx, n_idx].TFuel = detailedCells[j_detl[j_corresponding[i]]].fuelChannel.T[1, n_idx];
        elseif m_idx>1 then
          simplifiedBlock[i].simpFlowInfoPortIn[m_idx, n_idx].cpFuel= FCMedium.specificHeatCapacityCp(detailedCells[j_detl[j_corresponding[i]]].fuelChannel.Gas[m_idx-1, n_idx].state);
          simplifiedBlock[i].simpFlowInfoPortIn[m_idx, n_idx].TFuelRef = detailedCells[j_detl[j_corresponding[i]]].fuelChannel.Gas[m_idx-1, n_idx].T;
          simplifiedBlock[i].simpFlowInfoPortIn[m_idx, n_idx].TFuel = simplifiedBlock[i].simpFlowInfoPortOut[m_idx-1, n_idx].TFuel;
        end if;

        // Air channel
        simplifiedBlock[i].simpFlowInfoPortIn[m_idx, n_idx].mfAir = +detailedCells[j_detl[j_corresponding[i]]].airChannel.mfv[n_idx, m_idx];
        simplifiedBlock[i].simpFlowInfoPortIn[m_idx, n_idx].HfAirRef = +detailedCells[j_detl[j_corresponding[i]]].airChannel.mfv[n_idx, m_idx] * detailedCells[j_detl[j_corresponding[i]]].airChannel.hv[n_idx, m_idx];

        if n_idx==1 then
          simplifiedBlock[i].simpFlowInfoPortIn[m_idx, n_idx].cpAir = sensCpAir.cp;
          simplifiedBlock[i].simpFlowInfoPortIn[m_idx, n_idx].TAirRef = detailedCells[j_detl[j_corresponding[i]]].airChannel.T[1, m_idx];
          simplifiedBlock[i].simpFlowInfoPortIn[m_idx, n_idx].TAir = detailedCells[j_detl[j_corresponding[i]]].airChannel.T[1, m_idx];
        elseif n_idx>1 then
          simplifiedBlock[i].simpFlowInfoPortIn[m_idx, n_idx].cpAir = ACMedium.specificHeatCapacityCp(detailedCells[j_detl[j_corresponding[i]]].airChannel.Gas[n_idx-1, m_idx].state);
          simplifiedBlock[i].simpFlowInfoPortIn[m_idx, n_idx].TAirRef = detailedCells[j_detl[j_corresponding[i]]].airChannel.Gas[n_idx-1, m_idx].T;
          simplifiedBlock[i].simpFlowInfoPortIn[m_idx, n_idx].TAir = simplifiedBlock[i].simpFlowInfoPortOut[m_idx, n_idx-1].TAir;
        end if;

        // simpFlowInfoPortOut
        // Fuel channel
        simplifiedBlock[i].simpFlowInfoPortOut[m_idx, n_idx].mfFuel = -detailedCells[j_detl[j_corresponding[i]]].fuelChannel.mfv[m_idx+1, n_idx];
        simplifiedBlock[i].simpFlowInfoPortOut[m_idx, n_idx].cpFuel= FCMedium.specificHeatCapacityCp(detailedCells[j_detl[j_corresponding[i]]].fuelChannel.Gas[m_idx, n_idx].state);
        simplifiedBlock[i].simpFlowInfoPortOut[m_idx, n_idx].HfFuelRef = -detailedCells[j_detl[j_corresponding[i]]].fuelChannel.mfv[m_idx+1, n_idx] * detailedCells[j_detl[j_corresponding[i]]].fuelChannel.hv[m_idx+1, n_idx];

        if m_idx<nX then
          simplifiedBlock[i].simpFlowInfoPortOut[m_idx, n_idx].TFuelRef = detailedCells[j_detl[j_corresponding[i]]].fuelChannel.Gas[m_idx, n_idx].T;
          simplifiedBlock[i].simpFlowInfoPortOut[m_idx, n_idx].TFuel = simplifiedBlock[i].cellBlock.T[m_idx, n_idx];
        elseif m_idx==nX then
          simplifiedBlock[i].simpFlowInfoPortOut[m_idx, n_idx].TFuelRef = detailedCells[j_detl[j_corresponding[i]]].fuelChannel.T[nX, n_idx];
          simplifiedBlock[i].simpFlowInfoPortOut[m_idx, n_idx].TFuel = simplifiedBlock[i].cellBlock.dhT_xN.T[n_idx];
        end if;

        // Air channel
        simplifiedBlock[i].simpFlowInfoPortOut[m_idx, n_idx].mfAir = -detailedCells[j_detl[j_corresponding[i]]].airChannel.mfv[n_idx+1, m_idx];
        simplifiedBlock[i].simpFlowInfoPortOut[m_idx, n_idx].cpAir= ACMedium.specificHeatCapacityCp(detailedCells[j_detl[j_corresponding[i]]].airChannel.Gas[n_idx, m_idx].state);
        simplifiedBlock[i].simpFlowInfoPortOut[m_idx, n_idx].HfAirRef = -detailedCells[j_detl[j_corresponding[i]]].airChannel.mfv[n_idx+1, m_idx] * detailedCells[j_detl[j_corresponding[i]]].airChannel.hv[n_idx+1, m_idx];

        if n_idx<nY then
          simplifiedBlock[i].simpFlowInfoPortOut[m_idx, n_idx].TAirRef = detailedCells[j_detl[j_corresponding[i]]].airChannel.Gas[n_idx, m_idx].T;
          simplifiedBlock[i].simpFlowInfoPortOut[m_idx, n_idx].TAir = simplifiedBlock[i].cellBlock.T[m_idx, n_idx];
        elseif n_idx==nY then
          simplifiedBlock[i].simpFlowInfoPortOut[m_idx, n_idx].TAirRef = detailedCells[j_detl[j_corresponding[i]]].airChannel.T[nY, m_idx];
          simplifiedBlock[i].simpFlowInfoPortOut[m_idx, n_idx].TAir = simplifiedBlock[i].cellBlock.dhT_yN.T[m_idx];
        end if;

      end for;
    end for;
  end for;

  //Electrical connections - External
  connect(USimp.y, signalVoltage.v) annotation (Line(points={{-95,0},{-96.5,0},{
          -96.5,-2.33147e-15},{-108,-2.33147e-15}},
                            color={0,0,127}));
  connect(signalVoltage.n, detailedCells[1].pin_n) annotation (Line(points={{-120,
          -10},{-120,-18},{-72,-18},{-72,-32},{-30,-32},{-30,-32},{-19.2,-32}},
                                                                         color={
          0,0,255}));
  connect(pin_p, detailedCells[NdetailedCell].pin_p) annotation (Line(points={{-120,
          -54},{-38,-54},{-38,-40},{-19.2,-40}}, color={0,0,255}));
  connect(pin_n, signalVoltage.p) annotation (Line(points={{-120,24},{-120,10}},           color={0,0,255}));

  //Electrical connections - Inter cell (detailed cells)
  for j in 1:Ncell loop
    if not isSimpl[j] then
      // here: j_detl[j] is the detailed cell index corresponding to the complete cells index.
      if NdetailedCell > 1 and j_detl[j] <= NdetailedCell-1 then // More than one detailed cell in between
        // Electrical connections (Series/Internal)
        connect(detailedCells[j_detl[j]].pin_p, detailedCells[j_detl[j] + 1].pin_n);
      end if;
    end if;
  end for;

  // Connect z-direction dHT2 ports
  connect(plates[nIMPs + 2].dhT2_z1, dHT2_top) annotation (Line(points={{-7,69},
          {-8,69},{-8,78},{-36,78},{-36,99},{-69,99}},                                             color={0,0,0}));
  connect(plates[1].dhT2_z0, dHT2_bottom) annotation (Line(points={{-7,55},{-32,
          55},{-32,16},{-36,16},{-36,-128},{-49,-128},{-49,-139}},                                                           color={0,0,0}));

  // Connect endplates
  // Top plate
  if isSimpl[Ncell] then
    // Connect simplified cell to top plate
    connect(plates[nIMPs + 2].dhT2_z0, simplifiedBlock[j_verticalBlock[Ncell]].dHT2_z1);
  else
    // Connect detailed cell to top plate
    connect(plates[nIMPs + 2].dhT2_z0, detailedCells[j_detl[Ncell]].dHT2_z1);
  end if;

  // Bottom plate
  if isSimpl[1] then
    // Connect simplified cell to bottom plate below
    connect(plates[1].dhT2_z1, simplifiedBlock[j_verticalBlock[1]].dHT2_z0);
  else
    // Connect detailed cell to bottom plate below
    connect(plates[1].dhT2_z1, detailedCells[1].dHT2_z0);
  end if;

  // Consider possible cases: detl[j] to simp[j+1], simp[j] to detl[j+1], simp[1] to bottom plate, etc...
  for j in 1:Ncell-1 loop
    if isSimpl[j] then
      if mod(j, 10) == 0 and j > 10 then
        // Connect simplified block to intermediate plate above
        connect(simplifiedBlock[j_verticalBlock[j]].dHT2_z1, plates[integer(1+j/10)].dhT2_z0);
        if isSimpl[j + 1] then
          // Connect intermediate plate to simplified block above
          connect(plates[integer(1+j/10)].dhT2_z1, simplifiedBlock[j_verticalBlock[j+1]].dHT2_z0);
        else
          // Connect intermediate plate to detailed cell above
          connect(plates[integer(1+j/10)].dhT2_z1, detailedCells[j_detl[j] + 1].dHT2_z0);
        end if;
      elseif not isSimpl[j+1] and verticalBlockRangesActual[max(1, j_verticalBlock[j]), 2] - j == 0 then
        // Only connect upwards if simplified cell j is at the top of its vertical block (Simplified cells in between are not separate objects)
        // Connect simplified cell to detailed cell above (through interface)
        connect(simplifiedBlock[j_verticalBlock[j]].dHT2_z1, detailedCells[j_detl[j+1]].dHT2_z0);
      elseif verticalBlockRangesActual[max(1, j_verticalBlock[j]), 2] - j == 0 then
        // Only connect upwards if simplified cell j is at the top of its vertical block (simplified cells in between are not separate objects)
        // Connect simplified cell to simplified cell above
        connect(simplifiedBlock[j_verticalBlock[j]].dHT2_z1, simplifiedBlock[j_verticalBlock[j+1]].dHT2_z0);
      end if;
    elseif not isSimpl[j] then
      if mod(j, 10) == 0 and j > 10 then
        // Connect detailed cell to intermediate plate above
        connect(detailedCells[j_detl[j]].dHT2_z1, plates[integer(1+j/10)].dhT2_z0);
        if isSimpl[j+1] then
          // Connect intermediate plate to simplified block above
          connect(plates[integer(1+j/10)].dhT2_z1, simplifiedBlock[j_verticalBlock[j+1]].dHT2_z0);
        else
          // Connect intermediate cell to detailed cell above
          connect(plates[integer(1+j/10)].dhT2_z1, detailedCells[j_detl[j]+1].dHT2_z0);
        end if;
      elseif isSimpl[j+1] then
        // Connect detailed cell to simplified block above
        connect(detailedCells[j_detl[j]].dHT2_z1, simplifiedBlock[j_verticalBlock[j+1]].dHT2_z0);
      else
        //Connect detailed cell to detailed cell above
        connect(detailedCells[j_detl[j]].dHT2_z1, detailedCells[j_detl[j]+1].dHT2_z0);
      end if;
    end if;
  end for;

  // Gas flow connections
  connect(fuelIn, sensCpFuel.inlet) annotation (Line(points={{-120,48},{-96,48},{-96,40},{-88,40}},color={159,159,223}));
  connect(airIn, sensCpAir.inlet) annotation (Line(points={{-120,-76},{-120,-78},{-96,-78}},color={159,159,223}));
  connect(sensCpFuel.outlet, Fuel_Manifold.port_a) annotation (Line(points={{-76,40},{-68,40},{-68,30},{-70,30},{-70,17},{-66,17}},color={159,159,223}));
  connect(sensCpAir.outlet, Air_Manifold.port_a) annotation (Line(points={{-84,-78},
          {-84,-81},{-68,-81}},                                                                           color={159,159,223}));
  connect(Fuel_Manifold.ports_b[1:NdetailedCell], detailedCells.fuelIn) annotation (Line(points={{-56,17},
          {-42,17},{-42,-26.4},{-20.4,-26.4}},                                                                                                    color={0,127,255}));
  connect(Air_Manifold.ports_b[1:NdetailedCell], detailedCells.airIn) annotation (Line(points={{-56,-81},
          {-56,-82},{-30,-82},{-30,-62},{-20,-62}},                                                                                                 color={0,127,255}));
  connect(detailedCells.fuelOut, manifold_outFuel.ports_a[1:NdetailedCell]) annotation(Line(points={{20,-30},
          {20,36},{44.36,36}},                                                                                                    color={159,159,223}));
  connect(detailedCells.airout, manifold_outAir.ports_a[1:NdetailedCell]) annotation (Line(points={{20,-54},{30,-54},{30,-32},{28,-32},{28,-23},{44.36,-23}},color={159,159,223}));

  // x-direction
  // External x-direction connections
  connect(dHT_DHT2_xN.DHT2D_port,dHT2_xN)  annotation (Line(points={{78.4,-40},
          {108,-40},{108,-31},{121,-31}},
                                        color={0,0,0}));
  connect(dHT_DHT2_x0.DHT2D_port,dHT2_x0)  annotation (Line(points={{-70.4,-40},
          {-106,-40},{-106,-31},{-119,-31}},
                                          color={0,0,0}));

  // Connect simplified cells (externally) in x-direction
  for j in 1:Ncell loop
    if not isSimpl[j] then
      // Connect detailed cell inlet and outlet port in flow direction
      connect(detailedCells[j_detl[j]].dHT_x0, dHT_DHT2_x0.DHT_port[j]) annotation (Line(points={{-18.8,
              -50},{-52,-50},{-52,-40},{-61.2,-40}},    color={255,127,0}));
      connect(detailedCells[j_detl[j]].dHT_xN, dHT_DHT2_xN.DHT_port[j]) annotation (Line(points={{18.6,
              -42},{44.6,-42},{44.6,-40},{69.2,-40}},                                                                                    color={255,127,0}));
      connect(flowInterpolator1D.dHT[j], TdetailedToFVI1D[j_detl[j]].port);
    elseif isSimpl[j] then
      // Connect simplified cell ports. Blocks with simplified cells get connected to all corresponding ports (simplified cells lumped)
      connect(simplifiedBlock[j_verticalBlock[j]].dHT_x0[simpleCellHeatPortIndex[j]], dHT_DHT2_x0.DHT_port[j]) annotation (Line(points={{-10.2,8},
              {-50,8},{-50,-40},{-61.2,-40}}, color={255,127,0}));
      connect(simplifiedBlock[j_verticalBlock[j]].dHT_xN[simpleCellHeatPortIndex[j]], dHT_DHT2_xN.DHT_port[j]) annotation (Line(points={{10.2,8},
              {28,8},{28,-40},{69.2,-40}},                   color={255,127,0}));
      connect(simplifiedBlock[j_verticalBlock[j]].dHT_xN[simpleCellHeatPortIndex[j]], flowInterpolator1D.dHT[j]);
    end if;
  end for;

  // Plates x-direction
  connect(plates.dhT_x0, dHT_DHT2Plates_x0.DHT_port) annotation (Line(points={{-11,62},
          {-49.2,62}},                       color={255,127,0}));
  connect(dHT_DHT2Plates_x0.DHT2D_port, dHT2Plates_x0) annotation (Line(points={{-58.4,
          62},{-108,62},{-108,71},{-119,71}},        color={0,0,0}));
  connect(dHT_DHT2Plates_xN.DHT2D_port, dHT2Plates_xN) annotation (Line(points={{40.4,62},
          {110,62},{110,71},{121,71}},           color={0,0,0}));
  connect(plates.dhT_xN, dHT_DHT2Plates_xN.DHT_port)
    annotation (Line(points={{11,62},{31.2,62}}, color={255,127,0}));

  // Plates y - direction
  connect(plates.dhT_y0, dHT_DHT2Plates_y0.DHT_port) annotation (Line(points={{7,55},{
          7,50.9},{18,50.9},{18,46.8}},    color={255,127,0}));
  connect(plates.dhT_yN, dHT_DHT2Plates_yN.DHT_port)
    annotation (Line(points={{7,69},{18,69},{18,75.2}}, color={255,127,0}));
  connect(dHT_DHT2Plates_y0.DHT2D_port, dHT2Plates_y0) annotation (Line(points={{18,37.6},
          {18,22},{22,22},{22,10},{30,10},{30,-42},{28,-42},{28,-132},{33,-132},
          {33,-139}},                               color={0,0,0}));
  connect(dHT_DHT2Plates_yN.DHT2D_port, dHT2Plates_yN) annotation (Line(points={{18,84.4},
          {18,92},{33,92},{33,99}},           color={0,0,0}));

  // Cells y- direction
  for j in 1:Ncell loop
    if not isSimpl[j] then
      // Connect detailed cell inlet and outlet port in flow direction
      connect(detailedCells[j_detl[j]].dHT_y0, dHT_DHT2_y0.DHT_port[j]);
      connect(detailedCells[j_detl[j]].dHT_yN, dHT_DHT2_yN.DHT_port[j]);
    elseif isSimpl[j] then
      // Connect simplified cell ports. Blocks with simplified cells get connected to all corresponding ports (simplified cells lumped)
      connect(simplifiedBlock[j_verticalBlock[j]].dHT_y0, dHT_DHT2_y0.DHT_port[j]);
      connect(simplifiedBlock[j_verticalBlock[j]].dHT_yN, dHT_DHT2_yN.DHT_port[j]);
    end if;
  end for;

  // External y-direction
  connect(dHT_DHT2_yN.DHT2D_port, dHT2_yN) annotation (Line(points={{56,78.4},{
          56,99},{59,99}},      color={0,0,0}));
  connect(dHT_DHT2_y0.DHT2D_port, dHT2_y0) annotation (Line(points={{56,-100.4},
          {56,-139},{59,-139}},                                                                       color={0,0,0}));

  // Flow interpolator connections
  connect(flowInterpolator1D.fuelExit, fuelOut) annotation (Line(points={{95.6,10.8},
          {120,10.8},{120,48}}, color={159,159,223}));
  connect(flowInterpolator1D.airExit, airOut) annotation (Line(points={{95.6,1.2},
          {96,1.2},{96,-76},{118,-76}}, color={159,159,223}));
  connect(manifold_outFuel.port_b, flowInterpolator1D.fuelInlet) annotation (
      Line(points={{55,36},{54,36},{54,10.8},{74,10.8}}, color={0,127,255}));
  connect(manifold_outAir.port_b, flowInterpolator1D.airInlet) annotation (Line(
        points={{55,-23},{60,-23},{60,1.2},{74,1.2}}, color={0,127,255}));

  // Detailed cells input temperature
  connect(realTdetailedOut.y, TdetailedToFVI1D.T) annotation (Line(points={{109,
          -104},{100,-104},{100,-99.2},{84,-99.2}},
                                                  color={0,0,127}));

 annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-120,-140},
            {120,100}}),                                       graphics={ Rectangle(
          extent={{-120,100},{120,-140}},
          lineColor={28,108,200},
          fillColor={0,127,127},
          fillPattern=FillPattern.CrossDiag,
          lineThickness=0.5)}),                                 Diagram(
        coordinateSystem(preserveAspectRatio=false, extent={{-120,-140},{120,100}})),
    Documentation(info="<html>
<h2>SimplifiedStackCF</h2>

<p>
This model represents a crossflow fuel cell stack with a simplified approach. Out of a total of 
<code>Ncell</code> cells, <code>NdetailedCell</code> are modeled in full detail, while the remaining 
cells are approximated using simplified heat transfer blocks. The model includes thermal, 
fluid, and electrical interactions and supports interpolation of flow and voltage between 
detailed and simplified cells.
</p>
</html>"));
end SimplifiedStackCF;

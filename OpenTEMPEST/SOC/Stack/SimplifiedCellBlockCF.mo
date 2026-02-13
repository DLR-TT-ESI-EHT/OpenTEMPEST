within OpenTEMPEST.SOC.Stack;
model SimplifiedCellBlockCF
  "Model for CF simple cell block. Solid 2D block that represents the block of simplified cells with equivalent conductivities."
  import SI = Modelica.SIunits;
  parameter Integer nX=5    "number of CVs in the x-direction";
  parameter Integer nY=5    "number of CVs in the y-direction";
  parameter Integer nCellVertMult=1 "Number of cells this object is representing in z direction";
  parameter SI.Temperature TStart=750 "Starting temperature";

  // Dimensions from http://dx.doi.org/10.2139/ssrn.3987808
  parameter SI.Length lX = 0.1 "Length of unit cell" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lY = 0.1 "Width of unit cell" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lXpen = lX "Length of PEN" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lYpen = lY "Width of PEN" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lZpen = 3.425e-4 "Thickness of the PEN" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lXac = lX "Length of air channel" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lYac = lY "Width of air channel with ribs" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lZac = 1e-3 "Height of the air flow channel" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lXfc = lX "Total length of fuel channel" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lYfc = lY "Total width of all the fuel channel" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lZfc = 1e-3 "Height of the fuel flow channel" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lZsolid = 0.2e-3 "Thickness of IC unit" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length Bac = porAC*lY "Width of air channel without ribs" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length Bfc = porFC*lY "Width of fuel channel without ribs" annotation(Dialog(tab="Dimensions"));

  // PEN parameters
  parameter SI.Density rhoPEN=5900 "Density of PEN in kg/m3" annotation(Dialog(tab="PEN"));
  parameter SI.SpecificHeatCapacity cpPEN = 500 "Specific heat capacity of PEN in J/kgK" annotation(Dialog(tab="PEN"));
  parameter SI.SpectralEmissivity epsilonPEN = 0.8 "Emissivity of PEN" annotation(Dialog(tab="PEN"));
  parameter SI.ThermalConductivity lambdaPEN = 2.16 "Thermal conductivity of simplified PEN" annotation (Dialog(tab="PEN"));
  parameter SI.Temperature TPENRad=1118.15   "IC temperature assumption for radiative heat transfer" annotation(Dialog(tab="PEN"));

  // FC parameters
  parameter Real NuFCPEN = 12 "Nusselt number fuel channel on PEN side" annotation(Dialog(tab="Fuel Channel"));
  parameter Real NuFCIC = 9.86 "Nusselt number fuel channel on IC side" annotation(Dialog(tab="Fuel Channel"));
  parameter SI.ThermalConductivity lambdaFuel = 0.0935 "Thermal Conductivity of gas in fuel channel" annotation(Dialog(tab="Fuel Channel"));
  parameter Real porFC = 0.8 "Fuel channel porosity" annotation(Dialog(tab="Fuel Channel"));

  // AC parameters
  parameter Real pDrop(max=0.99) = 0.04 "Pressure drop as a factor of inlet pressure (between 0 and 0.99)" annotation(Dialog(tab="Air Channel"));
  parameter Real NuACPEN = 8.235 "Nusselt number air channel on PEN side" annotation(Dialog(tab="Air Channel"));
  parameter Real NuACIC = 7.54 "Nusselt number air channel on IC side" annotation(Dialog(tab="Air Channel"));
  parameter SI.ThermalConductivity lambdaAir = 72.98e-3 "Thermal conductivity of gas in air channel" annotation(Dialog(tab="Air Channel"));
  parameter Real porAC = 0.8 "Channel porosity" annotation(Dialog(tab="Air Channel"));

  // IC parameters
  parameter SI.Density rhoIC = 7700   "Density of IC" annotation(Dialog(tab="Interconnects"));
  parameter SI.SpecificHeatCapacity cpIC = 663.65 "Specific heat capacity of IC" annotation(Dialog(tab="Interconnects"));
  parameter SI.SpectralEmissivity epsilonIC = 0.1 "Emissivity of surface of IC" annotation(Dialog(tab="Interconnects"));
  parameter SI.ThermalConductivity lambdaIC = 23.77 "Thermal conductivity of interconnects" annotation(Dialog(tab="Interconnects"));
  parameter SI.Temperature TICRad=1088.15        "IC temperature assumption for radiative heat transfer" annotation(Dialog(tab="Interconnects"));

  // Thermal conductivities for the simplified block

  // Resistances in z-direction
  parameter SI.ThermalResistance RIC = lZsolid/(lambdaIC*lX*lY) "Resistance for interconnect thermal conduction" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RPEN = lZpen/(lambdaPEN*lXpen*lYpen) "Resistance for PEN thermal conduction" annotation(Dialog(tab="Simplifications"));

  parameter SI.ThermalResistance RCondRibsAC = lZac/(lambdaIC*(lYac - Bac)*lXac) "Resistance for ribs conduction air channel side" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RConvACPEN = 2*lZac/(NuACPEN*lambdaAir*lXac*Bac) "Air channel PEN side convection resistance" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RConvACIC = 2*lZac/(NuACIC*lambdaAir*lXac*Bac) "Air channel IC side convection resistance" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RConvAC = RConvACPEN + RConvACIC "Air channel total convection resistance" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RRadAC = (epsilonPEN + epsilonIC - epsilonPEN*epsilonIC)*(TPENRad-TICRad)/(Modelica.Constants.sigma*epsilonPEN*epsilonIC*Bac*lXac*(TPENRad^4 - TICRad^4)) "Resistance for radiation air channel side" annotation(Dialog(tab="Simplifications"));

  parameter SI.ThermalResistance RAC = 1/(1/RConvAC + 1/RRadAC) "Air channel total resistance without ribs conduction" annotation(Dialog(tab="Simplifications"));

  parameter SI.ThermalResistance RCondRibsFC = lZfc/(lambdaIC*(lYfc - Bfc)*lXfc) "Resistance for ribs conduction fuel channel side" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RConvFCPEN = 2*lZfc/(NuFCPEN*lambdaFuel*lXfc*Bfc) "Fuel channel PEN side convection" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RConvFCIC = 2*lZfc/(NuFCIC*lambdaFuel*lXfc*Bfc) "Fuel channel IC side convection" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RConvFC = RConvFCPEN + RConvFCIC "Fuel channel total convection resistance" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RRadFC = (epsilonPEN + epsilonIC - epsilonPEN*epsilonIC)*(TPENRad-TICRad)/(Modelica.Constants.sigma*epsilonPEN*epsilonIC*Bfc*lXfc*(TPENRad^4 - TICRad^4)) "Resistance for radiation fuel channel side" annotation(Dialog(tab="Simplifications"));

  parameter SI.ThermalResistance RFC = 1/(1/RConvFC + 1/RRadFC) "Fuel channel total resistance without ribs conduction" annotation(Dialog(tab="Simplifications"));

  parameter SI.ThermalResistance RTrans = 2*RIC + RAC + RPEN + RFC "Cell equivalent transversal resistance" annotation(Dialog(tab="Simplifications"));

  // Resistances in longitudinal direction
  parameter SI.ThermalResistance RcondIC = lX/(lambdaIC*lZsolid*lY) "Longitudinal conductivity resistance for interconnect" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RcondPEN = lX/(lambdaPEN*lZpen*lY) "Longitudinal conductivity resistance for PEN" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RcondFC = lX/(lambdaIC*lZsolid*(lYfc - Bfc)) "Longitudinal conductivity resistance for fuel channel ribs" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalResistance RcondAC = lX/(lambdaIC*lZsolid*(lYac - Bac)) "Longitudinal conductivity resistance for air channel ribs" annotation(Dialog(tab="Simplifications"));

  parameter SI.ThermalResistance RLong = 1/(2*1/RcondIC + 1/RcondPEN) "Cell equivalent longitudinal resistance without fuel and air channel ribs longitudinal conductivity" annotation(Dialog(tab="Simplifications"));

  // Equivalent thermal conductivities for the simplified block
  parameter SI.ThermalConductivity kCellLong = lX/(RLong*lZCell*lY) "Cell equivalent longitudinal conductivity" annotation(Dialog(tab="Simplifications"));
  parameter SI.ThermalConductivity kCellTrans = lZCell/(RTrans*lX*lY) "Cell equivalent transversal conductivity" annotation(Dialog(tab="Simplifications"));

  // Simplification Parameters
  parameter SI.Length lZBlock = nCellVertMult*(2*lZsolid + lZfc + lZac + lZpen) "Total height of simplified block" annotation(Dialog(tab="Simplifications"));
  parameter SI.Length lZCell = (2*lZsolid + lZfc + lZac + lZpen) "Total height of a cell" annotation(Dialog(tab="Simplifications"));

  parameter SI.Density rhoBlock = (rhoPEN*lZpen + 2*rhoIC*lZsolid)/lZCell "Total density of a cell without ribs in AC/FC" annotation(Dialog(tab="Simplifications"));
  parameter SI.SpecificHeatCapacity cpBlock = (cpPEN*rhoPEN*lZpen + 2*cpIC*rhoIC*lZsolid)/(rhoPEN*lZpen + 2*rhoIC*lZsolid) "Specific heat capacity of a cell without ribs in AC/FC" annotation(Dialog(tab="Simplifications"));

  SI.EnthalpyFlowRate dH[nX, nY];
  SI.Power PEl[nX, nY];
  SI.Voltage UCalc[nX, nY];

  Electrochem.Interfaces.SimpFlowInfoPort simpFlowInfoPortOut[nX,nY]
    annotation (Placement(transformation(extent={{64,20},{90,46}}),
        iconTransformation(extent={{64,20},{90,46}})));
  Electrochem.Interfaces.ReactionHeatInfoPort reactionHeatInfoPort[nX,nY]
    annotation (Placement(transformation(extent={{-6,20},{20,46}}),
        iconTransformation(extent={{-6,20},{20,46}})));
  Electrochem.Interfaces.SimpFlowInfoPort simpFlowInfoPortIn[nX,nY] annotation
    (Placement(transformation(extent={{-76,20},{-50,46}}), iconTransformation(
          extent={{-76,20},{-50,46}})));

  OpenTEMPEST.Heat.HeatSource2DNonUniformFV heatSource2DNonUniformFV(nX=nX, nY=
        nY) annotation (Placement(transformation(extent={{-72,-52},{-50,-30}})));
  ThermoPower.Thermal.DHTVolumes dHT_x0[nCellVertMult](each N=nY) annotation (
      Placement(transformation(extent={{-112,-10},{-92,10}}),
        iconTransformation(extent={{-112,-10},{-92,10}})));
  ThermoPower.Thermal.DHTVolumes dHT_xN[nCellVertMult](each N=nY) annotation (
      Placement(transformation(extent={{92,-10},{112,10}}), iconTransformation(
          extent={{92,-10},{112,10}})));
  ThermoPower.Thermal.DHTVolumes dHT_y0(N=nX) annotation (Placement(transformation(extent={{30,-106},
            {50,-86}}), iconTransformation(extent={{30,-106},{50,-86}})));
  ThermoPower.Thermal.DHTVolumes dHT_yN(N=nX) annotation (Placement(transformation(extent={{30,88},
            {50,108}}), iconTransformation(extent={{30,88},{50,108}})));
  Heat.DHTVolumes2D dHT2_z1(i=nX, j=nY) annotation (Placement(transformation(
          extent={{-80,74},{-60,94}}), iconTransformation(extent={{-80,74},{-60,
            94}})));
  Heat.DHTVolumes2D dHT2_z0(i=nX, j=nY) annotation (Placement(transformation(
          extent={{-80,-86},{-60,-66}}), iconTransformation(extent={{-80,-86},{
            -60,-66}})));
  OpenTEMPEST.Heat.ThermalCollector1D thermalCollector1D_xN(m=nCellVertMult, N=
        nY) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={40,0})));
  OpenTEMPEST.Heat.ThermalCollector1D thermalCollector1D_x0(m=nCellVertMult, N=
        nY) annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=-90,
        origin={-40,0})));
  OpenTEMPEST.Heat.ThermalConductor1D tC1_x0[nCellVertMult](each N=nY, each G=
        99^3)
    annotation (Placement(transformation(extent={{-82,-10},{-62,10}})));
  OpenTEMPEST.Heat.ThermalConductor1D tC1_xN[nCellVertMult](each N=nY, each G=
        99^3) annotation (Placement(transformation(extent={{64,-10},{84,10}})));

  OpenTEMPEST.Heat.Solid2D cellBlock(
    redeclare package SolidMat = TEMPEST.Solid.Material.Custom,
    nX=nX,
    nY=nY,
    Tstartbar=TStart,
    lX=lX,
    lY=lY,
    lZ=lZBlock,
    kCustom_long=kCellLong,
    kCustom_trans=kCellTrans,
    rhoCustom=rhoBlock,
    cpCustom=cpBlock)
    annotation (Placement(transformation(extent={{-20,-18},{18,18}})));

equation

  // Pass information to simplified block
  PEl = -UCalc .* reactionHeatInfoPort.I;
  UCalc[:,:] = (reactionHeatInfoPort[:,:].Vid .- (reactionHeatInfoPort[:,:].Vid .- reactionHeatInfoPort[:,:].Vop).*((0.27266 .+ 356130*exp(-0.01264*cellBlock.T[:,:]))./(0.27266 .+ 356130*exp(-0.01264*reactionHeatInfoPort[:,:].TPEN))));

  dH = PEl      .+ simpFlowInfoPortIn.HfAirRef .+ (simpFlowInfoPortIn.TAir .- simpFlowInfoPortIn.TAirRef) .* simpFlowInfoPortIn.cpAir .* simpFlowInfoPortIn.mfAir
                .+ simpFlowInfoPortOut.HfAirRef .+ (simpFlowInfoPortOut.TAir .- simpFlowInfoPortOut.TAirRef) .* simpFlowInfoPortOut.cpAir .* simpFlowInfoPortOut.mfAir
                .+ simpFlowInfoPortIn.HfFuelRef .+ (simpFlowInfoPortIn.TFuel .- simpFlowInfoPortIn.TFuelRef) .* simpFlowInfoPortIn.cpFuel .* simpFlowInfoPortIn.mfFuel
                .+ simpFlowInfoPortOut.HfFuelRef .+ (simpFlowInfoPortOut.TFuel .- simpFlowInfoPortOut.TFuelRef) .* simpFlowInfoPortOut.cpFuel .* simpFlowInfoPortOut.mfFuel;

  // Heat information passed to the simplified block
  heatSource2DNonUniformFV.power = nCellVertMult * dH;

  // x-direction

  // y-direction
  connect(dHT_y0, cellBlock.dhT_y0) annotation (Line(points={{40,-96},{40,-78},
          {12.3,-78},{12.3,-12.6}},color={255,127,0}));
  connect(dHT_yN, cellBlock.dhT_yN) annotation (Line(points={{40,98},{34,98},{34,
          54},{12.3,54},{12.3,12.6}}, color={255,127,0}));

  // z-direction
  connect(cellBlock.dhT2_z0,dHT2_z0)  annotation (Line(points={{-14.3,-12.6},{-26,
          -12.6},{-26,-76},{-70,-76}},
                           color={0,0,0}));
  connect(cellBlock.dhT2_z1,dHT2_z1)  annotation (Line(points={{-14.3,12.6},{-26,
          12.6},{-26,84},{-70,84}},
                         color={0,0,0}));

  // Internal heat connection
  connect(heatSource2DNonUniformFV.wall, cellBlock.dhT2_int) annotation (Line(
        points={{-61,-44.3},{-61,-48},{-1,-48},{-1,0}}, color={0,0,0}));
  connect(thermalCollector1D_x0.dHT1, cellBlock.dhT_x0) annotation (Line(points
        ={{-30,-1.88738e-15},{-25.95,-1.88738e-15},{-25.95,0},{-21.9,0}}, color
        ={255,127,0}));
  connect(cellBlock.dhT_xN, thermalCollector1D_xN.dHT1) annotation (Line(points
        ={{19.9,0},{24.95,0},{24.95,1.88738e-15},{30,1.88738e-15}}, color={255,127,
          0}));
  connect(thermalCollector1D_xN.dHT0, tC1_xN.dHT0) annotation (Line(points={{50,
          -1.88738e-15},{57.2,-1.88738e-15},{57.2,0},{64.4,0}}, color={255,127,0}));
  connect(tC1_xN.dHT1, dHT_xN)  annotation (Line(points={{83.6,0},{102,0}}, color={255,127,0}));
  connect(tC1_x0.dHT1, thermalCollector1D_x0.dHT0) annotation (Line(points={{-62.4,
          0},{-56.2,0},{-56.2,1.88738e-15},{-50,1.88738e-15}}, color={255,127,0}));
  connect(dHT_x0, tC1_x0.dHT0)
    annotation (Line(points={{-102,0},{-81.6,0}}, color={255,127,0}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={28,108,200},
          fillColor={50,173,200},
          fillPattern=FillPattern.CrossDiag,
          lineThickness=0.5)}),                                  Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    experiment(__Dymola_Algorithm="Cvode"));
end SimplifiedCellBlockCF;

within OpenTEMPEST.BOP.Reactors.MeOH;
model MeOHEnergyBalance
  "Methanol reactor model, balancing component and heat streams, based on defined inflow and fixed conversion rate"
  import SI = Modelica.SIunits;
  replaceable package Medium = OpenTEMPEST.Medium.Fuel_MethanolReactor;

  extends Modelica.Icons.UnderConstruction;

  parameter Real conversion = 0.96 "ratio C/H>=0.25; const conversion rate nMeOH/(nCO+nCO2)in";
  parameter Real maxH2conversion= 0.99 "value between 0...1; maximum H2 conversion";
  parameter SI.Temperature TReactor=250 + 273.15 "[K] temperature of Purge gas (unequal reactor temperature)";
  parameter SI.Temperature TPurge=250 + 273.15 "[K] temperature of Purge gas (unequal reactor temperature)";
  parameter SI.Temperature TProduct=40 + 273.15 "[K] temperature of condensated product flow";
  parameter SI.Pressure pReactor=1 "pressure gas Purge after flush/sparator (unequal reator pressure)";
  parameter Real fac = 1 "factor for accoutning for additional heat sinks";
  parameter SI.SpecificEnthalpy dhvMeOH=1165000 " [J/kg] latent heat of evaporization Methanol at respective T,p";
            //https://www.engineeringtoolbox.com/methanol-methyl-alcohol-properties-CH3OH-d_2031.html
//   parameter SI.SpecificEnthalpy dhvWater=2256000 " [J/kg] latent heat of evaporization H2O at respective T,p";
            //https://www.engineeringtoolbox.com/methanol-methyl-alcohol-properties-CH3OH-d_2031.html
  SI.MolarFlowRate nfMeOH " molar flow rate of methanol";
  SI.SpecificEnthalpy hPurge " specific enthalpy of Purge stream";
  SI.SpecificEnthalpy hProduct " specific enthalpy of Product stream";
  Real dGwgs;
  Real Kwgs;
  SI.MassFlowRate mfPurge "mass flow rate of outgoing gas";
  SI.MassFlowRate mfProduct "mass flow rate of outgoing product";
  SI.MolarFlowRate nfDry[Medium.nXi] "molar flow rates of species in dry inflow";
  SI.MolarFlowRate nfPurge[Medium.nXi] "molar flow rates of species in Purge";
  Real XPurge[Medium.nXi] "mass fractions in gas outlet";
  Real XIn[Medium.nXi] "mass fractions inflow"; // (start={0.117773596,0.05271881,0.261454681,0.55718662,0.003588371,0.007277923,1e-6}) ;
  Real YPurge[Medium.nXi]; //(start = {0.30143255,7.99E-05,0.003377404,0.09457149,0.60049284,4.58E-05,0}) "molar fractions in Purge";
  Real YDry[Medium.nXi]; //(start = {0.69590056,0,0.07502401,0.22907542,0,0,0})  "molar fractions of dry inflow (water removed)";
  //   SI.MolarFlowRate nfProduct[Medium.nXi];
  //   Real YProduct[Medium.nXi] "molar fractions in Product stream";
  SI.MolarMass MMDry "mean molar mass of inflow";
  SI.MolarMass MMPurge "mean molar mass of Purge gas";
  SI.MolarMass MMProduct " mean molar mass of Product fluid";
  SI.MolarFlowRate nfPurgeTot "total molar flow of Purge gas";
  SI.MolarFlowRate nfDryTot "total molar flow of inflow";
  Real ratio "ratio that relates H2 consumption with CO and CO2 fractions at inlet";

  ThermoPower.Gas.FlangeA infl(redeclare package Medium =
        OpenTEMPEST.Medium.Fuel_MethanolReactor) annotation (Placement(
        transformation(extent={{-128,-40},{-100,-12}}), iconTransformation(
          extent={{-128,-40},{-100,-12}})));
  ThermoPower.Gas.FlangeB outflPurge(redeclare package Medium =
        OpenTEMPEST.Medium.Fuel_MethanolReactor) annotation (Placement(
        transformation(extent={{100,0},{128,28}}), iconTransformation(extent={{
            100,0},{128,28}})));

  Modelica.Blocks.Interfaces.RealOutput outMfProduct annotation (Placement(transformation(
        origin={100,0},
        extent={{-16,-16},{16,16}},
        rotation=0), iconTransformation(extent={{102,-38},{126,-14}})));
  Modelica.Blocks.Sources.RealExpression realWProduct(y=-mfProduct) annotation (Placement(transformation(extent={{54,-8},{66,8}})));

  ThermoPower.Thermal.HT hT annotation (Placement(transformation(extent={{-16,60},
             {16,92}}), iconTransformation(extent={{-12,34},{12,58}})));

  ThermoPower.Gas.SourceMassFlow sourceMFPurge(
    redeclare package Medium = OpenTEMPEST.Medium.Fuel_MethanolReactor,
    use_in_w0=true,
    use_in_T=true,
    use_in_X=true)
    annotation (Placement(transformation(extent={{40,30},{60,50}})));
  Modelica.Blocks.Sources.RealExpression realWPurge(y=mfPurge)
                                                    annotation (Placement(transformation(extent={{24,40},{38,54}})));
  Modelica.Blocks.Sources.RealExpression realTPurge(y=TPurge)
                                                    annotation (Placement(transformation(extent={{30,50},{44,66}})));
//   RealVector realXPurge(n=7, y=XPurge) annotation (Placement(transformation(extent={{36,62},{50,78}})));

  inner ThermoPower.System system annotation (Placement(transformation(extent={{-134,56},{-114,76}})));
  ThermoPower.Gas.SinkPressure sinkPressure(redeclare package Medium =
        OpenTEMPEST.Medium.Fuel_MethanolReactor, use_in_p0=true)
    annotation (Placement(transformation(extent={{-72,-10},{-52,10}})));
  ThermoPower.Gas.SensP sensP(redeclare package Medium =
        OpenTEMPEST.Medium.Fuel_MethanolReactor)
    annotation (Placement(transformation(extent={{70,34},{90,54}})));

protected
  SI.SpecificEnthalpy hH2=Modelica.Media.IdealGases.Common.Functions.h_T(
      T=TPurge,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.H2);
  SI.SpecificEnthalpy hCH4=Modelica.Media.IdealGases.Common.Functions.h_T(
      T=TPurge,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.CH4);
  SI.SpecificEnthalpy hCO2=Modelica.Media.IdealGases.Common.Functions.h_T(
      T=TPurge,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.CO2);
  SI.SpecificEnthalpy hCO=Modelica.Media.IdealGases.Common.Functions.h_T(
      T=TPurge,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.CO);
  SI.SpecificEnthalpy hH2O=Modelica.Media.IdealGases.Common.Functions.h_T(
      T=TPurge,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.H2O);
  SI.SpecificEnthalpy hN2=Modelica.Media.IdealGases.Common.Functions.h_T(
      T=TPurge,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.N2);
  SI.SpecificEnthalpy hMeOH=Modelica.Media.IdealGases.Common.Functions.h_T(
      T=TPurge,
      exclEnthForm=false,
      refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
      data=Modelica.Media.IdealGases.Common.SingleGasesData.CH3OH);

equation
  XPurge =Modelica.Media.Interfaces.PartialMixtureMedium.moleToMassFractions(
    YPurge, OpenTEMPEST.Medium.Fuel_MethanolReactor.MMX);
  //   XProduct = Modelica.Media.Interfaces.PartialMixtureMedium.moleToMassFractions(YProduct, TEMPEST.Medium.Fuel_MethanolReactor.MMX);
  MMPurge = sum(YPurge[:]* Medium.MMX[:]);
  MMDry    = sum(YDry[:]   * Medium.MMX[:]);
  MMProduct = Medium.MMX[7]; //sum(YProduct[:] * Medium.MMX[:]);
  mfPurge   = nfPurgeTot*MMPurge;
  mfProduct = nfMeOH*MMProduct;
  nfPurgeTot = sum(nfPurge[:]);
  nfDryTot   = sum(nfDry[1:4]) + nfDry[6]+nfDry[7];
  for i in 1:Medium.nXi loop
    nfDry[i] =XIn[i]*infl.m_flow/OpenTEMPEST.Medium.Fuel_MethanolReactor.MMX[i];
    YDry[i] = nfDry[i]/nfDryTot;
    YPurge[i] = nfPurge[i]/nfPurgeTot;
  end for;

// MeOH production rate:
  ratio = nfDry[3]/(nfDry[3]+nfDry[4])*3  + nfDry[4]/(nfDry[3]+nfDry[4])*2;
  if ratio*conversion*(nfDry[3]+nfDry[4])<= nfDry[1]*maxH2conversion then
    nfMeOH = conversion*(nfDry[3]+nfDry[4]);
  else
    nfMeOH = nfDry[1]/ratio *maxH2conversion;
  end if;

// Species balance (7 eq's for 7 species):
  nfPurge[1] + nfPurge[5] +2*nfMeOH + 2*nfPurge[7] = nfDry[1] + nfDry[5] + 2*nfDry[7];  //H2 balance
  nfPurge[4] + nfPurge[3] + nfMeOH + nfPurge[7] = nfDry[3] + nfDry[4] + nfDry[7];   // C balance
  nfPurge[5] + 2*nfPurge[3] + nfPurge[4] + nfMeOH + nfPurge[7] = 2*nfDry[3] + nfDry[4] + nfDry[5]+ nfDry[7]; //O balance
  nfPurge[2] = nfDry[2];
  nfPurge[6] = nfDry[6];
  nfPurge[7] = 0;
  Kwgs = pReactor*YPurge[4]*pReactor*YPurge[5] / (pReactor*YPurge[1]*pReactor*YPurge[3]);

  dGwgs= 32.1153*(TReactor) - 3.5211E4; // Marius' shortcut from NASA
  Kwgs = exp(-dGwgs/Modelica.Constants.R/TReactor);

// Energy balance:
  -hT.Q_flow= (infl.m_flow*inStream(infl.h_outflow) - mfPurge*hPurge - mfProduct*(hProduct + dhvMeOH))*fac; //- mfWater*(hWater + dhvWater)
  hProduct = -7459000 + 81.2*(TProduct-298.15); // H0 + cp*dT
  hPurge   = hH2*XPurge[1] + hCH4*XPurge[2] + hCO2*XPurge[3] + hCO*XPurge[4] + hH2O*XPurge[5] + hN2*XPurge[6] + hMeOH*XPurge[7];

// Boundary condition:
  XIn = inStream(infl.Xi_outflow);

//connections
  for i in 1:Medium.nXi loop
    sourceMFPurge.in_X[i]=XPurge[i];
  end for;

  connect(realWPurge.y, sourceMFPurge.in_w0) annotation (Line(points={{38.7,47},{44,47},{44,45}}, color={0,0,127}));
  connect(realTPurge.y, sourceMFPurge.in_T) annotation (Line(points={{44.7,58},{50,58},{50,45}}, color={0,0,127}));
//   connect(realXPurge.y, sourceMFPurge.in_X) annotation (Line(points={{50.7,70},{56,70},{56,45}}, color={0,0,127}));
  connect(realWProduct.y, outMfProduct) annotation (Line(points={{66.6,0},{100,0}}, color={0,0,127}));
  connect(sourceMFPurge.flange, outflPurge) annotation (Line(points={{60,40},
          {88,40},{88,14},{114,14}},                                                    color={159,159,223}));

  connect(infl, sinkPressure.flange) annotation (Line(points={{-114,-26},{-92,
          -26},{-92,0},{-72,0}},                           color={159,159,
          223}));
  connect(sourceMFPurge.flange, sensP.flange)
    annotation (Line(points={{60,40},{80,40}}, color={159,159,223}));
  connect(sensP.p, sinkPressure.in_p0) annotation (Line(points={{87,50},{90,
          50},{90,96},{-68.45,96},{-68.45,5.95}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,80}}),                                   graphics={
          Rectangle(
          extent={{-100,34},{100,-86}},
          lineColor={28,108,200},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),              Text(
          extent={{-64,6},{62,-54}},
          lineColor={255,255,255},
          textString="MeOH process")}),                          Diagram(
        coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
            80}}),                                   graphics={Text(
          extent={{-128,76},{38,-20}},
          lineColor={238,46,47},
          textString="conversion rate should be dependent on inflow not fixed")}),
    experiment(StopTime=10000, __Dymola_Algorithm="Cvode"),
    __Dymola_experimentSetupOutput(equidistant=false),
    __Dymola_experimentFlags(
      Advanced(GenerateVariableDependencies=false, OutputModelicaCode=true),
      Evaluate=false,
      OutputCPUtime=true,
      OutputFlatModelica=false));
end MeOHEnergyBalance;

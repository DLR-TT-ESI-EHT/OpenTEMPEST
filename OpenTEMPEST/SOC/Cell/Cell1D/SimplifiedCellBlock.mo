within OpenTEMPEST.SOC.Cell.Cell1D;
model SimplifiedCellBlock
  import SI = Modelica.SIunits;
  parameter Integer N(min=3)=5 "number of CVs in each layer";
  //parameter TEMPEST.ECReactorModels.Cell.baseModelsFV_SOC.solidMaterialOptions matOpt=TEMPEST.ECReactorModels.Cell.baseModelsFV_SOC.solidMaterialOptions.Crofer22APU "options for interconnect material";
  parameter Integer nCellVertMult=1 "number of cells this object is representing in z direction";

  parameter SI.Temperature TStart "Starting temperature";
  parameter SI.Length lX=11.2e-2 "Length of active part of cell";
  parameter SI.Length lY=7.2e-2 "Width of active part of cell";
  parameter SI.Length lZ=0.5e-3 "Thickness of IC over the active part of cell";
  parameter SI.ThermalConductivity k_solidCustomLong=0.2812 "Thermal Conductivity of solid along x-direction" annotation (Dialog(tab="Interconnect"));
  parameter SI.ThermalConductivity k_solidCustomTrans=0.2812 "Thermal Conductivity of solid along z-direction" annotation (Dialog(tab="Interconnect"));
  parameter SI.SpecificHeatCapacity cp_solidCustom=463.8 "Specific heat capacity of solid parallel to windows" annotation (Dialog(tab="Interconnect"));
  parameter SI.Density rho_solidCustom(displayUnit="kg/m3") = 1330 "Density of solid parallel to windows" annotation (Dialog(tab="Interconnect"));

  parameter SI.Density rho_custom = 1 "Density of custom material" annotation (Dialog(enable=matOpt ==
          TEMPEST.ECReactorModels.Cell.BaseModelsFV.solidMaterialOptions.Other));
  parameter SI.SpecificHeatCapacity cp_custom = 1 "Specific Heat capacity of custom material" annotation (Dialog(enable=matOpt ==
          TEMPEST.ECReactorModels.Cell.BaseModelsFV.solidMaterialOptions.Other));
  parameter SI.ThermalConductivity k_custom = 1 "Thermal condusctivity of custom material" annotation (Dialog(enable=matOpt ==
          TEMPEST.ECReactorModels.Cell.BaseModelsFV.solidMaterialOptions.Other));

  OpenTEMPEST.Heat.Solid1D solid1D(
    redeclare package SolidMat = OpenTEMPEST.Solid.Material.Custom,
    N=N,
    Tstartbar=TStart,
    lX=lX,
    lY=lY,
    lZ=lZ*nCellVertMult,
    kCustom_long=k_solidCustomLong,
    kCustom_trans=k_solidCustomTrans,
    rhoCustom=rho_solidCustom,
    cpCustom=cp_solidCustom)
    annotation (Placement(transformation(extent={{-16,-14},{16,16}})));
  ThermoPower.Thermal.HT hTactive_x0[nCellVertMult]
    annotation (Placement(transformation(extent={{-116,-10},{-96,10}}),
        iconTransformation(extent={{-110,-4},{-96,10}})));
  ThermoPower.Thermal.HT hTactive_x1[nCellVertMult]
    annotation (Placement(transformation(extent={{88,-10},{108,10}}),
        iconTransformation(extent={{94,-4},{108,10}})));
  ThermoPower.Thermal.DHTVolumes dHTactive_y0(each N=N) annotation (Placement(
        transformation(extent={{-42,86},{46,104}}),   iconTransformation(extent={{-40,92},
            {46,104}})));
  ThermoPower.Thermal.DHTVolumes dHTactive_y1(each N=N) annotation (Placement(
        transformation(extent={{-40,-108},{48,-90}}), iconTransformation(extent={{-38,
            -102},{48,-90}})));
  ThermoPower.Thermal.DHTVolumes dHTactive_z0(each N=N) annotation (Placement(
        transformation(extent={{42,-70},{70,-44}}),   iconTransformation(extent={{44,-70},
            {68,-50}})));
  ThermoPower.Thermal.DHTVolumes dHTactive_z1(each N=N) annotation (Placement(
        transformation(extent={{-70,44},{-42,70}}),   iconTransformation(extent={{-66,52},
            {-42,70}})));
   ThermoPower.Thermal.HeatSource1DNonUniformFV heatSource1DNonUniformFV(Nw=N)
           annotation (Placement(transformation(extent={{-86,-42},{-66,-22}})));
  Electrochem.Interfaces.SimpFlowInfoPort simpFlowInfoPortIn[N] annotation (
      Placement(transformation(extent={{-62,-82},{-42,-62}}),
        iconTransformation(extent={{-62,-82},{-42,-62}})));
  Electrochem.Interfaces.ReactionHeatInfoPort reactionHeatInfoPort[N]
    annotation (Placement(transformation(extent={{12,-80},{32,-60}}),
        iconTransformation(extent={{12,-80},{32,-60}})));
  Electrochem.Interfaces.SimpFlowInfoPort simpFlowInfoPortOut[N] annotation (
      Placement(transformation(extent={{86,-80},{106,-60}}), iconTransformation(
          extent={{86,-80},{106,-60}})));

  SI.EnthalpyFlowRate dH[N];
  SI.Power PEl[N];
  SI.Voltage UCalc[N];
  Modelica.Thermal.HeatTransfer.Components.ThermalCollector thermalCollector(m=nCellVertMult) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={52,0})));
  Modelica.Thermal.HeatTransfer.Components.ThermalCollector thermalCollector1(m=nCellVertMult) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-58,0})));

  Modelica.Thermal.HeatTransfer.Components.ThermalConductor tcX1[nCellVertMult](each G=99^3)
                                                                                          annotation (Placement(transformation(extent={{-96,-10},{-76,10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor tcX2[nCellVertMult](each G=99^3) annotation (Placement(transformation(extent={{68,-10},{88,10}})));

equation

  for k in 1:N loop
    PEl[k] = -UCalc[k]*reactionHeatInfoPort[k].I;
    UCalc[k] = (reactionHeatInfoPort[k].Vid - (reactionHeatInfoPort[k].Vid - reactionHeatInfoPort[k].Vop)*((0.27266 + 356130*exp(-0.01264*solid1D.T[k]))/(0.27266 + 356130*exp(-0.01264*reactionHeatInfoPort[k].TPEN))));//*(140*realPassLeakage.y + 1));
    dH[k] = PEl[k]
        + simpFlowInfoPortIn[k].HfAirRef + (simpFlowInfoPortIn[k].TAir-simpFlowInfoPortIn[k].TAirRef) * simpFlowInfoPortIn[k].cpAir * simpFlowInfoPortIn[k].mfAir
        + simpFlowInfoPortOut[k].HfAirRef + (simpFlowInfoPortOut[k].TAir-simpFlowInfoPortOut[k].TAirRef) * simpFlowInfoPortOut[k].cpAir * simpFlowInfoPortOut[k].mfAir
        + simpFlowInfoPortIn[k].HfFuelRef + (simpFlowInfoPortIn[k].TFuel-simpFlowInfoPortIn[k].TFuelRef) * simpFlowInfoPortIn[k].cpFuel * simpFlowInfoPortIn[k].mfFuel
        + simpFlowInfoPortOut[k].HfFuelRef + (simpFlowInfoPortOut[k].TFuel-simpFlowInfoPortOut[k].TFuelRef) * simpFlowInfoPortOut[k].cpFuel * simpFlowInfoPortOut[k].mfFuel;
  end for;

  heatSource1DNonUniformFV.power=nCellVertMult*dH[:];
  connect(heatSource1DNonUniformFV.wall, solid1D.dhT_int) annotation (Line(points={{-76,-35},{-76,-42},{0,-42},{0,1.6}}, color={255,127,0}));

  connect(solid1D.dhT_z1, dHTactive_z1) annotation (Line(points={{-11.2,11.5},{-11.2,
          57},{-56,57}}, color={255,127,0}));

  connect(solid1D.dhT_y0, dHTactive_y0) annotation (Line(points={{11.2,-9.5},{11.2,
          -14},{38,-14},{38,64},{2,64},{2,95}}, color={255,127,0}));
  connect(solid1D.dhT_y1, dHTactive_y1) annotation (Line(points={{11.2,11.5},{
          11.2,20},{-36,20},{-36,-64},{4,-64},{4,-99}},
                                                   color={255,127,0}));
  connect(solid1D.dhT_z0, dHTactive_z0) annotation (Line(points={{-11.2,-9.5},{-11.2,
          -57},{56,-57}}, color={255,127,0}));

  connect(thermalCollector1.port_b, solid1D.hT_x0) annotation (Line(points={{-48,0},
          {-32,0},{-32,1},{-17.6,1}},                                                                               color={191,0,0}));
  connect(solid1D.hT_xN, thermalCollector.port_b) annotation (Line(points={{17.6,1},{33.8,1},{33.8,0},{42,0}}, color={191,0,0}));

  connect(hTactive_x0, tcX1.port_a) annotation (Line(points={{-106,0},{-96,0}}, color={191,0,0}));
  connect(thermalCollector1.port_a, tcX1.port_b) annotation (Line(points={{-68,0},{-76,0}}, color={191,0,0}));
  connect(hTactive_x1, tcX2.port_b) annotation (Line(points={{98,0},{88,0}}, color={191,0,0}));
  connect(tcX2.port_a, thermalCollector.port_a) annotation (Line(points={{68,0},{62,0}}, color={191,0,0}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={28,108,200},
          fillColor={28,108,200},
          fillPattern=FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(revisions="<html>
</html>", info="<html>
<h2>Simplified 1D Cell Block</h2>

<p>
This model represents a simplified 1D block of cells and associated thermal volumes. It discretizes 
the block of cells into N control volumes along the x-direction and optionally represents 
multiple cells stacked in the z-direction.
</p>
</html>"));
end SimplifiedCellBlock;

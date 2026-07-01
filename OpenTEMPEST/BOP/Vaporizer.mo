within OpenTEMPEST.BOP;
model Vaporizer
  "Evaporator model. Input is mass flow in water. Output is water vapor as gas mix. All extensive variables are depending on the mass flow."
 import SI = Modelica.SIunits;
 parameter Boolean thermalInput = false;
 parameter SI.AbsolutePressure pStart = 101325;
  replaceable package MediumInput = Modelica.Media.Water.StandardWater constrainedby
    Modelica.Media.Interfaces.PartialTwoPhaseMedium
    annotation(choicesAllMatching = true);
  replaceable package MediumOutput = Medium.Fuel_CH4          constrainedby
    Modelica.Media.Interfaces.PartialMedium
   annotation(choicesAllMatching = true);

   parameter SI.MassFlowRate mfNom = 8/1000 "Massflowrate";

  ThermoPower.Water.FlangeA waterInlet(redeclare package Medium = MediumInput)
                                           annotation (Placement(transformation(
          rotation=0, extent={{-134,-7},{-112,15}}), iconTransformation(extent={
            {-100,-100},{-48,-48}})));
  ThermoPower.Gas.FlangeB gasOutlet1(redeclare package Medium = MediumOutput)
                           annotation (Placement(transformation(rotation=0,
          extent={{106,-4},{126,16}}), iconTransformation(extent={{50,50},{100,100}})));
protected
  parameter SI.MassFlowRate mfDesign = 8/1000;
  parameter Real powerDesign = mfDesign*(4.18*(100-20) + 2.26e3)*1e3*1.3;
public
  ThermoPower.Water.SinkPressure sinkPressure(
    redeclare package Medium = MediumInput,
    p0=100000,
    use_in_p0=true)
    annotation (Placement(transformation(extent={{10,-2},{30,18}})));
  ThermoPower.Water.DrumEquilibrium drumEquilibrium1(
    redeclare package Medium = MediumInput,
    Vd=mfNom/mfDesign*0.01,
    Mm=mfNom/mfDesign*1,
    cm=450,
    allowFlowReversal=false,
    pstart=pStart,
    Vlstart=0.005,
    noInitialPressure=true)
    annotation (Placement(transformation(extent={{-68,-8},{-48,12}})));
  ThermoPower.Thermal.HeatSource1DFV heatSource1DFV(Nw=5) annotation (Placement(
        transformation(
        extent={{-10,10},{10,-10}},
        rotation=90,
        origin={-2,-42})));
  ThermoPower.Thermal.HT_DHTVolumes hT_DHTVolumes(N=5)
    annotation (Placement(transformation(extent={{-56,-46},{-42,-32}})));
  ThermoPower.Gas.SensW sensW(redeclare package Medium = MediumInput)
    annotation (Placement(transformation(extent={{-36,14},{-16,-6}})));
  ThermoPower.Gas.SourceMassFlow sourceMassFlow(
    redeclare package Medium = MediumOutput,
    p0=100000,
    Xnom={1e-5,1e-5,1e-5,1e-5,1,1e-5},
    use_in_w0=true,
    use_in_T=true)
    annotation (Placement(transformation(extent={{42,18},{62,-2}})));
  ThermoPower.Water.SensT sensT(redeclare package Medium = MediumInput)
    annotation (Placement(transformation(extent={{-18,14},{2,-6}})));
  ThermoPower.Thermal.MetalWallFV metalWallFV(
    Nw=5,
    M=mfNom/mfDesign*10,
    cm=450,
    WallRes=true,
    UA_ext=mfNom/mfDesign*1500,
    Tstartbar=1.1*(45.668*ln(pStart) - 426.44) + 273.15)
                 annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-26,-42})));
  ThermoPower.Water.SensW sensW1(redeclare package Medium = MediumInput)
    annotation (Placement(transformation(extent={{-98,-12},{-78,8}})));
  ThermoPower.Gas.SensP sensP(redeclare package Medium = MediumOutput)
    annotation (Placement(transformation(extent={{104,2},{84,22}})));
  Modelica.Blocks.Continuous.LimPID PIDPessure(
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=0.198712,
    Ti=3.67,
    yMax=mfNom/mfDesign*powerDesign,
    yMin=0) annotation (Placement(transformation(extent={{48,-30},{28,-50}})));
  ThermoPower.Electrical.Load load(Pnom(displayUnit="kW") = 5000, usePowerInput=
       true) annotation (Placement(transformation(extent={{10,-94},{30,-74}})));
  ThermoPower.Electrical.PowerConnection port annotation (Placement(
        transformation(rotation=0, extent={{20,-120},{40,-100}})));
  ThermoPower.Thermal.DHTVolumes dHTVolumes(N=5)
                                            if thermalInput
    annotation (Placement(transformation(extent={{-34,-120},{-14,-100}})));
  ThermoPower.Water.ValveLin valveLin(redeclare package Medium =
        Modelica.Media.Water.WaterIF97_pT, Kv=mfNom/200000)
    annotation (Placement(transformation(extent={{-46,14},{-34,26}})));
  Modelica.Blocks.Continuous.LimPID PIDValue(
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=100,
    Ti=3.3,
    yMax=1,
    yMin=0.1/10000)
    annotation (Placement(transformation(extent={{-62,38},{-42,58}})));
  ThermoPower.Water.SensP sensP1(redeclare package Medium =
        Modelica.Media.Water.WaterIF97_pT)
    annotation (Placement(transformation(extent={{-46,-4},{-26,-24}})));
  Modelica.Blocks.Sources.Ramp Valve2(
    height=0,
    duration=0,
    offset=500000,
    startTime=0)
    annotation (Placement(transformation(extent={{92,-50},{72,-30}})));
  ThermoPower.Water.SensP sensP2(redeclare package Medium =
        Modelica.Media.Water.WaterIF97_pT)
    annotation (Placement(transformation(extent={{-8,12},{12,32}})));
equation
  connect(hT_DHTVolumes.HT_port,drumEquilibrium1. wall)
    annotation (Line(points={{-57.4,-39},{-58,-39},{-58,-7}},color={191,0,0}));
  connect(sensW.w,sourceMassFlow. in_w0) annotation (Line(points={{-19,-2},{-16,
          -2},{-16,-12},{46,-12},{46,3}},color={0,0,127}));
  connect(sensW.outlet,sensT. inlet) annotation (Line(points={{-20,8},{-14,8}},
                           color={159,159,223}));
  connect(sensT.outlet,sinkPressure. flange)
    annotation (Line(points={{-2,8},{10,8}},               color={0,0,255}));
  connect(sensT.T,sourceMassFlow. in_T) annotation (Line(points={{0,-2},{10,-2},
          {10,-10},{52,-10},{52,3}}, color={0,0,127}));
  connect(hT_DHTVolumes.DHT_port,metalWallFV. int) annotation (Line(points={{
          -41.3,-39},{-36,-39},{-36,-42},{-29,-42}}, color={255,127,0}));
  connect(sensW1.outlet,drumEquilibrium1. feed) annotation (Line(points={{-82,-6},
          {-76,-6},{-76,-2.4},{-67,-2.4}},
                                        color={0,0,255}));
  connect(sensP.p, sinkPressure.in_p0) annotation (Line(points={{87,18},{86,18},
          {86,32},{16,32},{16,16.4}}, color={0,0,127}));
  connect(sensW1.inlet, waterInlet) annotation (Line(points={{-94,-6},{-104,-6},
          {-104,4},{-123,4}},
                           color={0,0,255}));
  connect(sensP.flange, gasOutlet1) annotation (Line(points={{94,8},{106,8},{106,
          6},{116,6}}, color={159,159,223}));
  connect(heatSource1DFV.wall, metalWallFV.ext)
    annotation (Line(points={{-5,-42},{-22.9,-42}}, color={255,127,0}));
  connect(port,load. port) annotation (Line(points={{30,-110},{30,-75.4},{20,
          -75.4}}, color={0,0,255}));
  connect(load.referencePower, heatSource1DFV.power) annotation (Line(points={{16.7,
          -84},{16,-84},{16,-42},{2,-42}},   color={0,0,127}));
  connect(metalWallFV.ext, dHTVolumes) annotation (Line(points={{-22.9,-42},{
          -18,-42},{-18,-110},{-24,-110}},
                                       color={255,127,0}));
  connect(drumEquilibrium1.steam, valveLin.inlet) annotation (Line(points={{-51.2,
          9.2},{-45.6,9.2},{-45.6,20},{-46,20}},       color={0,0,255}));
  connect(valveLin.outlet, sensW.inlet) annotation (Line(points={{-34,20},{-34,
          14},{-32,14},{-32,8}}, color={0,0,255}));
  connect(sensW.w, PIDValue.u_m) annotation (Line(points={{-19,-2},{-18,-2},{
          -18,32},{-52,32},{-52,36}}, color={0,0,127}));
  connect(sensW1.w, PIDValue.u_s) annotation (Line(points={{-80,4},{-74,4},{-74,
          48},{-64,48}}, color={0,0,127}));
  connect(sensP1.flange, valveLin.inlet) annotation (Line(points={{-36,-10},{
          -46,-10},{-46,10},{-45.6,10},{-45.6,20},{-46,20}}, color={0,0,255}));
  connect(sensP1.p, PIDPessure.u_m)
    annotation (Line(points={{-28,-20},{38,-20},{38,-28}}, color={0,0,127}));
  connect(Valve2.y, PIDPessure.u_s)
    annotation (Line(points={{71,-40},{60,-40},{60,-40},{50,-40}},
                                                 color={0,0,127}));
  connect(PIDPessure.y, heatSource1DFV.power)
    annotation (Line(points={{27,-40},{14,-40},{14,-42},{2,-42}},
                                                color={0,0,127}));
  connect(PIDValue.y, valveLin.cmd) annotation (Line(points={{-41,48},{-38,48},
          {-38,24.8},{-40,24.8}}, color={0,0,127}));
  connect(port, port) annotation (Line(
      points={{30,-110},{30,-110}},
      color={0,0,255},
      thickness=0.5));
  connect(dHTVolumes, dHTVolumes)
    annotation (Line(points={{-24,-110},{-24,-110}}, color={255,127,0}));
  connect(sensT.outlet, sensP2.flange)
    annotation (Line(points={{-2,8},{0,8},{0,18},{2,18}}, color={0,0,255}));
  connect(sensP2.flange, sinkPressure.flange) annotation (Line(points={{2,18},{4,
          18},{4,8},{10,8}},                color={0,0,255}));
  connect(sourceMassFlow.flange, sensP.flange) annotation (Line(points={{62,8},
          {78,8},{78,8},{94,8}}, color={159,159,223}));
  annotation (Icon(graphics={
        Ellipse(
          extent={{-78,80},{82,-80}},
          lineColor={128,128,128},
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-58,0},{-58,-6},{-56,-16},{-50,-30},{-42,-42},{-36,-46},{-30,
              -50},{-20,-56},{-14,-58},{-6,-60},{-4,-60},{2,-60},{8,-60},{14,-58},
              {24,-56},{32,-52},{38,-48},{44,-42},{50,-36},{54,-28},{60,-18},{62,
              -8},{62,0},{-58,0}},
          lineColor={128,128,128},
          fillPattern=FillPattern.Solid,
          fillColor={0,0,255}),
        Polygon(
          points={{-58,0},{-56,16},{-48,34},{-34,48},{-24,54},{-14,58},{-4,60},{
              2,60},{12,60},{22,56},{32,52},{38,48},{48,40},{54,30},{58,22},{60,
              14},{62,6},{62,0},{-58,0}},
          lineColor={128,128,128},
          fillColor={170,255,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-34,-80},{42,-100}},
          lineColor={128,128,128},
          fillColor={255,128,0},
          fillPattern=FillPattern.Solid)}));
end Vaporizer;

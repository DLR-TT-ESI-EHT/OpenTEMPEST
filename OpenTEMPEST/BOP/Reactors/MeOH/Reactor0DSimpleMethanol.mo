within OpenTEMPEST.BOP.Reactors.MeOH;
model Reactor0DSimpleMethanol
  "Simple 0D reactor model with SMR and WGS reaction"
  extends Modelica.Icons.UnderConstruction;

  import SI = Modelica.SIunits;
  parameter SI.AbsolutePressure pStart "Reactor pressure for initialization";
  parameter SI.Temperature TStart "Reactor temperature for initialization";
  parameter SI.MassFraction xStart[Medium.nX] = Medium.reference_X "Mass fraction for intialization";
  parameter SI.Volume VReactor "total free reactor volume for gas flow";
  parameter Real VolRatio "ratio of free volume to the catalyst+substrat volume";
  parameter SI.SpecificHeatCapacity c "heat capacity of the catalyst+substrat";
  //parameter SI.Length tau "thickness of catalyst surface";
  parameter SI.Area A "Total catalyst surface area in contact with gas";
  parameter SI.CoefficientOfHeatTransfer alphaCat "heat transfer coeff from catalyst to gas";
  parameter SI.ThermalConductance GOuter "Average thermal conduction from catalyst to the outer heat port";
  parameter SI.Density rho "density of catalyst+substrat";

  replaceable package Medium = TEMPEST.Medium.Gas_MethanolReactor;
                                                              // #FB
  parameter ThermoPower.Choices.Init.Options initOpt = ThermoPower.Choices.Init.Options.noInit "initialization options";
  parameter Boolean noInitialPressure=false
    "Remove initial equation on pressure"
    annotation (choices(checkBox=true));
  ThermoPower.Gas.FlangeA infl(redeclare package Medium = Medium) "Inlet of reactor" annotation (
    Placement(transformation(extent={{-100,-20},{-60,20}}),      iconTransformation(extent={{-20,
            -100},{20,-60}})));
  ThermoPower.Gas.FlangeB outfl(redeclare package Medium = Medium) "Outlet of reactor" annotation (
    Placement(transformation(extent={{60,-20},{100,20}}),      iconTransformation(extent={{-18,60},
            {20,100}})));
  ThermoPower.Gas.Plenum plenum(redeclare package Medium = Medium, V=VReactor,
    pstart=pStart,
    Tstart=TStart,
    Xstart=xStart,
    noInitialPressure=noInitialPressure)
    annotation (Placement(transformation(extent={{10,-20},{50,20}})));
  Modelica.Thermal.HeatTransfer.Components.Convection convectionCatGas
    annotation (Placement(transformation(extent={{-10,34},{10,54}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor conductionHTCat(G=GOuter)
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={-30,70})));
  ThermoPower.Thermal.HT hT annotation (Placement(transformation(extent={{-40,
            84},{-20,104}}), iconTransformation(extent={{-10,-10},{10,10}})));
equation

  convectionCatGas.Gc  = alphaCat * A;

  connect(plenum.outlet, outfl)
    annotation (Line(points={{50,0},{80,0}},  color={159,159,223}));
  connect(outfl, outfl)
    annotation (Line(points={{80,0},{80,0}}, color={159,159,223}));
  connect(convectionCatGas.fluid, plenum.thermalPort)
    annotation (Line(points={{10,44},{30,44},{30,14}}, color={191,0,0}));
  connect(hT, conductionHTCat.port_a)
    annotation (Line(points={{-30,94},{-30,80}},      color={191,0,0}));
  connect(reactorNodeTest.Fuele_r, infl) annotation (Line(points={{-42.8,1.64},
          {-62.4,1.64},{-62.4,0},{-80,0}}, color={159,159,223}));
  connect(reactorNodeTest.Fuelr_mg, plenum.inlet) annotation (Line(points={{
          -17.2,2},{-4,2},{-4,0},{10,0}}, color={159,159,223}));
  connect(convectionCatGas.solid, reactorNodeTest.hT) annotation (Line(points={
          {-10,44},{-20,44},{-20,14.42},{-29.84,14.42}}, color={191,0,0}));
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}})),
    Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}),
        graphics={Bitmap(
          extent={{-60,-80},{60,80}},
          imageSource="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",
          fileName=
              "modelica://OpenTEMPEST/../../../../../../../../Desktop/Reaktor.png")}));
end Reactor0DSimpleMethanol;

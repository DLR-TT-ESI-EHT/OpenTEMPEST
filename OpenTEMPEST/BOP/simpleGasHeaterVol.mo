within OpenTEMPEST.BOP;
model simpleGasHeaterVol
  extends Flow.BaseClasses.simpleGasHeaterBase;

  import SI = Modelica.SIunits;
  parameter ThermoPower.Units.HydraulicResistance R=1 "Hydraulic resistance";
  parameter Boolean includePressureLoss=false;
  parameter Boolean noInitialPressure=false
    "Remove initial equation on pressure"
    annotation (Dialog(tab="Initialisation"),choices(checkBox=true));

  ThermoPower.Thermal.HT ht annotation (Placement(transformation(extent={{-10,54},
            {10,74}}), iconTransformation(extent={{-20,20},{20,60}})));

  ThermoPower.Gas.Plenum plenum(
    redeclare package Medium = medium,
    V=ACS*l,
    pstart=pStart,
    Tstart=TStart,
    Xstart=xStart,
    noInitialPressure=noInitialPressure)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
  Modelica.Thermal.HeatTransfer.Components.Convection convection annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={0,30})));

  Modelica.Blocks.Continuous.FirstOrder firstOrder(T=6)
    annotation (Placement(transformation(extent={{50,20},{30,40}})));
  ThermoPower.Gas.PressDropLin pressDropLin(redeclare package Medium = medium,
      R=R)
 if includePressureLoss
    annotation (Placement(transformation(extent={{22,-10},{42,10}})));
  ThermoPower.Gas.SensW sensW(redeclare package Medium = medium) if not
    includePressureLoss
    annotation (Placement(transformation(extent={{22,-30},{42,-10}})));

equation

  if not useAlphaIn then

    eta = medium.dynamicViscosity(plenum.gas.state); // 18.2e-6;//
    cp =  medium.specificHeatCapacityCp(plenum.gas.state); //  1;//
    lambda = medium.thermalConductivity(plenum.gas.state); // 0.0262; //0.0262; //

  else
    eta = -1;
    cp = -1;
    lambda = -1;
  end if;

  firstOrder.u = alpha .* A;

  connect(ht, convection.solid) annotation (Line(points={{0,64},{0,52},{1.77636e-15,
          52},{1.77636e-15,40}}, color={191,0,0}));
  connect(convection.fluid, plenum.thermalPort) annotation (Line(points={{-1.77636e-15,
          20},{-1.77636e-15,13.5},{0,13.5},{0,7}}, color={191,0,0}));
  connect(firstOrder.y, convection.Gc)
    annotation (Line(points={{29,30},{10,30}}, color={0,0,127}));
  connect(plenum.inlet, inlet)
    annotation (Line(points={{-10,0},{-90,0}}, color={159,159,223}));
   connect(plenum.outlet, pressDropLin.inlet)
    annotation (Line(points={{10,0},{22,0}}, color={159,159,223}));
  connect(pressDropLin.outlet, outlet) annotation (Line(points={{42,0},{54,0},{54,
          10},{90,10}}, color={159,159,223}));
  connect(sensW.inlet, plenum.outlet) annotation (Line(points={{26,-24},{16,-24},
          {16,0},{10,0}}, color={159,159,223}));
  connect(sensW.outlet, outlet) annotation (Line(points={{38,-24},{64,-24},{64,
          10},{90,10}}, color={159,159,223}));

  annotation (Icon(graphics={Rectangle(
          extent={{-60,20},{60,-20}},
          lineColor={28,108,200},
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid),
        Text(extent={{-100,-30},{100,-54}},  textString="%name")}));
end simpleGasHeaterVol;

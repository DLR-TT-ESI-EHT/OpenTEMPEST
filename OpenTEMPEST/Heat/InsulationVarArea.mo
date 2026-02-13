within OpenTEMPEST.Heat;
model InsulationVarArea
  extends Modelica.Icons.ObsoleteModel; //should be replaced by Solid0DVarArea
  import SI = Modelica.SIunits;
  parameter SI.ThermalConductivity lambda "thermalConductivity";
  parameter SI.SpecificHeatCapacity c "SpecificHeatCapacity";
  parameter SI.Density rho "Density";
  parameter SI.Area area1 "area1";
  parameter SI.Area area2 "area2";
  parameter SI.Length depth "depth";
  parameter SI.Temperature TStart "Start temperature";
  ThermoPower.Thermal.HT ht1 annotation (Placement(transformation(extent={{-12,40},
            {8,60}}),    iconTransformation(extent={{-12,40},{8,60}})));
  ThermoPower.Thermal.HT ht2 annotation (Placement(transformation(extent={{-10,-60},
            {10,-40}}),  iconTransformation(extent={{-10,-60},{10,-40}})));
  OpenTEMPEST.Heat.conductiveHT conductiveHT1(
    lambda=lambda,
    L=depth/2,
    A=area1*3/4 + area2*1/4)
    annotation (Placement(transformation(extent={{-10,34},{10,14}})));
  OpenTEMPEST.Heat.conductiveHT conductiveHT2(
    lambda=lambda,
    L=depth/2,
    A=area1*1/4 + area2*3/4)
    annotation (Placement(transformation(extent={{-10,-32},{10,-12}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor(C=c*rho*
        depth*1/2*(area1 + area2),                                                                 T(start=TStart))
    annotation (Placement(transformation(extent={{-10,-6},{10,14}})));
  ThermoPower.Thermal.HT htC annotation (Placement(transformation(extent={{-12,-8},
            {8,12}}),    iconTransformation(extent={{-12,-8},{8,12}})));
equation
  connect(heatCapacitor.port, conductiveHT1.side1)
    annotation (Line(points={{0,-6},{0,21}},           color={191,0,0}));
  connect(heatCapacitor.port, conductiveHT2.side1)
    annotation (Line(points={{0,-6},{0,-19}},
                                            color={191,0,0}));
  connect(conductiveHT1.side2, ht1) annotation (Line(points={{0,27.1},{0,50},
          {-2,50}},             color={191,0,0}));
  connect(ht2, conductiveHT2.side2) annotation (Line(points={{0,-50},{0,
          -25.1}},           color={191,0,0}));
  connect(ht1, ht1) annotation (Line(points={{-2,50},{-4,50},{-4,50},{-2,50}},
        color={191,0,0}));
  connect(heatCapacitor.port, htC)
    annotation (Line(points={{0,-6},{12,-6},{12,2},{-2,2}}, color={191,0,0}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Polygon(
          points={{-60,40},{60,40},{80,-40},{-80,-40},{-60,40}},
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid,
          pattern=LinePattern.None)}),                           Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end InsulationVarArea;

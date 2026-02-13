within OpenTEMPEST.Heat;
model Solid0D

  replaceable package SolidMat = TEMPEST.Solid.SolidMatBase  annotation(choicesAllMatching = true);
  SolidMat.BaseProperties Solid(T=T, each kCustom_trans=kCustom_trans, each kCustom_long=kCustom_long, each rhoCustom=rhoCustom, each cpCustom=cpCustom);

  import SI = Modelica.SIunits;

  // Dimensions
  parameter SI.Length lX = 1 "Length of Solid" annotation(Dialog(group="Dimensions"));
  parameter SI.Length lY = 1 "Width of Solid" annotation(Dialog(group="Dimensions"));
  parameter SI.Length lZ = 1 "Height of Solid" annotation(Dialog(group="Dimensions"));

  // Initial Values
  parameter SI.Temperature TStart=773.15 annotation(Dialog(group="Initialisation"));
  parameter Boolean force_der_T_Start=false annotation(Dialog(group="Initialisation"));
  parameter SI.TemperatureSlope der_T_Start=0 annotation(Dialog(group="Initialisation"));

  // Block Volume
  parameter Real dV = lX*lY*lZ "Volume of a Solid";

  // CV Values
  SI.Temperature T(start=TStart);

  // from Solid
  SI.ThermalConductivity k_trans = Solid.k_trans "Effective Thermal Conductivity across layers";
  SI.ThermalConductivity k_long= Solid.k_long "Effective Thermal Conductivity in plane of layers (=k_trans for homogeneous materials)";
  SI.Density rho = Solid.rho "Density";
  SI.SpecificHeatCapacity cp = Solid.cp "Specific Heat Capacity";

  parameter SI.ThermalConductivity kCustom_trans = 1 "Thermal Conductivity across layers" annotation(Dialog(group="Custom Material Only"));
  parameter SI.ThermalConductivity kCustom_long = 1 "Thermal Conductivity in plane of layers (=k_trans for homogeneous materials)" annotation(Dialog(group="Custom Material Only"));
  parameter SI.Density rhoCustom=1   "Density" annotation(Dialog(group="Custom Material Only"));
  parameter SI.SpecificHeatCapacity cpCustom = 1 "Specific Heat Capacity" annotation(Dialog(group="Custom Material Only"));

  // Thermal ports
  ThermoPower.Thermal.HT hT_x1 annotation (Placement(transformation(extent={{100,-10},
            {120,10}}), iconTransformation(extent={{100,-10},{120,10}})));
  ThermoPower.Thermal.HT hT_z1 annotation (Placement(transformation(extent={{-60,60},
            {-40,80}}),     iconTransformation(extent={{-60,60},{-40,80}})));
  ThermoPower.Thermal.HT hT_z0 annotation (Placement(transformation(extent={{-60,-80},
            {-40,-60}}), iconTransformation(extent={{-60,-80},{-40,-60}})));
  ThermoPower.Thermal.HT hT_y0 annotation (Placement(transformation(extent={{40,-80},
            {60,-60}}),    iconTransformation(extent={{40,-80},{60,-60}})));
  ThermoPower.Thermal.HT hT_y1 annotation (Placement(transformation(extent={{40,60},
            {60,80}}),  iconTransformation(extent={{40,60},{60,80}})));
  ThermoPower.Thermal.HT hT_int annotation (Placement(transformation(extent={{-10,-10},
            {10,10}}), iconTransformation(extent={{-10,-10},{10,10}})));

  ThermoPower.Thermal.HT hT_x0 annotation (Placement(transformation(extent={{-120,
            -10},{-100,10}}),
                         iconTransformation(extent={{-120,-10},{-100,10}})));
initial equation
  if force_der_T_Start then
  der(T)=der_T_Start;
  else
  T = TStart;
  end if;

equation

  // Energy Balance
  rho*cp*dV*der(T) = hT_x0.Q_flow + hT_x1.Q_flow + hT_z1.Q_flow +  hT_z0.Q_flow + hT_y1.Q_flow + hT_y0.Q_flow + hT_int.Q_flow;

  // Boundary Conditions
  hT_z0.Q_flow = (k_trans/(0.5*lZ))*lY*lX * (hT_z0.T - T); // Conductive Flow up
  hT_z1.Q_flow = (k_trans/(0.5*lZ))*lY*lX * (hT_z1.T - T); // Conductive Flow down
  hT_y0.Q_flow = (k_long/(0.5*lY))*lZ*lX * (hT_y0.T - T); // Conductive Flow Left
  hT_y1.Q_flow = (k_long/(0.5*lY))*lZ*lX * (hT_y1.T - T); // Conductive Flow Right
  hT_x0.Q_flow = (k_long/(0.5*lX))*lY*lZ * (hT_x0.T -T); // Conductive Flow In
  hT_x1.Q_flow = (k_long/(0.5*lX))*lY*lZ * (hT_x1.T -T); // Conductive Flow Out
  hT_int.T = T; // Internal HT

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Rectangle(
          extent={{-100,-20},{100,-60}},
          lineColor={28,108,200},
          fillColor={116,116,116},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-100,60},{100,20}},
          lineColor={28,108,200},
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-100,20},{100,-20}},
          lineColor={162,29,33},
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-58,62},{-40,38}},
          lineColor={238,46,47},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="z1"),
        Text(
          extent={{-58,-36},{-40,-66}},
          lineColor={238,46,47},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="z0"),
        Text(
          extent={{40,62},{62,40}},
          lineColor={238,46,47},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="y1"),
        Text(
          extent={{82,10},{98,-14}},
          lineColor={238,46,47},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="x1"),
        Text(
          extent={{42,-38},{62,-60}},
          lineColor={238,46,47},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="y0"),
        Text(
          extent={{-98,12},{-82,-10}},
          lineColor={238,46,47},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="x0"),
          Text(
            extent={{-134,-12},{138,-50}},
            lineColor={255,255,255},
          textString="%name")}),                                    Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end Solid0D;

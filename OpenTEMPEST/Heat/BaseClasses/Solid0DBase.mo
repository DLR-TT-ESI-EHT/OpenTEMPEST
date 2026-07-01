within OpenTEMPEST.Heat.BaseClasses;
partial model Solid0DBase

  import SI = Modelica.SIunits;
  replaceable package SolidMat = OpenTEMPEST.Solid.SolidMatBase  annotation(choicesAllMatching = true);

  // Initial Values
  parameter SI.Temperature TStartbar=773.15 "Uniform initial temperature of the solid" annotation(Dialog(group="Initialisation"));
  parameter Boolean force_der_T_Start=false annotation(Dialog(group="Initialisation"));
  parameter SI.TemperatureSlope der_T_Start=0 annotation(Dialog(group="Initialisation"));

  // Dimensions
  parameter SI.Length lX = 1 "Total length of solid" annotation(Dialog(group="Dimensions"));
  parameter SI.Length lY = 1 "Total width of solid" annotation(Dialog(group="Dimensions"));
  parameter SI.Length lZ = 1 "Total thickness of solid" annotation(Dialog(group="Dimensions"));

  parameter Real dV = lX*lY*lZ "Volume of the solid"  annotation(Dialog(group="Dimensions Extra"));
  parameter SI.Area Ax = lY*lZ "Cross sectional area in x-plane"  annotation(Dialog(group="Dimensions Extra"));
  parameter SI.Area Ay = lX*lZ "Cross sectional area in y-plane"  annotation(Dialog(group="Dimensions Extra"));
  parameter SI.Area Az = lX*lY "Cross sectional area in z-plane"  annotation(Dialog(group="Dimensions Extra"));

  SI.Temperature T(start=TStartbar) "Temperature";

  SolidMat.BaseProperties Solid(
    T=T,
    kCustom_trans=kCustom_trans,
    kCustom_long=kCustom_long,
    rhoCustom=rhoCustom,
    cpCustom=cpCustom) "Material properties object";

  // Solid properties
  SI.ThermalConductivity k_trans = Solid.k_trans "Effective Thermal Conductivity across layers";
  SI.ThermalConductivity k_long= Solid.k_long "Effective Thermal Conductivity in plane of layers (=k_trans for homogeneous materials)";
  SI.Density rho = Solid.rho "Density";
  SI.SpecificHeatCapacity cp = Solid.cp "Specific heat capacity";

  SI.HeatFlowRate Qext "User-defined heat source or sink for each CV";

  // Custom material properties
  parameter SI.ThermalConductivity kCustom_trans = 1 "Thermal Conductivity across layers" annotation(Dialog(group="Custom Material Only"));
  parameter SI.ThermalConductivity kCustom_long = 1 "Thermal Conductivity in plane of layers (=k_trans for homogeneous materials)" annotation(Dialog(group="Custom Material Only"));
  parameter SI.Density rhoCustom = 1   "Density" annotation(Dialog(group="Custom Material Only"));
  parameter SI.SpecificHeatCapacity cpCustom = 1 "Specific Heat Capacity" annotation(Dialog(group="Custom Material Only"));

  // Thermal ports
  ThermoPower.Thermal.HT hT_xN annotation (Placement(transformation(extent={{100,-10},
            {120,10}}), iconTransformation(extent={{100,-10},{120,10}})));
  ThermoPower.Thermal.HT hT_zN annotation (Placement(transformation(extent={{-60,60},
            {-40,80}}),     iconTransformation(extent={{-60,60},{-40,80}})));
  ThermoPower.Thermal.HT hT_z0 annotation (Placement(transformation(extent={{-60,-80},
            {-40,-60}}), iconTransformation(extent={{-60,-80},{-40,-60}})));
  ThermoPower.Thermal.HT hT_y0 annotation (Placement(transformation(extent={{40,-80},
            {60,-60}}),    iconTransformation(extent={{40,-80},{60,-60}})));
  ThermoPower.Thermal.HT hT_yN annotation (Placement(transformation(extent={{40,60},
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
  T = TStartbar;
  end if;

equation

  // Energy Balance
  rho*cp*dV*der(T) = hT_x0.Q_flow + hT_xN.Q_flow + hT_zN.Q_flow +  hT_z0.Q_flow + hT_yN.Q_flow + hT_y0.Q_flow + hT_int.Q_flow;

  // Boundary Conditions
  hT_z0.Q_flow = (k_trans/(0.5*lZ))*lY*lX * (hT_z0.T - T); // Conductive Flow up
  hT_zN.Q_flow = (k_trans/(0.5*lZ))*lY*lX * (hT_zN.T - T); // Conductive Flow down
  hT_y0.Q_flow = (k_long/(0.5*lY))*lZ*lX * (hT_y0.T - T); // Conductive Flow Left
  hT_yN.Q_flow = (k_long/(0.5*lY))*lZ*lX * (hT_yN.T - T); // Conductive Flow Right
  hT_x0.Q_flow = (k_long/(0.5*lX))*lY*lZ * (hT_x0.T -T); // Conductive Flow In
  hT_xN.Q_flow = (k_long/(0.5*lX))*lY*lZ * (hT_xN.T -T); // Conductive Flow Out
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
          textString="zN"),
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
          textString="yN"),
        Text(
          extent={{82,10},{98,-14}},
          lineColor={238,46,47},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="xN"),
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
        coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<h2>Solid0DBase</h2>
<p>
The model represents a zero-dimensional (lumped) solid with a single control volume.
It computes the transient temperature of the solid based on conduction through all 
six faces (x0, xN, y0, yN, z0, zN), as well as a user-defined external heat flux. 
The model supports anisotropic and custom material properties.
</p>
</html>"));
end Solid0DBase;

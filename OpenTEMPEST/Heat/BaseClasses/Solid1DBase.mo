within OpenTEMPEST.Heat.BaseClasses;
partial model Solid1DBase

  import SI = Modelica.SIunits;
  replaceable package SolidMat = OpenTEMPEST.Solid.SolidMatBase annotation(choicesAllMatching = true);

  parameter Integer N(min=3) = 5 "Number of CVs in the solid";

  // Initial Values
  parameter SI.Temperature Tstartbar = 1073.15 "Uniform initial temperature of the solid" annotation(Dialog(group="Initialisation"));

  // Dimensions
  parameter SI.Length lX = 1  "Total length of solid" annotation(Dialog(group="Dimensions"));
  parameter SI.Length lY = 1 "Total width of solid" annotation(Dialog(group="Dimensions"));
  parameter SI.Length lZ = 1 "Total thickness of solid" annotation(Dialog(group="Dimensions"));

  parameter SI.Length dx =  lX /N "x-Length of a CV" annotation(Dialog(group="Dimensions Extra"));
  parameter SI.Volume dV = Ax*dx "Volume of a CV" annotation(Dialog(group="Dimensions Extra"));
  parameter SI.Area Ax = lY*lZ "Cross sectional area in x plane" annotation(Dialog(group="Dimensions Extra"));
  parameter SI.Area Ay = dx*lZ "Cross sectional area in y-plane"  annotation(Dialog(group="Dimensions Extra"));
  parameter SI.Area Az = dx*lY "Cross sectional area in z-plane"  annotation(Dialog(group="Dimensions Extra"));

  SI.Temperature T[N] "Temperature in each CV";

  SolidMat.BaseProperties Solid[N](
    T=T,
    each kCustom_trans=kCustom_trans,
    each kCustom_long=kCustom_long,
    each rhoCustom=rhoCustom,
    each cpCustom=cpCustom) "Material properties object for each CV";

  // Solid Properties
  SI.ThermalConductivity k_trans[N] = Solid[:].k_trans "Thermal conductivities in CV";
  SI.ThermalConductivity k_long[N] = Solid[:].k_long "Thermal conductivities in CV";
  SI.Density rho[N] = Solid[:].rho "Density of CV";
  SI.SpecificHeatCapacity cp[N] = Solid[:].cp "heat capacity of CV";

  SI.ThermalConductivity kv_long[N-1] "Harmonic mean of thermal conductivities at internal CV boundaries";  // https://doi.org/10.1016/S1018-3639(18)30628-7

  parameter SI.ThermalConductivity kCustom_trans = 1 annotation(Dialog(group="Custom Material Only"));
  parameter SI.ThermalConductivity kCustom_long = kCustom_trans annotation(Dialog(group="Custom Material Only"));
  parameter SI.Density rhoCustom=1   annotation(Dialog(group="Custom Material Only"));
  parameter SI.SpecificHeatCapacity cpCustom = 1 annotation(Dialog(group="Custom Material Only"));

  SI.HeatFlowRate Qext[N];

  // Thermal ports
  ThermoPower.Thermal.HT hT_x0 annotation (Placement(transformation(extent={{-120,-10},{-100,10}}),
                         iconTransformation(extent={{-120,-10},{-100,10}})));
  ThermoPower.Thermal.HT hT_xN annotation (Placement(transformation(extent={{100,-10},
            {120,10}}), iconTransformation(extent={{100,-10},{120,10}})));
  ThermoPower.Thermal.DHTVolumes dhT_z1(N=N) annotation (Placement(transformation(extent={{-80,60},
            {-60,80}}), iconTransformation(extent={{-80,60},{-60,80}})));
  ThermoPower.Thermal.DHTVolumes dhT_z0(N=N) annotation (Placement(transformation(extent={{-80,-80},
            {-60,-60}}), iconTransformation(extent={{-80,-80},{-60,-60}})));
  ThermoPower.Thermal.DHTVolumes dhT_y0(N=N) annotation (Placement(transformation(extent={{60,-80},
            {80,-60}}),iconTransformation(extent={{60,-80},{80,-60}})));
  ThermoPower.Thermal.DHTVolumes dhT_y1(N=N) annotation (Placement(transformation(extent={{60,60},
            {80,80}}),  iconTransformation(extent={{60,60},{80,80}})));
  ThermoPower.Thermal.DHTVolumes dhT_int(N=N)   annotation (Placement(transformation(extent={{-10,-6},
            {10,14}}), iconTransformation(extent={{-10,-6},{10,14}})));

initial equation
  T = fill(Tstartbar, N);

equation
  kv_long[:]  = 2.*k_long[1:N-1] .*k_long[2:N] ./(k_long[1:N-1]  .+ k_long[2:N]);

  // Energy Balance
   rho[1]*cp[1]*dV*der(T[1])                = hT_x0.Q_flow + kv_long[1]*Ax*(T[2] - T[1])/dx + dhT_z1.Q[1] +  dhT_z0.Q[1] + dhT_y1.Q[1] + dhT_y0.Q[1] + dhT_int.Q[1] + Qext[1];
   rho[2:N-1].*cp[2:N-1].*dV.*der(T[2:N-1]) = kv_long[1:N-2].*Ax.*(T[1:N-2] .- T[2:N-1])./dx .+ kv_long[2:N-1].*Ax.*(T[3:N] .- T[2:N-1])./dx .+ dhT_z1.Q[2:N-1] .+ dhT_z0.Q[2:N-1] .+ dhT_y1.Q[2:N-1] .+ dhT_y0.Q[2:N-1] .+ dhT_int.Q[2:N-1] + Qext[2:N-1];
   rho[N]*cp[N]*dV*der(T[N])                = hT_xN.Q_flow + kv_long[N-1]*Ax*(T[N-1] - T[N])/dx + dhT_z1.Q[N] + dhT_z0.Q[N] + dhT_y1.Q[N] + dhT_y0.Q[N] + dhT_int.Q[N] + Qext[N];

  // Boundary Conditions
  dhT_z0.Q[:] = (k_trans./(0.5*lZ))*lY*dx .* (dhT_z0.T[:] .- T[:]); // Conductive Flow up
  dhT_z1.Q[:] = (k_trans./(0.5*lZ))*lY*dx .* (dhT_z1.T[:] .- T[:]); // Conductive Flow down
  dhT_y0.Q[:] = (k_long./(0.5*lY))*lZ*dx .* (dhT_y0.T[:] .- T[:]); // Conductive Flow Left
  dhT_y1.Q[:] = (k_long./(0.5*lY))*lZ*dx .* (dhT_y1.T[:] .- T[:]); // Conductive Flow Right
  hT_x0.Q_flow = (k_long[1]/(0.5*dx))*Ax * (hT_x0.T - T[1]); // Conductive Flow In
  hT_xN.Q_flow = (k_long[N]/(0.5*dx))*Ax * (hT_xN.T - T[N]); // Conductive Flow Out
  dhT_int.T[:] = T[:];

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
          extent={{-78,62},{-58,42}},
          lineColor={244,125,35},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="z1"),
        Text(
          extent={{-78,-38},{-58,-62}},
          lineColor={244,125,35},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="z0"),
        Text(
          extent={{62,-34},{80,-64}},
          lineColor={244,125,35},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="y0"),
        Text(
          extent={{62,62},{80,38}},
          lineColor={244,125,35},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="y1"),
        Text(
          extent={{-98,12},{-82,-12}},
          lineColor={238,46,47},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="x0"),
        Text(
          extent={{82,12},{98,-10}},
          lineColor={238,46,47},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="xN"),
          Text(
            extent={{-140,-16},{140,-56}},
            lineColor={255,255,255},
          textString="%name")}),                                    Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    Documentation(revisions="<html>
</html>

", info="<html>
<h2>Solid1DBase</h2>
<p>
The model represents a one-dimensional solid discretized into control volumes (CVs)
in the x-direction. It computes the transient temperature distribution 
within the solid based on conduction in the x-direction, as well as heat exchange 
in the y- and z-directions through thermal ports and user-defined external 
heat fluxes. The model supports anisotropic and custom material properties.
</p>
</html>"));
end Solid1DBase;

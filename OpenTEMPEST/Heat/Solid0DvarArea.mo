within OpenTEMPEST.Heat;
model Solid0DvarArea
  "solid 1D volume with area variation along the x-axis (truncated pyramid with base area in y-z-plane)"

  import SI = Modelica.SIunits;
  replaceable package SolidMat = TEMPEST.Solid.SolidMatBase annotation(choicesAllMatching = true);

  SolidMat.BaseProperties Solid(T=T, kCustom_trans=kCustom, kCustom_long=kCustom, rhoCustom=rhoCustom, cpCustom=cpCustom);

  // Initial Values
  parameter SI.Temperature TstartIn=773.15 annotation(Dialog(group="Initialisation"));
  parameter SI.Temperature TstartOut=773.15 annotation(Dialog(group="Initialisation"));
  parameter SI.Temperature Tstart = 773.15 annotation(Dialog(group="Initialisation"));

  // Dimensions
  parameter SI.Length lX = 1 "Total length of solid" annotation(Dialog(group="Dimensions"));

  parameter SI.Area A0 = 1 "CS area at x/L=0" annotation(Dialog(group="Dimensions"));
  parameter SI.Area A1 = 1 "CS area at x/L=1" annotation(Dialog(group="Dimensions"));

  SI.Volume dV "Volume of CV";

  SI.Temperature T "Average Temperature in CV";

  // Solid Properties
  SI.ThermalConductivity k = Solid.k_long "Thermal Conductivity";
  SI.Density rho = Solid.rho "Density";
  SI.SpecificHeatCapacity cp = Solid.cp "Specific Heat Capacity";

  parameter SI.ThermalConductivity kCustom = 1 annotation(Dialog(group="Custom Material Only"));
  parameter SI.Density rhoCustom = 1 annotation(Dialog(group="Custom Material Only"));
  parameter SI.SpecificHeatCapacity cpCustom = 1 annotation(Dialog(group="Custom Material Only"));

  // Thermal ports
  ThermoPower.Thermal.HT hT_x0 annotation (Placement(transformation(extent={{-120,
            -10},{-100,10}}),
                         iconTransformation(extent={{-120,-10},{-100,10}})));
  ThermoPower.Thermal.HT hT_x1 annotation (Placement(transformation(extent={{100,-10},
            {120,10}}), iconTransformation(extent={{100,-10},{120,10}})));

initial equation
  T = Tstart;

equation
  //cross sectional area in x-plane
  dV = ( A0 + A1 + sqrt(A0 * A1)) *lX/3;

  // Energy Balance
  rho*cp*dV*der(T) = hT_x0.Q_flow + hT_x1.Q_flow;

  // Boundary Conditions
  hT_x0.Q_flow = (k/(0.5*lX))*A0 * (hT_x0.T -T); // Conductive Flow In
  hT_x1.Q_flow = (k/(0.5*lX))*A1 * (hT_x1.T -T); // Conductive Flow Out

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Polygon(
          points={{-100,40},{100,80},{100,-80},{-100,-40},{-100,12},{-100,40}},
          lineColor={28,108,200},
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-98,12},{-82,-12}},
          lineColor={238,46,47},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="x0"),
        Text(
          extent={{80,12},{98,-12}},
          lineColor={238,46,47},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="xN"),
          Text(
            extent={{-139,23},{139,-23}},
            lineColor={255,255,255},
          textString="%name",
          origin={43,-1},
          rotation=90)}),                                           Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    Documentation(revisions="<html>
<ul>
<li><i>03 Aug 2021</i> by <a href=\"faisal.sedeqi@dlr.de\">Faisal Sedeqi</a>:<br>Modified Energy balance to use harmonic mean of k at interface between CV.</li>
<li><i>27 Jul 2021</i> by <a href=\"hans.wiggenhauser@dlr.de\">Hans Wiggenhauser</a> and <a href=\"faisal.sedeqi@dlr.de\">Faisal Sedeqi</a>:<br>First release. </li>
</ul>
</html>", info="<html>
<p>truncated pyramid volume with heat transfer along x-direction. Pyramid base area is y0*z0.</p>
</html>"));
end Solid0DvarArea;

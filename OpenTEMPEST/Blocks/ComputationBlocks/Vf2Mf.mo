within OpenTEMPEST.Blocks.ComputationBlocks;
model Vf2Mf "Converts a volume flow into a mass flow"
  import SI = Modelica.SIunits;

  parameter SI.Temperature TRef = 273.15;
  parameter SI.Pressure pRef = 101325;
  parameter SI.MolarMass M;

  Modelica.Blocks.Interfaces.RealInput u
    annotation (Placement(transformation(extent={{-80,-20},{-40,20}})));
  Modelica.Blocks.Interfaces.RealOutput y
    annotation (Placement(transformation(extent={{40,-20},{80,20}})));

equation
  y = pRef*u*M/Modelica.Constants.R/TRef/1000/60;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
          extent={{-40,40},{40,-40}},
          lineColor={28,108,200},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid)}), Diagram(coordinateSystem(
          preserveAspectRatio=false)),
    Documentation(revisions="<html>
<ul>
<li>&quot;before versioning&quot; by Marius Tomberg: <br>Created</li>
</ul>
</html>"));
end Vf2Mf;

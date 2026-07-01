within OpenTEMPEST.Heat;
model Solid1D
  extends OpenTEMPEST.Heat.BaseClasses.Solid1DBase;

equation

  Qext = zeros(N);

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

"));
end Solid1D;

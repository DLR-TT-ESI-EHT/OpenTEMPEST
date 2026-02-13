within OpenTEMPEST.Heat;
model Solid2D
  extends OpenTEMPEST.Heat.BaseClasses.Solid2DBase;

equation
  Qext[:] = zeros(nX, nY);

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Rectangle(
          extent={{-100,60},{100,20}},
          lineColor={28,108,200},
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-100,-20},{100,-60}},
          lineColor={28,108,200},
          fillColor={116,116,116},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-100,20},{100,-20}},
          lineColor={162,29,33},
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-76,64},{-62,40}},
          lineColor={255,240,0},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="z1"),
        Text(
          extent={{-78,-38},{-60,-62}},
          lineColor={255,240,0},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="z0"),
        Text(
          extent={{62,-34},{80,-62}},
          lineColor={244,125,35},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="y0"),
        Text(
          extent={{62,62},{80,40}},
          lineColor={244,125,35},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="yN"),
          Text(
            extent={{-138,-18},{142,-58}},
            lineColor={255,255,255},
          textString="%name"),
        Text(
          extent={{-96,12},{-78,-10}},
          lineColor={244,125,35},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="x0"),
        Text(
          extent={{82,10},{98,-10}},
          lineColor={244,125,35},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="xN")}),                                         Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    Documentation(revisions="<html>
<ul>
<li><i>27 Aug 2021</i> by <a href=\"faisal.sedeqi@dlr.de\">Faisal Sedeqi</a>:<br>Modified Conduction Between CVs including Direction Dependent Conductivity. </li>
<li><i>03 Aug 2021</i> by <a href=\"faisal.sedeqi@dlr.de\">Faisal Sedeqi</a>:<br>Changed port names. Will change conduction between CVs next.</li>
<li><i>23 Jul 2021</i> by <a href=\"hans.wiggenhauser@dlr.de\">Hans Wiggenhauser</a>:<br>finished and tested.</li>
<li><i>20 Jul 2021</i> by <a href=\"faisal.sedeqi@dlr.de\">Faisal Sedeqi</a>:<br>created, edge energy balance missinge. </li>
</ul>
</html>"));
end Solid2D;

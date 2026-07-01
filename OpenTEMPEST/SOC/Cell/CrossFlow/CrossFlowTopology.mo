within OpenTEMPEST.SOC.Cell.CrossFlow;
model CrossFlowTopology

    parameter Integer nX(min=1)=5;
    parameter Integer nY(min=1)=5;
    parameter Integer nSpecies=2;
    parameter Boolean includeVarStream=true "True if VariablesStream is used";

  Heat.DHTVolumes2D dHTT_side1(i=nX, j=nY) "nX * nY" annotation (Placement(
        transformation(extent={{-80,40},{80,60}}), iconTransformation(extent={{
            -80,40},{80,60}})));
  Heat.DHTVolumes2D dHTT_side2(i=nY, j=nX) "nY * nX" annotation (Placement(
        transformation(extent={{-80,-60},{80,-40}}), iconTransformation(extent=
            {{-80,-60},{80,-40}})));

  OpenTEMPEST.SOC.Electrochem.Interfaces.VariablesStream PEN_side1[nX,nY](each
      nspecies=nSpecies) if includeVarStream "nX*nY" annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-90,50})));
  OpenTEMPEST.SOC.Electrochem.Interfaces.VariablesStream PEN_side2[nY,nX](each
      nspecies=nSpecies) if includeVarStream "nY*nX" annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-90,-50})));

equation

  if includeVarStream then
    for i in 1:nX loop
      for j in 1:nY loop
        // Clockwise rotation
        dHTT_side1.T[i,j] = dHTT_side2.T[j, nX+1-i];
        dHTT_side1.Q[i,j] = -dHTT_side2.Q[j, nX+1-i];
        connect(PEN_side1[i,j], PEN_side2[j, nX+1-i]);
      end for;
    end for;
  else
    for i in 1:nX loop
      for j in 1:nY loop
        // Counter clockwise rotation
         dHTT_side1.T[i,j] = dHTT_side2.T[nY+1-j, i];
         dHTT_side1.Q[i,j] =-dHTT_side2.Q[nY+1-j, i];
      end for;
    end for;
  end if;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,40},{100,-40}},
          lineColor={28,108,200},
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid), Line(
          points={{23,16},{-17,16},{-17,-16},{-9,-8},{-17,-16},{-23,-10}},
          color={238,46,47},
          thickness=1,
          smooth=Smooth.Bezier,
          origin={-23,6},
          rotation=0)}),                                         Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<h2>CrossFlowTopology</h2>
<p>
Structural model that maps the discretised temperature and heat-flux fields between two
heat-transfer surfaces arranged in a cross-flow configuration. The model performs no
thermodynamic calculations; it only rewrites array indices to couple the two channel sides.
</p>
<ul>
<li>Side 1 is discretised on an <code>nX &times; nY</code> grid; side 2 on a transposed <code>nY &times; nX</code> grid.</li>
<li>When <code>includeVarStream = true</code>, a <b>clockwise</b> 90&deg; index rotation is applied and the
    <code>VariablesStream</code> connectors (<code>PEN_side1</code>, <code>PEN_side2</code>) are instantiated and connected.</li>
<li>When <code>includeVarStream = false</code>, a <b>counter-clockwise</b> rotation is used and the
    <code>VariablesStream</code> connectors are removed entirely (conditional components).</li>
<li>Heat-flux sign is inverted across sides to enforce energy conservation:
    <code>Q_side1[i,j] = -Q_side2[j, nX+1-i]</code>.</li>
</ul>
<h3>Index Mapping</h3>
<ul>
<li><b>CW</b> (<code>includeVarStream = true</code>):&nbsp;
    <code>T_side1[i,j] = T_side2[j, nX+1-i]</code></li>
<li><b>CCW</b> (<code>includeVarStream = false</code>):&nbsp;
    <code>T_side1[i,j] = T_side2[nY+1-j, i]</code></li>
</ul>
</html>
"));
end CrossFlowTopology;

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
          rotation=0)}),                                         Diagram(coordinateSystem(preserveAspectRatio=false)));
end CrossFlowTopology;

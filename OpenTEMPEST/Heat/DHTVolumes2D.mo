within OpenTEMPEST.Heat;
connector DHTVolumes2D "Distributed Heat Terminal"
  parameter Integer i=1 "number of volumes in direction 1";
  parameter Integer j=1 "number of volumes in direction 2";
//   parameter Integer N[i,j] "Number of volumes";
  Modelica.SIunits.Temperature T[i,j] "Temperature at the volumes";
  flow Modelica.SIunits.Power Q[i,j] "Heat flow at the volumes";
  annotation (
          Diagram(coordinateSystem(preserveAspectRatio=false)),
          Icon(graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillColor={255,240,0},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-46,42},{54,-42}},
          lineColor={0,0,0},
          textString="2D")}));
end DHTVolumes2D;

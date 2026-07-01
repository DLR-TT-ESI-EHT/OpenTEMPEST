within OpenTEMPEST.Heat.BaseClasses;
partial model HeatTransferNodeBase

  parameter Integer N = 5;

  ThermoPower.Thermal.DHTVolumes coldWall(N=N)
    annotation (Placement(transformation(extent={{-40,40},{40,60}})));
  ThermoPower.Thermal.DHTVolumes hotWall(N=N)
    annotation (Placement(transformation(extent={{-40,-40},{40,-20}})));

equation

   hotWall.Q = -coldWall.Q;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
          extent={{-100,40},{100,-20}},
          lineColor={28,108,200},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid), Text(
          extent={{-26,22},{28,-6}},
          lineColor={28,108,200},
          fillColor={255,255,255},
          fillPattern=FillPattern.None,
          textString="HT Node")}),                               Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end HeatTransferNodeBase;

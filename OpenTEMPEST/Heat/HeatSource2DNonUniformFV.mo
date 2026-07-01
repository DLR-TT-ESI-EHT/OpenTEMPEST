within OpenTEMPEST.Heat;
model HeatSource2DNonUniformFV
  "Distributed Heat Flow Source for Finite Volume models with non-uniformly distributed flow"
  extends ThermoPower.Icons.HeatFlow;
  parameter Integer nX = 1 "Number of volumes x-direction";
  parameter Integer nY = 1 "Number of volumes y-direction";
  OpenTEMPEST.Heat.DHTVolumes2D wall(i=nX, j=nY) annotation (Placement(
        transformation(extent={{-40,-40},{40,-20}}, rotation=0)));
  Modelica.Blocks.Interfaces.RealInput power[nX, nY] annotation (Placement(
        transformation(
        origin={0,40},
        extent={{-20,-20},{20,20}},
        rotation=270)));
equation
  wall.Q = -power;
  annotation (
    Diagram(graphics),
    Icon(graphics={Text(
          extent={{-100,-44},{100,-68}},
          lineColor={191,95,0},
          textString="%name")}),
    Documentation(info=""));
end HeatSource2DNonUniformFV;

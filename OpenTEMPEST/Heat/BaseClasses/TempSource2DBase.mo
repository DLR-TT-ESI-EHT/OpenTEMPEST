within OpenTEMPEST.Heat.BaseClasses;
partial model TempSource2DBase
  "Base Model Distributed Temperature Source for 2D Finite Volume models"
  extends ThermoPower.Icons.HeatFlow;
  parameter Integer nX = 1 "Number of volumes x-direction";
  parameter Integer nY = 1 "Number of volumes y-direction";

  OpenTEMPEST.Heat.DHTVolumes2D wall(i=nX, j=nY) annotation (Placement(
        transformation(extent={{-40,-40},{40,-20}}, rotation=0)));

  annotation (
    Icon(graphics={Text(
          extent={{-100,-44},{100,-68}},
          lineColor={191,95,0},
          textString="%name")}),
    Documentation(info=""));
end TempSource2DBase;

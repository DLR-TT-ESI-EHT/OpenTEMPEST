within OpenTEMPEST.Heat.BaseClasses;
partial model TempSource2DBase "Base model for distributed temperature sources for 2D finite volume models"
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

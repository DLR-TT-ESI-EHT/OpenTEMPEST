within OpenTEMPEST.Heat;
model TempSource2DNonUniform
  "Non-uniformly Distributed Temperature Source for 2D Finite Volume models"
  extends OpenTEMPEST.Heat.BaseClasses.TempSource2DBase;

  Modelica.Blocks.Interfaces.RealInput temperature[nX, nY] annotation (Placement(
         transformation(
         origin={0,40},
         extent={{-20,-20},{20,20}},
         rotation=270)));

equation

wall.T = temperature;

end TempSource2DNonUniform;

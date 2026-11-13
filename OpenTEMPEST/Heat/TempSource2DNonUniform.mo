within OpenTEMPEST.Heat;
model TempSource2DNonUniform "Non-uniformly distributed temperature source for 2D finite volume models"
  extends OpenTEMPEST.Heat.BaseClasses.TempSource2DBase;

  Modelica.Blocks.Interfaces.RealInput temperature[nX, nY] annotation (Placement(
         transformation(
         origin={0,40},
         extent={{-20,-20},{20,20}},
         rotation=270)));

equation

wall.T = temperature;

end TempSource2DNonUniform;

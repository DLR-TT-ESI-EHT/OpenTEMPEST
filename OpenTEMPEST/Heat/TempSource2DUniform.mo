within OpenTEMPEST.Heat;
model TempSource2DUniform
  "Uniformly Distributed Temperature Source for 2D Finite Volume models"
  extends OpenTEMPEST.Heat.BaseClasses.TempSource2DBase;

  Modelica.Blocks.Interfaces.RealInput temperature( unit="K", displayUnit = "degC") "Temperature [K]" annotation (Placement(
         transformation(
         origin={0,40},
         extent={{-20,-20},{20,20}},
         rotation=270)));

equation

wall.T = fill(temperature, nX, nY);

end TempSource2DUniform;

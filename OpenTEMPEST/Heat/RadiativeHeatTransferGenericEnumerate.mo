within OpenTEMPEST.Heat;
model RadiativeHeatTransferGenericEnumerate
  "0D radiative heat transfer with A2>=A1 & A1 is convex or flat"
  import SI = Modelica.SIunits;
  parameter SI.SpectralEmissivity epsilon1 "emissivity of surface 1";
  parameter SI.SpectralEmissivity epsilon2 "emissivity of surface 2";
  parameter SI.Area A1 "Area of surface 1, used for special case";
  parameter SI.Area A2 "Area of surface 2";
  //parameter Boolean specialCase1=false "A1 = A2 = inf & surfaces are flat and parallel";
  //parameter Boolean specialCase2=false "A2>>A1";

  replaceable parameter Enumerations.RadiationSpecialCases specialCase=
      Enumerations.RadiationSpecialCases.case3 "Radiation Special Case";

  ThermoPower.Thermal.HT side1 annotation (
    Placement(transformation(extent = {{-26, 2}, {8, 36}}), iconTransformation(extent = {{-8, 20}, {8, 36}})));
  ThermoPower.Thermal.HT side2 annotation (
    Placement(transformation(extent = {{-26, -2}, {8, -36}}), iconTransformation(extent = {{-8, -20}, {8, -36}})));
equation

  if specialCase == Enumerations.RadiationSpecialCases.case1 then
    side1.Q_flow = Modelica.Constants.sigma * A1 * 1/(1/epsilon1+1/epsilon2-1) * ( sign(side1.T) * side1.T ^ 4 - sign(side2.T) * side2.T ^ 4);
  elseif specialCase == Enumerations.RadiationSpecialCases.case2 then
    side1.Q_flow = Modelica.Constants.sigma * A1 * epsilon1 * ( sign(side1.T) * side1.T ^ 4 - sign(side2.T) * side2.T ^ 4);
  else
    side1.Q_flow = Modelica.Constants.sigma * A1 / (1/epsilon1+A1/A2*(1/epsilon2-1)) * ( sign(side1.T) * side1.T ^ 4 - sign(side2.T) * side2.T ^ 4);
  end if;

  side1.Q_flow = -side2.Q_flow;

  annotation (
    Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}})),
    Icon(coordinateSystem(preserveAspectRatio=false,   extent={{-100,-100},
            {100,100}}),                                                                     graphics={  Rectangle(extent = {{-80, 20}, {80, -20}}, lineColor = {0, 0, 0}, fillColor = {135, 135, 135},
            fillPattern =                                                                                                   FillPattern.CrossDiag)}));
end RadiativeHeatTransferGenericEnumerate;

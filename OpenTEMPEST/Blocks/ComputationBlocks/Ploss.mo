within OpenTEMPEST.Blocks.ComputationBlocks;
model Ploss "Power loss of the power electronics"
extends Modelica.Blocks.Interfaces.SISO;

  parameter Modelica.SIunits.Power Pdesign = 2377;
  parameter Real a = 0.048212553240941;
  parameter Real b = 0.009846461934665;
  parameter Real c = 0.011514564183595;

//protected
  Real PlossRel;
  Real x;
  Modelica.SIunits.Power P;
  Modelica.SIunits.Power Ploss;

equation

  u = P;
  y = Ploss;
  x = P/Pdesign;

  PlossRel = a * x^2 + b*x + c;
  PlossRel = Ploss/Pdesign;

end Ploss;

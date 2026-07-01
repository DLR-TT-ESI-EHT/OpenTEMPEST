within OpenTEMPEST.Blocks.Functions;
function g0_T
  "Calculating the gibbs energy at given Temperature and Standard pressure of a gas
  
  Example: G_h2o = TEMPEST.Blocks.Functions.g0_T(T=T, data=Modelica.Media.IdealGases.Common.SingleGasesData.H2O);
  "
  extends Modelica.Icons.Function;
  import SI = Modelica.SIunits;
  input Modelica.Media.IdealGases.Common.DataRecord data "Ideal gas data";
  input SI.Temperature T "Temperature";
  output SI.SpecificGibbsFreeEnergy g;
algorithm
  g := Modelica.Media.IdealGases.Common.Functions.h_T(data, T = T, exclEnthForm = false) - T * Modelica.Media.IdealGases.Common.Functions.s0_T(data, T = T);
end g0_T;

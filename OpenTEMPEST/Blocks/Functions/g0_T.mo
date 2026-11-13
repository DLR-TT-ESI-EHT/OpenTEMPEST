within OpenTEMPEST.Blocks.Functions;
function g0_T "Gibbs energy of a gas at temperature T and standard pressure"
  extends Modelica.Icons.Function;
  import SI = Modelica.SIunits;
  input Modelica.Media.IdealGases.Common.DataRecord data "Ideal gas data";
  input SI.Temperature T "Temperature";
  output SI.SpecificGibbsFreeEnergy g;
algorithm
  g := Modelica.Media.IdealGases.Common.Functions.h_T(data, T = T, exclEnthForm = false) - T * Modelica.Media.IdealGases.Common.Functions.s0_T(data, T = T);
  annotation (Documentation(info="<html>
<p>Returns the molar Gibbs energy of a single gas at temperature <code>T</code> and standard pressure.</p>
<p>Example:</p>
<pre>G_h2o = g0_T(T=T, data=Modelica.Media.IdealGases.Common.SingleGasesData.H2O);</pre>
</html>"));
end g0_T;

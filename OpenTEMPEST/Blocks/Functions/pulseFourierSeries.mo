within OpenTEMPEST.Blocks.Functions;
function pulseFourierSeries
  input Real u;
  input Real frequency;
  input Real duty;
  input Integer n=25;
  input Real t;
  output Real y;
algorithm
  y:= duty + (2/Modelica.Constants.pi) * sum( (1 ./linspace(1,n,n)) .* Modelica.Math.sin(Modelica.Constants.pi*linspace(1,n,n)*duty) .* Modelica.Math.cos(frequency*linspace(1,n,n)*t));

end pulseFourierSeries;

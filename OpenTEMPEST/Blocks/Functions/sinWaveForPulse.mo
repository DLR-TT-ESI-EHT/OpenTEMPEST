within OpenTEMPEST.Blocks.Functions;
function sinWaveForPulse
  input Real u;
  input Real duration;
  input Real duty;
  output Real y;

algorithm
  y:= 0.5*(cos(2*Modelica.Constants.pi/duration) - cos(2*Modelica.Constants.pi*u/duration)) - (1-duty);

end sinWaveForPulse;

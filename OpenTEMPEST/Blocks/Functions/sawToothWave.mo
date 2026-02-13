within OpenTEMPEST.Blocks.Functions;
function sawToothWave
  input Real u;
  input Real duration;
  input Real duty;
  output Real y;
algorithm
  y:= u/duration - floor(0.5+u/duration) - (0.5-duty);

end sawToothWave;

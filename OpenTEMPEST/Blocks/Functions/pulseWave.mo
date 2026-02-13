within OpenTEMPEST.Blocks.Functions;
function pulseWave
  input Real u;
  input Real duration;
  input Real duty;
  output Real y;

algorithm

  y := smoothHeaviside(Functions.sawToothWave(
    u,
    duration,
    duty));

end pulseWave;

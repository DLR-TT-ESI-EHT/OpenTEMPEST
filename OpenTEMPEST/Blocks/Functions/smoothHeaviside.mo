within OpenTEMPEST.Blocks.Functions;
function smoothHeaviside
  input Real u;
  output Real y;
algorithm

  y:=0.5 + 0.5*Modelica.Math.tanh(10000*u);

//   if u>=0 then
//     y:=1;
//   else
//     y:=0;
//   end if;

end smoothHeaviside;

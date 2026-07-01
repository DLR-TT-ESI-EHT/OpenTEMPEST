within OpenTEMPEST.Blocks.Functions;
function polySolid
  "
  Function to calculate the thermophysical properties of a solid, for a given set of coefficients
  "

  input Real Tref;
  input Real c[4];
  output Real y;

algorithm

  y :=c[1]*Tref.^3 .+ c[2]*Tref.^2 .+ c[3]*Tref .+ c[4];

end polySolid;

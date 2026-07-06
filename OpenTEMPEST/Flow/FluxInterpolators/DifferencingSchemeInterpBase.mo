within OpenTEMPEST.Flow.FluxInterpolators;
partial function DifferencingSchemeInterpBase
  input Integer N "Number of points in U";
  input Real U[N] "Field variable values";
  input Real U0 "Ghost cell value (upstream)";
  output Real Y[N] "Upwind differenced values";

end DifferencingSchemeInterpBase;

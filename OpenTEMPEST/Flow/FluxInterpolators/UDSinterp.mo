within OpenTEMPEST.Flow.FluxInterpolators;
function UDSinterp
  "UDS: Upwind Differencing Scheme (1st Order)"
  extends DifferencingSchemeInterpBase;
algorithm
  // Guard: Ensure we have at least two elements to apply the scheme
  assert(N >= 2, "N must be at least 2 for LUDS to be applied");
  for i in 1:N loop
    Y[i] := U[i];
  end for;

end UDSinterp;

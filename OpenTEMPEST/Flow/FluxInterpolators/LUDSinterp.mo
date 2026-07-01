within OpenTEMPEST.Flow.FluxInterpolators;
function LUDSinterp
  "LUDS: Linear Upwind Differencing Scheme (2nd Order)"
  extends DifferencingSchemeInterpBase;
algorithm
  // Guard: Ensure we have at least two elements to apply the scheme
  assert(N >= 2, "N must be at least 2 for LUDS to be applied");

  // First internal point (uses ghost cell U0)
  Y[1] := (3*U[1] - U0) / 2;

  // Main loop: linear upwind differencing for the rest
  for i in 2:N loop
    Y[i] := (3*U[i] - U[i-1]) / 2;
  end for;

end LUDSinterp;

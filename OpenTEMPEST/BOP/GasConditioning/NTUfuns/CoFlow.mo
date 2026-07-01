within OpenTEMPEST.BOP.GasConditioning.NTUfuns;
function CoFlow
  extends NTUbase;

algorithm
  eta:=max(1e-5,(1 - exp(-NTU*(1+Cr)))/(1+Cr));

end CoFlow;

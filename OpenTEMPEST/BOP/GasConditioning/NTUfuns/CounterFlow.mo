within OpenTEMPEST.BOP.GasConditioning.NTUfuns;
function CounterFlow
  extends NTUbase;

algorithm
  eta:=max(1e-5,if abs(Cr-1)>1e-2 then (1 - exp(-NTU*(1-Cr)))/(1-Cr*exp(-NTU*(1-Cr))) else NTU/(1+NTU));

end CounterFlow;

within OpenTEMPEST.BOP.GasConditioning.HexEtaTypes;
model NTU_FixedUA
  extends EtaBase;

  parameter Real UA=1;
  Real NTU = UA/(C_min + 1e-12);
  replaceable function NTUmethod =
      NTUfuns.CounterFlow                                                                     constrainedby
    NTUfuns.NTUbase
  annotation(choicesAllMatching=True);

equation

  eta = NTUmethod(NTU, C_r);

end NTU_FixedUA;

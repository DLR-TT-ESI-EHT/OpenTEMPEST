within OpenTEMPEST.BOP.GasConditioning.HexEtaTypes;
model ConstantEta
  extends EtaBase;
  parameter Real etaNom = 0.99;

equation
  eta = etaNom;

end ConstantEta;

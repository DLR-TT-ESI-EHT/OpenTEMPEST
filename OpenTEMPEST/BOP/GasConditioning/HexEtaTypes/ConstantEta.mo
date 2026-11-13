within OpenTEMPEST.BOP.GasConditioning.HexEtaTypes;
model ConstantEta "Constant heat exchanger effectiveness"
  extends EtaBase;
  parameter Real etaNom = 0.99;

equation
  eta = etaNom;

end ConstantEta;

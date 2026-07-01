within OpenTEMPEST.BOP.GasConditioning.HexEtaTypes;
model NTU_VarUA
  extends EtaBase;

  parameter Real b=1;

  parameter Real mf_NomA=1;
  parameter Real mf_NomB=1;
  parameter Real hA_RefA = 1;
  parameter Real hA_RefB = 1;

  parameter Real TnomA=1000;
  parameter Real TnomB=500;
  parameter Real PnomA=1e5;
  parameter Real PnomB=1e5;
  parameter Real XnomA[MediumA.nXi] = MediumA.reference_X;
  parameter Real XnomB[MediumB.nXi] = MediumB.reference_X;

  Real lambda_RefA = MediumA.thermalConductivity(MediumA.setState_pTX(PnomA, TnomA, XnomA));
  Real lambda_RefB = MediumB.thermalConductivity(MediumB.setState_pTX(PnomB, TnomB, XnomB));
  Real cp_RefA = MediumA.specificHeatCapacityCp(MediumA.setState_pTX(PnomA, TnomA, XnomA));
  Real cp_RefB = MediumB.specificHeatCapacityCp(MediumB.setState_pTX(PnomB, TnomB, XnomB));
  Real mu_RefA = MediumA.dynamicViscosity(MediumA.setState_pTX(PnomA, TnomA, XnomA));
  Real mu_RefB = MediumB.dynamicViscosity(MediumB.setState_pTX(PnomB, TnomB, XnomB));

  Real UA;
  Real NTU = UA/(C_min + 1e-12);
  replaceable function NTUmethod =
      NTUfuns.CounterFlow                                                                     constrainedby
    NTUfuns.NTUbase
  annotation(choicesAllMatching=True);

equation

  1/UA = 1/(hA_RefA*(mf_A/mf_NomA)^b*(mu_A/mu_RefA)^(1/3-b)*(cp_RefA/cp_A)^(1/3)*(lambda_A/lambda_RefA)^(2/3) + 1e-9) +
         1/(hA_RefB*(mf_B/mf_NomB)^b*(mu_B/mu_RefB)^(1/3-b)*(cp_RefB/cp_B)^(1/3)*(lambda_B/lambda_RefB)^(2/3) + 1e-9);
  eta = NTUmethod(NTU, C_r);

end NTU_VarUA;

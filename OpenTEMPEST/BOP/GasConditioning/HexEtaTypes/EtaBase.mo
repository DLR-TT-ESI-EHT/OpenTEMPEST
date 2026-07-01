within OpenTEMPEST.BOP.GasConditioning.HexEtaTypes;
partial model EtaBase

  Real eta;
  input Real C_min;
  input Real C_max;
  input Real mf_A;
  input Real mf_B;
  input Real lambda_A;
  input Real lambda_B;
  input Real cp_A;
  input Real cp_B;
  input Real mu_A;
  input Real mu_B;

  Real C_r = C_min/(C_max+1e-12);

  replaceable package MediumA = Modelica.Media.Interfaces.PartialMedium annotation(Dialog(enable=true, showStartAttribute = false));
  replaceable package MediumB = Modelica.Media.Interfaces.PartialMedium annotation(Dialog(enable=false, showStartAttribute = false));

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end EtaBase;

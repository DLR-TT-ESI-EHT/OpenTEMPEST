within OpenTEMPEST.SOC.Electrochem.Components;
model BV_Steam "Steam Electrolysis with B-V, Ohmic ASR and Diffusion"
  extends OpenTEMPEST.SOC.Electrochem.Components.BV_ElectrochemBase;

import SI = Modelica.SIunits;

  SI.MolarEnthalpy H_h2o, H_o2, H_h2;
  SI.CurrentDensity Jo_FE;
  SI.MolarEnthalpy dHr;

  replaceable Real E_fe = 95160;
  replaceable Real a = 0.04;
  replaceable Real b = 0.18;
  replaceable Real gammaFE = 1.52e5;

equation
  J_H = J;
  J_C = 0;

  // Nernst
  // Gibbs energy Calculations
  delGr = go_T + Modelica.Constants.R*Tpen*(Modelica.Math.log(abs(yF[5])) - Modelica.Math.log(abs(yF[1])) - 0.5*Modelica.Math.log(abs(yA[1])) - 0.5*Modelica.Math.log(P_A/p0));
  go_T = (0.05354*(Tpen) - 245.9767)*1e3;  // linear relationship from NASA polynomials

  // Activation Anode
  //Butler-Volmer:
  J = Jo_FE*(Modelica.Math.exp(((alphaFE*2*Modelica.Constants.F)/(Modelica.Constants.R*Tpen))*UactFE) - Modelica.Math.exp(((-(1 - alphaFE)*2*Modelica.Constants.F)/(Modelica.Constants.R*Tpen))*UactFE));
  Jo_FE = gammaFE*Tpen*((P_F/p0*abs(yF_tpb[1]))^a)*((P_F/p0*abs(yF_tpb[5]))^b)*Modelica.Math.exp(-E_fe/(Modelica.Constants.R*Tpen));

  // q_electroChem
  H_h2o =Medium.Fuel_CH4.MMX[5]*
    Modelica.Media.IdealGases.Common.Functions.h_T(
    T=Tpen,
    exclEnthForm=false,
    refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
    data=Modelica.Media.IdealGases.Common.SingleGasesData.H2O);
  H_o2 =Medium.Air_Medium.MMX[1]*
    Modelica.Media.IdealGases.Common.Functions.h_T(
    T=Tpen,
    exclEnthForm=false,
    refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
    data=Modelica.Media.IdealGases.Common.SingleGasesData.O2);
  H_h2 =Medium.Fuel_CH4.MMX[1]*
    Modelica.Media.IdealGases.Common.Functions.h_T(
    T=Tpen,
    exclEnthForm=false,
    refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
    data=Modelica.Media.IdealGases.Common.SingleGasesData.H2);
  dHr = H_h2o - 0.5*H_o2 - H_h2;
  q_electroChem = -r*dHr;

  //Diffusion overpotential
  UDiffFE = -Modelica.Constants.R*Tpen/(2*Modelica.Constants.F)*(Modelica.Math.log(abs(yF_tpb[1])) + Modelica.Math.log(abs(yF[5])) - Modelica.Math.log(abs(yF_tpb[5])) - Modelica.Math.log(abs(yF[1])));

  if diffusionActive then
    yF_tpb[1] = yF[1] - J*Modelica.Constants.R*Tpen / (2*Modelica.Constants.F*P_F) * (1/Deff_h2_k + 1/Deff_h2h2o)*tauFE - (Deff_h2h2o - Deff_h2n2) * Deff_h2on2 / ((Deff_h2on2 - Deff_h2n2) * Deff_h2h2o) * yF[6] * (Modelica.Math.exp(J * Modelica.Constants.R * Tpen / (2 * Modelica.Constants.F * P_F) * (1/Deff_h2n2 - 1/Deff_h2on2) * tauFE) - 1) - (Deff_h2h2o - Deff_h2co) * Deff_h2oco / ((Deff_h2oco - Deff_h2co) * Deff_h2h2o) * yF[4] * (Modelica.Math.exp(J * Modelica.Constants.R * Tpen / (2 * Modelica.Constants.F * P_F) * (1 / Deff_h2co - 1 / Deff_h2oco) * tauFE) - 1) - (Deff_h2h2o - Deff_h2co2) * Deff_h2oco2 / ((Deff_h2oco2 - Deff_h2co2) * Deff_h2h2o) * yF[3] * (Modelica.Math.exp(J * Modelica.Constants.R * Tpen / (2 * Modelica.Constants.F * P_F) * (1 / Deff_h2co2 - 1 / Deff_h2oco2) * tauFE) - 1) - (Deff_h2h2o - Deff_h2ch4) * Deff_h2och4 / ((Deff_h2och4 - Deff_h2ch4) * Deff_h2h2o) * yF[2] * (Modelica.Math.exp(J * Modelica.Constants.R * Tpen / (2 * Modelica.Constants.F * P_F) * (1 / Deff_h2ch4 - 1 / Deff_h2och4) * tauFE) - 1);
    yF_tpb[5] = yF[5] + J*Modelica.Constants.R*Tpen / (2*Modelica.Constants.F*P_F) * (1/Deff_h2o_k + 1/Deff_h2h2o)*tauFE + (Deff_h2h2o - Deff_h2on2) * Deff_h2n2 / ((Deff_h2on2 - Deff_h2n2) * Deff_h2h2o) * yF[6] * (Modelica.Math.exp(J * Modelica.Constants.R * Tpen / (2 * Modelica.Constants.F * P_F) * (1/Deff_h2n2 - 1/Deff_h2on2) * tauFE) - 1) + (Deff_h2h2o - Deff_h2oco) * Deff_h2co / ((Deff_h2oco - Deff_h2co) * Deff_h2h2o) * yF[4] * (Modelica.Math.exp(J * Modelica.Constants.R * Tpen / (2 * Modelica.Constants.F * P_F) * (1 / Deff_h2co - 1 / Deff_h2oco) * tauFE) - 1) + (Deff_h2h2o - Deff_h2oco2) * Deff_h2co2 / ((Deff_h2oco2 - Deff_h2co2) * Deff_h2h2o) * yF[3] * (Modelica.Math.exp(J * Modelica.Constants.R * Tpen / (2 * Modelica.Constants.F * P_F) * (1 / Deff_h2co2 - 1 / Deff_h2oco2) * tauFE) - 1) + (Deff_h2h2o - Deff_h2och4) * Deff_h2ch4 / ((Deff_h2och4 - Deff_h2ch4) * Deff_h2h2o) * yF[2] * (Modelica.Math.exp(J * Modelica.Constants.R * Tpen / (2 * Modelica.Constants.F * P_F) * (1 / Deff_h2ch4 - 1 / Deff_h2och4) * tauFE) - 1);
    yF_tpb[4] = yF[4]*Modelica.Math.exp( (Modelica.Constants.R*Tpen*tauFE/(2*Modelica.Constants.F*P_F))*(J*((1/Deff_h2co)-(1/Deff_h2oco))));
    yF_tpb[3] = yF[3]*Modelica.Math.exp( (Modelica.Constants.R*Tpen*tauFE/(2*Modelica.Constants.F*P_F))*(J*((1/Deff_h2co2)-(1/Deff_h2oco2))));
    yF_tpb[6] = yF[6]*Modelica.Math.exp( (Modelica.Constants.R*Tpen*tauFE/(2*Modelica.Constants.F*P_F))*(J*((1/Deff_h2n2)-(1/Deff_h2on2))));
    yF_tpb[2] = yF[2]*Modelica.Math.exp( (Modelica.Constants.R*Tpen*tauFE/(2*Modelica.Constants.F*P_F))*(J*((1/Deff_h2ch4)-(1/Deff_h2och4))));
  else
    yF_tpb = yF;
  end if;

end BV_Steam;

within OpenTEMPEST.SOC.Electrochem.Components;
model BV_Co "Co-Electrolysis with B-V, Ohmic ASR and Diffusion"
  extends OpenTEMPEST.SOC.Electrochem.Components.BV_ElectrochemBase;

import SI = Modelica.SIunits;

  replaceable Real E_fe = 95160;
  replaceable Real E_fe_C = 1.12*111700;
  replaceable Real a = 0.04;
  replaceable Real b = 0.18;
  replaceable Real a_C = 0.04;
  replaceable Real b_C = 0.18;
  replaceable Real gammaFE = 1.52e5;
  replaceable Real gammaFE_C = 6.63e5;

  SI.MolarEnthalpy H_o2, H_co2, H_co, H_h2o, H_h2;
  SI.CurrentDensity Jo_FEC, Jo_FEH;
  SI.MolarEnthalpy dHr_H, dHr_C;

  Units.MolarFlux r_H=J_H/(2*Modelica.Constants.F) "reaction rate of H2/H2O electrochemical reactions mol/m^2 S";
  Units.MolarFlux r_C=J_C/(2*Modelica.Constants.F) "reaction rate of CO/CO2 electrochemical reactions mol/m^2 S";

equation
  // Nernst
    // Gibbs energy Calculations
    delGr = 0.5*(go_T + Modelica.Constants.R*Tpen*(Modelica.Math.log(abs(yF[3])) + Modelica.Math.log(abs(yF[5])) - Modelica.Math.log(abs(yF[1])) - Modelica.Math.log(abs(yF[4])) - Modelica.Math.log(abs(yA[1])) - Modelica.Math.log(abs(P_A/p0))));
    go_T = (0.141*Tpen - 529.045)*1e3; // linear relationship from NASA polynomials

  // Activation Anode
    //UactFE = Modelica.Constants.R*Tpen/Modelica.Constants.F*Modelica.Math.asinh(J_C/(2*Jo_FEC));
    J_C = Jo_FEC*(Modelica.Math.exp(((alphaFE*2*Modelica.Constants.F)/(Modelica.Constants.R*Tpen))*UactFE) - Modelica.Math.exp(((-(1-alphaFE)*2*Modelica.Constants.F)/(Modelica.Constants.R*Tpen))*UactFE));
    Jo_FEC = gammaFE_C*Tpen*((P_F/p0*abs(yF_tpb[4]))^a_C)*((P_F/p0*abs(yF_tpb[3]))^b_C)*Modelica.Math.exp(-E_fe_C/(Modelica.Constants.R*Tpen));
    Jo_FEH = gammaFE*Tpen*((P_F/p0*abs(yF_tpb[1]))^a)*((P_F/p0*abs(yF_tpb[5]))^b)*Modelica.Math.exp(-E_fe/(Modelica.Constants.R*Tpen));
    J_H = Jo_FEH*(Modelica.Math.exp(((alphaFE*2*Modelica.Constants.F)/(Modelica.Constants.R*Tpen))*UactFE) - Modelica.Math.exp(((-(1-alphaFE)*2*Modelica.Constants.F)/(Modelica.Constants.R*Tpen))*UactFE));

  // q_electroChem
    H_h2o =Medium.Fuel_CH4.MMX[5]*
    Modelica.Media.IdealGases.Common.Functions.h_T(
    T=Tpen,
    exclEnthForm=false,
    refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
    data=Modelica.Media.IdealGases.Common.SingleGasesData.H2O);
    H_o2  =Medium.Air_Medium.MMX[1]*
    Modelica.Media.IdealGases.Common.Functions.h_T(
    T=Tpen,
    exclEnthForm=false,
    refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
    data=Modelica.Media.IdealGases.Common.SingleGasesData.O2);
    H_h2  =Medium.Fuel_CH4.MMX[1]*
    Modelica.Media.IdealGases.Common.Functions.h_T(
    T=Tpen,
    exclEnthForm=false,
    refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
    data=Modelica.Media.IdealGases.Common.SingleGasesData.H2);
    H_co2 =Medium.Fuel_CH4.MMX[3]*
    Modelica.Media.IdealGases.Common.Functions.h_T(
    T=Tpen,
    exclEnthForm=false,
    refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
    data=Modelica.Media.IdealGases.Common.SingleGasesData.CO2);
    H_co  =Medium.Fuel_CH4.MMX[4]*
    Modelica.Media.IdealGases.Common.Functions.h_T(
    T=Tpen,
    exclEnthForm=false,
    refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
    data=Modelica.Media.IdealGases.Common.SingleGasesData.CO);
    dHr_C = H_co2 - 0.5 * H_o2 - H_co;
    dHr_H = H_h2o - 0.5 * H_o2 - H_h2;
    //dHr = dHr_C + dHr_H;
    J = J_H + J_C;
    q_electroChem = -r_H*dHr_H - r_C*dHr_C;

  //Diffusion overpotential
    UDiffFE = -Modelica.Constants.R*Tpen/(4*Modelica.Constants.F)*(Modelica.Math.log(abs(yF[5])) + Modelica.Math.log(abs(yF[3])) + Modelica.Math.log(abs(yF_tpb[1])) + Modelica.Math.log(abs(yF_tpb[4])) - (Modelica.Math.log(abs(yF_tpb[5])) + Modelica.Math.log(abs(yF_tpb[3])) + Modelica.Math.log(abs(yF[1])) + Modelica.Math.log(abs(yF[4]))));

  if diffusionActive then
    yF_tpb[1]  = (yF[1] - J_H*( (Modelica.Constants.R*Tpen*tauFE)/(2*Modelica.Constants.F*P_F)) *( ((1/Deff_h2_k) + (1/Deff_h2h2o)) + ((1/Deff_h2ch4) - (1/Deff_h2h2o))*yF_tpb[2] + ((1/Deff_h2co2) - (1/Deff_h2h2o))*yF_tpb[3] + ((1/Deff_h2co) - (1/Deff_h2h2o))*yF_tpb[4] + ((1/Deff_h2n2) - (1/Deff_h2h2o))*yF_tpb[6]))/(1 + J_C*( (Modelica.Constants.R*Tpen*tauFE)/(2*Modelica.Constants.F*P_F))*((1/Deff_h2co2) - (1/Deff_h2co)));
    yF_tpb[5]  = (yF[5] + J_H*( (Modelica.Constants.R*Tpen*tauFE)/(2*Modelica.Constants.F*P_F)) *( ((1/Deff_h2o_k) + (1/Deff_h2h2o)) + ((1/Deff_h2och4) - (1/Deff_h2h2o))*yF_tpb[2] + ((1/Deff_h2oco2) - (1/Deff_h2h2o))*yF_tpb[3] + ((1/Deff_h2oco) - (1/Deff_h2h2o))*yF_tpb[4] + ((1/Deff_h2on2) - (1/Deff_h2h2o))*yF_tpb[6]))/(1 + J_C*( (Modelica.Constants.R*Tpen*tauFE)/(2*Modelica.Constants.F*P_F))*((1/Deff_h2oco2) - (1/Deff_h2oco)));
    yF_tpb[4]  = (yF[4] - J_C*( (Modelica.Constants.R*Tpen*tauFE)/(2*Modelica.Constants.F*P_F)) *( ((1/Deff_co_k) + (1/Deff_coco2)) + ((1/Deff_h2co) - (1/Deff_coco2))*yF_tpb[1] + ((1/Deff_coch4) - (1/Deff_coco2))*yF_tpb[2] + ((1/Deff_h2oco) - (1/Deff_coco2))*yF_tpb[5] + ((1/Deff_con2) - (1/Deff_coco2))*yF_tpb[6]))/(1 + J_H*( (Modelica.Constants.R*Tpen*tauFE)/(2*Modelica.Constants.F*P_F))*((1/Deff_h2oco) - (1/Deff_h2co)));
    yF_tpb[3]  = (yF[3] + J_C*( (Modelica.Constants.R*Tpen*tauFE)/(2*Modelica.Constants.F*P_F)) *( ((1/Deff_co2_k) + (1/Deff_coco2)) + ((1/Deff_h2co2) - (1/Deff_coco2))*yF_tpb[1] + ((1/Deff_co2ch4) - (1/Deff_coco2))*yF_tpb[2] + ((1/Deff_h2oco2) - (1/Deff_coco2))*yF_tpb[5] + ((1/Deff_co2n2) - (1/Deff_coco2))*yF_tpb[6]))/(1 + J_H*( (Modelica.Constants.R*Tpen*tauFE)/(2*Modelica.Constants.F*P_F))*((1/Deff_h2oco2) - (1/Deff_h2co2)));
    yF_tpb[6] = yF[6]*Modelica.Math.exp( (Modelica.Constants.R*Tpen*tauFE/(2*Modelica.Constants.F*P_F))*(J_H*((1/Deff_h2n2)-(1/Deff_h2on2)) + J_C*((1/Deff_con2)-(1/Deff_co2n2))));
    yF_tpb[2] = yF[2]*Modelica.Math.exp((Modelica.Constants.R*Tpen*tauFE/(2*Modelica.Constants.F*P_F))*(J_H*((1/Deff_h2ch4)-(1/Deff_h2och4)) + J_C*((1/Deff_coch4)-(1/Deff_co2ch4))));
  else
    yF_tpb = yF;
  end if;
end BV_Co;

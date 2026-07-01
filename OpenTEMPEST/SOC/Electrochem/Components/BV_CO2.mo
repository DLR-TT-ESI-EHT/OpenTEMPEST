within OpenTEMPEST.SOC.Electrochem.Components;
model BV_CO2 "CO2 Electrolysis with B-V, Ohmic ASR and Diffusion"
   extends OpenTEMPEST.SOC.Electrochem.Components.BV_ElectrochemBase;

import SI = Modelica.SIunits;

  replaceable Real E_fe_C = 1.12*111700 "Activation energy fuel electrode for CO/CO2";
  replaceable Real a_C = 0.04;
  replaceable Real b_C = 0.18;
  replaceable Real gammaFE_C = 6.63e5 "Pre-exponential factor fuel electrode CO/CO2 Jo";

  SI.MolarEnthalpy H_o2, H_co2, H_co;
  SI.CurrentDensity Jo_FE "Fuel electrode exchange current density";
  SI.MolarEnthalpy dHr "Enthalpy change of reaction";

equation
  J_C = J;
  J_H = 0;

  // Nernst
    // Gibbs energy Calculations
    delGr = go_T + Modelica.Constants.R*Tpen*(Modelica.Math.log(abs(yF[3])) - Modelica.Math.log(abs(yF[4])) - 0.5*Modelica.Math.log(abs(yA[1])) - 0.5*Modelica.Math.log(abs(P_A/p0)));
    go_T = (0.0875*(Tpen) - 283.0683)*1e3; // linear relationship from NASA polynomials

  // Activation Anode
    J = Jo_FE*(Modelica.Math.exp(((alphaFE*2*Modelica.Constants.F)/(Modelica.Constants.R*Tpen))*UactFE) - Modelica.Math.exp(((-(1 - alphaFE)*2*Modelica.Constants.F)/(Modelica.Constants.R*Tpen))*UactFE));
    //UactFE = Modelica.Constants.R*Tpen/Modelica.Constants.F*Modelica.Math.asinh(J/(2*Jo_FE));
    Jo_FE = gammaFE_C*Tpen*((P_F/p0*yF_tpb[4])^a_C)*((P_F/p0*yF_tpb[3])^b_C)*Modelica.Math.exp(-E_fe_C/(Modelica.Constants.R*Tpen));

  // q_electroChem
    H_co2 =Medium.Fuel_CH4.MMX[3]*
    Modelica.Media.IdealGases.Common.Functions.h_T(
    T=Tpen,
    exclEnthForm=false,
    refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
    data=Modelica.Media.IdealGases.Common.SingleGasesData.CO2);
    H_o2  =Medium.Air_Medium.MMX[1]*
    Modelica.Media.IdealGases.Common.Functions.h_T(
    T=Tpen,
    exclEnthForm=false,
    refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
    data=Modelica.Media.IdealGases.Common.SingleGasesData.O2);
    H_co  =Medium.Fuel_CH4.MMX[4]*
    Modelica.Media.IdealGases.Common.Functions.h_T(
    T=Tpen,
    exclEnthForm=false,
    refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
    data=Modelica.Media.IdealGases.Common.SingleGasesData.CO);
    dHr = H_co2 - 0.5 * H_o2 - H_co;
    q_electroChem = -r*dHr;

  //Diffusion overpotential
    UDiffFE = -Modelica.Constants.R*Tpen/(2*Modelica.Constants.F)*(Modelica.Math.log(abs(yF_tpb[4])) + Modelica.Math.log(abs(yF_tpb[3])) - Modelica.Math.log(abs(yF_tpb[3])) - Modelica.Math.log(abs(yF_tpb[4])));

  if diffusionActive then
    // verify if analytical solution here is correct
      yF_tpb[1] = yF[1]*Modelica.Math.exp( (Modelica.Constants.R*Tpen*tauFE/(2*Modelica.Constants.F*P_F))*(J*((1/Deff_h2co)-(1/Deff_h2co2))));
      yF_tpb[5] = yF[5]*Modelica.Math.exp( (Modelica.Constants.R*Tpen*tauFE/(2*Modelica.Constants.F*P_F))*(J*((1/Deff_h2oco)-(1/Deff_h2oco2))));
      yF_tpb[4] = (yF[4] - J*( (Modelica.Constants.R*Tpen*tauFE)/(2*Modelica.Constants.F*P_F)) *( ((1/Deff_co_k) + (1/Deff_coco2)) + ((1/Deff_h2co) - (1/Deff_coco2))*yF_tpb[1] + ((1/Deff_coch4) - (1/Deff_coco2))*yF_tpb[2] + ((1/Deff_h2oco) - (1/Deff_coco2))*yF_tpb[5] + ((1/Deff_con2) - (1/Deff_coco2))*yF_tpb[6]));
      yF_tpb[3] = (yF[3] + J*( (Modelica.Constants.R*Tpen*tauFE)/(2*Modelica.Constants.F*P_F)) *( ((1/Deff_co2_k) + (1/Deff_coco2)) + ((1/Deff_h2co2) - (1/Deff_coco2))*yF_tpb[1] + ((1/Deff_co2ch4) - (1/Deff_coco2))*yF_tpb[2] + ((1/Deff_h2oco2) - (1/Deff_coco2))*yF_tpb[5] + ((1/Deff_co2n2) - (1/Deff_coco2))*yF_tpb[6]));
      yF_tpb[6] = yF[6]*Modelica.Math.exp( (Modelica.Constants.R*Tpen*tauFE/(2*Modelica.Constants.F*P_F))*(J*((1/Deff_con2)-(1/Deff_co2n2))));
      yF_tpb[2] = yF[2]*Modelica.Math.exp((Modelica.Constants.R*Tpen*tauFE/(2*Modelica.Constants.F*P_F))*(J*((1/Deff_coch4)-(1/Deff_co2ch4))));
  else
    yF_tpb = yF;
  end if;
  annotation (Documentation(info="<html>
  <body>
    <p>References for parameters:</p>
    <ul>
      <li>DOI: 10.1149/1945-7111/ad5e01</li>
    </ul>
  </body>
</html>"));
end BV_CO2;

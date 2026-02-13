within OpenTEMPEST.SOC.Electrochem.Components;
model ASR_Steam
  extends OpenTEMPEST.SOC.Electrochem.Components.ASR_ElectrochemBase(redeclare replaceable model
                        ASRobj = OpenTEMPEST.SOC.Electrochem.ASR.ASR_log (A=
            6464.7, B=6.8416));

    import SI = Modelica.SIunits;
    SI.MolarEnthalpy H_h2o, H_o2, H_h2;

    Real go_T;

equation
  J_H = J;
  J_C = 0;
  H_h2o =Medium.Fuel_CH4.MMX[5]*Modelica.Media.IdealGases.Common.Functions.h_T(
    T=Tpen,
    exclEnthForm=false,
    refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
    data=Modelica.Media.IdealGases.Common.SingleGasesData.H2O);
  H_o2 =Medium.Air_Medium.MMX[1]*Modelica.Media.IdealGases.Common.Functions.h_T(
    T=Tpen,
    exclEnthForm=false,
    refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
    data=Modelica.Media.IdealGases.Common.SingleGasesData.O2);
  H_h2 =Medium.Fuel_CH4.MMX[1]*Modelica.Media.IdealGases.Common.Functions.h_T(
    T=Tpen,
    exclEnthForm=false,
    refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
    data=Modelica.Media.IdealGases.Common.SingleGasesData.H2);
  dHr = H_h2o - 0.5*H_o2 - H_h2;

  // Gibbs energy Calculations
  go_T = (0.05354*(Tpen) - 245.9767)*1e3;
  // linear relationship from NASA polynomials
  delGr = go_T + Modelica.Constants.R*Tpen*(Modelica.Math.log(abs(yF[5])) - Modelica.Math.log(abs(yF[1])) - 0.5*Modelica.Math.log(abs(yA[1])) - 0.5*Modelica.Math.log(P_A/p0));
  //J/mol
  r = J/2/Modelica.Constants.F;

end ASR_Steam;

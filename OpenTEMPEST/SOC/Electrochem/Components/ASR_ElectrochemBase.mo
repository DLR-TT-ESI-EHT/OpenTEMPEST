within OpenTEMPEST.SOC.Electrochem.Components;
partial model ASR_ElectrochemBase "Base model for electrochemistry with only ASR approach"
  extends OpenTEMPEST.SOC.Electrochem.Components.ElectrochemBase;

  import SI = Modelica.SIunits;

  replaceable model ASRobj = OpenTEMPEST.SOC.Electrochem.ASR.ASR_Base
                                            constrainedby
    OpenTEMPEST.SOC.Electrochem.ASR.ASR_Base                                                                         annotation(choicesAllMatching=true);

  ASRobj asrobj(Tpen=Tpen);
  SI.Voltage Uideal;
  SI.SpecificGibbsFreeEnergy delGr "Specific Gibbs free energy change of reaction";
  SI.MolarEnthalpy dHr "Molar enthalpy of reaction";
equation
  // Electrochemical Heat Production
  q_electroChem = -r.*dHr;

  // Voltage calculation
  Uideal = (-delGr/(2*Modelica.Constants.F));
  Uop = Uideal - J*asrobj.ASR*1e-4;

end ASR_ElectrochemBase;

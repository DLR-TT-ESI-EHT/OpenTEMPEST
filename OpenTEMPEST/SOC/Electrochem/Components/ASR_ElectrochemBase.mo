within OpenTEMPEST.SOC.Electrochem.Components;
partial model ASR_ElectrochemBase
  extends OpenTEMPEST.SOC.Electrochem.Components.ElectrochemBase;

  import SI = Modelica.SIunits;

  replaceable model ASRobj = OpenTEMPEST.SOC.Electrochem.ASR.ASR_Base
                                            constrainedby
    OpenTEMPEST.SOC.Electrochem.ASR.ASR_Base                                                                         annotation(choicesAllMatching=true);

  ASRobj asrobj(Tpen=Tpen);
  SI.Voltage Uideal;
  Real delGr;
  Real dHr;
equation
  // Electrochemical Heat Production
  q_electroChem = -r.*dHr;

  // Voltage calculation
  Uideal = (-delGr/(2*Modelica.Constants.F));
  Uop = Uideal - J*asrobj.ASR*1e-4;

end ASR_ElectrochemBase;

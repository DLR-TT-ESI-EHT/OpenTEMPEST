within OpenTEMPEST.SOC.Electrochem.ASR;
model ASR_Linear "Linear ASR relationship: ASR = A*Tpen - B"
  extends ASR_Base;

equation

    ASR = A*Tpen - B;

end ASR_Linear;

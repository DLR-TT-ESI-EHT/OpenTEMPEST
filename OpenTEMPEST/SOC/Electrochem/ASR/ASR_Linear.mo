within OpenTEMPEST.SOC.Electrochem.ASR;
model ASR_Linear "Linear ASR relationship: ASR /(ohm cm^2) = A*Tpen - B"
  extends ASR_Base;

equation

    ASR*1e4 = A*Tpen - B;

end ASR_Linear;

within OpenTEMPEST.SOC.Electrochem.ASR;
model ASR_Linear "Linear ASR correlation, ASR/(Ohm.cm2) = A*T - B"
  extends ASR_Base;

equation

    ASR*1e4 = A*Tpen - B;

end ASR_Linear;

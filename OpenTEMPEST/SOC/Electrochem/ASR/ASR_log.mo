within OpenTEMPEST.SOC.Electrochem.ASR;
model ASR_log "Logarithmic ASR correlation, ln(ASR/(Ohm.cm2)) = A/T - B"
  extends ASR_Base;

equation

  Modelica.Math.log(ASR*1e4) = (A/Tpen) - B;

end ASR_log;

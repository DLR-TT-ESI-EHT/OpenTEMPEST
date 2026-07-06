within OpenTEMPEST.SOC.Electrochem.ASR;
model ASR_log "Logarithmic ASR relationship: ln(ASR/(ohm cm^2)) = A/Tpen - B"
  extends ASR_Base;

equation

  Modelica.Math.log(ASR) = (A/Tpen) - B;

end ASR_log;

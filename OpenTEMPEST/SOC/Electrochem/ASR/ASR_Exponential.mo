within OpenTEMPEST.SOC.Electrochem.ASR;
model ASR_Exponential
  "Exponential ASR relation: ASR/(ohm cm^2) = A*exp(B*(Tpen-273.15))"
  extends ASR_Base;

equation

    ASR*1e4 = A*Modelica.Math.exp(B*(Tpen-273.15));

end ASR_Exponential;

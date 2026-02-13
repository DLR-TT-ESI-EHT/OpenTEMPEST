within OpenTEMPEST.SOC.Electrochem.ASR;
model ASR_Exponential
  "Exponential ASR relation: ASR/(ohm cm^2) = A*exp(B*(Tpen-273.15))"
  extends ASR_Base;

equation

    ASR = A*Modelica.Math.exp(B*(Tpen-273.15));

end ASR_Exponential;

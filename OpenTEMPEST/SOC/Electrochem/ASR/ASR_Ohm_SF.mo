within OpenTEMPEST.SOC.Electrochem.ASR;
model ASR_Ohm_SF "Ohmic ASR for SF - No dimensions"

  extends OpenTEMPEST.SOC.Electrochem.ASR.ASR_OhmBase;

  parameter Real Yo = 27.266e-2 "[Ohm*cm^2]";

equation

    ASR = (Yo +A*exp(B*Tpen))/100^2;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end ASR_Ohm_SF;

within OpenTEMPEST.SOC.Electrochem.ASR;
model ASR_Ohm_ESC "Ohmic ASR for ESC - No dimensions"

  extends OpenTEMPEST.SOC.Electrochem.ASR.ASR_OhmBase;

  parameter Real Yo = 27.266e-2 "[Ohm*cm^2]";

equation

    ASR = (Yo +A*exp(B*Tpen))/100^2;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    Documentation(revisions="<html>
</html>", info="<html>
<html><body>
DOI: 10.1016/j.ijhydene.2018.12.168
</html>"));
end ASR_Ohm_ESC;

within OpenTEMPEST.SOC.Electrochem.ASR;
model ASR_Ohm_ESC "Ohmic ASR of an electrolyte supported cell (ESC)"

  extends OpenTEMPEST.SOC.Electrochem.ASR.ASR_OhmBase;

  parameter Real Yo = 27.266e-2 "[Ohm*m^2]";

equation

    ASR*1e4 = Yo +A*exp(B*Tpen);

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    Documentation(revisions="<html>
</html>", info="<html>
<p>Ohmic ASR as a function of the PEN temperature. Parameters from <a href=\"https://doi.org/10.1016/j.ijhydene.2018.12.168\">doi:10.1016/j.ijhydene.2018.12.168</a>.</p>
</html>"));
end ASR_Ohm_ESC;

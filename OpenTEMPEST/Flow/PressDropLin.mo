within OpenTEMPEST.Flow;
model PressDropLin "Linear pressure drop with pressure difference output"
  extends ThermoPower.Gas.PressDropLin;

  Real dp = inlet.p - outlet.p;

end PressDropLin;

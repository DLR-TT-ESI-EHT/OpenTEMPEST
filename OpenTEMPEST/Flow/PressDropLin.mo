within OpenTEMPEST.Flow;
model PressDropLin "added dp calculation in thermopower pressure drop lin"
  extends ThermoPower.Gas.PressDropLin;

  Real dp = inlet.p - outlet.p;

end PressDropLin;

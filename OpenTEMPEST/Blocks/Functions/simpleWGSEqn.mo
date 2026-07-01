within OpenTEMPEST.Blocks.Functions;
function simpleWGSEqn

  extends Modelica.Icons.Function;
  import SI=Modelica.SIunits;

  input Real T "Temperature";
  input SI.MassFraction xIn[Medium.nX] "Mass composition";

  output SI.MassFraction xEq[Medium.nX]  "Mass composition in equilibrium";

protected
  SI.MoleFraction yEq[Medium.nX];
  SI.MoleFraction yIn[Medium.nX];
  Real C; // (* number of C atoms *)
  Real H; // (* number of H atoms *)
  Real O; // (* number of O atoms *)
  Real DeltaG;
  Real K_C;
  Real ha;
  Real hb;
  Real hc;
  Real n_out_H2;
  Real n_out_CO2;
  Real n_out_CO;
  Real n_out_H2O;
  Real n_out_N2;
  Real n_dot_out;
  Boolean C_NG;
  Boolean errorFlag;

protected
  replaceable package Medium = OpenTEMPEST.Medium.Fuel_CH4;

algorithm

  yIn :=Medium.massToMoleFractions(xIn, Medium.MMX);
  yIn := yIn/sum(yIn);

  // element sums
  C :=yIn[2] + yIn[3] + yIn[4];
  H :=2*yIn[1] + 4*yIn[2] + 2*yIn[5];
  O :=2*yIn[3] + yIn[4] + yIn[5];

  // reaction constant
  if T<500+273.15 then
    errorFlag :=true;
  else
    errorFlag :=false;
  end if;
  DeltaG :=32.1153*(T) - 3.5211E4;
  K_C :=exp(-DeltaG/Modelica.Constants.R/T);

  // help variables to solve the root function
  ha := K_C - 1;
  hb := K_C*(O - H - 2.0*C) - (O - 1.0/2.0*H - C);
  hc := 1.0/2.0*C*H*K_C - 1.0/2.0*H*K_C*(O - 1.0/2.0*H - C);

  // outlet molar masses
  if hb*hb - 4.0*ha*hc < 0 or ha == 0 then
    errorFlag := true;
  else
    n_out_H2 := (-hb - sqrt(hb*hb - 4.0*ha*hc))/2.0/ha;
    n_out_CO2 := n_out_H2 - 1.0/2.0*H - C + O;
    n_out_CO := C - n_out_CO2;
    n_out_H2O := 1.0/2.0*H - n_out_H2;
    n_out_N2 := yIn[6];
    n_dot_out := n_out_H2 + n_out_CO2 + n_out_CO + n_out_H2O + n_out_N2;
  end if;

  // outlet molar fraction
  if not errorFlag then
    yEq[3] := n_out_CO2/n_dot_out;
    yEq[1] := n_out_H2/n_dot_out;
    yEq[5] := n_out_H2O/n_dot_out;
    yEq[6] := n_out_N2/n_dot_out;
    yEq[4] := n_out_CO/n_dot_out;
    yEq[2] := 0;
  else
    yEq[:] := fill(-1, Medium.nX);
  end if;

  xEq :=Medium.moleToMassFractions(yEq, Medium.MMX);

    annotation (
    Documentation(info="<html>
    Calculated the water gas shift equilibrium (CH4 is consumed).
</html>",
        revisions="<html>
<ul>
<li><i>19-01-2021</i>
by <a href=\"mailto:marius.tomberg@dlr.de\">Marius Tomberg</a>:<br>
       Changed variable type from real to molar and mass fractions</li>
</ul>
</html>"));
end simpleWGSEqn;

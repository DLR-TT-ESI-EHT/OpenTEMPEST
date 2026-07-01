within OpenTEMPEST.Blocks.ComputationBlocks;
model WGS_MSR_Equilibrium
  "Model computing equilibrium considering WGS and MSR reactions, for above 500°C"

//   extends Modelica.Icons.Function;
  import SI = Modelica.SIunits;

//   input Real T "Temperature";
//   input SI.MassFraction xIn[Medium.nX] "Inlet composition, mass fractions";
//   input SI.Pressure p "pressure in Pa";
//
//   output SI.MassFraction xEq[Medium.nX]  "Composition in equilibrium, mass fractions";

  Modelica.Blocks.Interfaces.RealInput in_composition[Medium.nX] "mass fractions" annotation (
      Placement(transformation(
        origin={-50,28},
        extent={{-10,-10},{10,10}},
        rotation=0), iconTransformation(extent={{-60,12},{-40,32}})));
  Modelica.Blocks.Interfaces.RealInput in_T "Temperature in °C" annotation (Placement(
        transformation(
        origin={-50,-4},
        extent={{10,-10},{-10,10}},
        rotation=180),iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-50,-4})));
  Modelica.Blocks.Interfaces.RealInput in_p "pressure in Pa" annotation (Placement(transformation(
        origin={-50,-26},
        extent={{10,-10},{-10,10}},
        rotation=180), iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-50,-26})));

  Modelica.Blocks.Interfaces.RealOutput out_composition[Medium.nX] annotation (
      Placement(transformation(
        origin={70,0},
        extent={{-10,-10},{10,10}},
        rotation=0), iconTransformation(extent={{60,-10},{80,10}})));

  SI.MoleFraction yEq[Medium.nX];
  SI.MoleFraction yIn[Medium.nX];
  Real C; // (* number of C atoms *)
  Real H; // (* number of H atoms *)
  Real O; // (* number of O atoms *)
  Real DeltaG_WGS;
  Real DeltaG_MSR;
  Real K_WGS;
  Real K_MSR;
  Real nOut_H2;
  Real nOut_CH4;
  Real nOut_CO2;
  Real nOut_CO;
  Real nOut_H2O;
  Real nOut_N2;
  Real nOut;
//   Boolean C_NG;
  Boolean errorFlag;

protected
  package Medium = OpenTEMPEST.Medium.Fuel_CH4;

equation

  yIn = Medium.massToMoleFractions(in_composition[1:Medium.nXi], Medium.MMX);
//   yIn = yIn/sum(yIn);

  // element sums
  C =yIn[2] + yIn[3] + yIn[4];
  H =2*yIn[1] + 4*yIn[2] + 2*yIn[5];
  O =2*yIn[3] + yIn[4] + yIn[5];

  // reaction constant
  if in_T<500+273.15 then
    errorFlag =true;
  else
    errorFlag =false;
  end if;
  DeltaG_WGS =32.1153*(in_T) - 3.5211E4; // linear approximation
  K_WGS =exp(-DeltaG_WGS/Modelica.Constants.R/in_T);
  DeltaG_MSR =-252.642810968035*in_T + 225215.698063031; // linear approximation
  K_MSR =1.01325^2*10^10*exp(-DeltaG_MSR/Modelica.Constants.R/in_T);

  // set of equations to be solved:
  K_WGS*(nOut_CO * nOut_H2O) = (nOut_H2 * nOut_CO2);
  K_MSR*(nOut_H2O * nOut_CH4) = (nOut_H2^3 * nOut_CO) * in_p^2/nOut^2;
  C = nOut_CH4 + nOut_CO2 + nOut_CO;
  H = 2*nOut_H2 + 4*nOut_CH4 + 2*nOut_H2O;
  O = 2*nOut_CO2 + nOut_CO + nOut_H2O;
  nOut_N2 =  yIn[6];
  nOut =  nOut_H2 + nOut_CH4 + nOut_CO2 + nOut_CO + nOut_H2O + nOut_N2;

  // outlet molar fraction
  if not errorFlag then
    yEq[3] =  nOut_CO2/nOut;
    yEq[1] =  nOut_H2/nOut;
    yEq[5] =  nOut_H2O/nOut;
    yEq[6] =  nOut_N2/nOut;
    yEq[4] =  nOut_CO/nOut;
    yEq[2] =  nOut_CH4/nOut;
  else
    yEq[:] =  fill(-1, Medium.nX);
  end if;
//    xEq = Medium.moleToMassFractions(yEq, Medium.MMX);
  out_composition[1:Medium.nXi] =  Medium.moleToMassFractions(yEq, Medium.MMX);

    annotation (
    Documentation(info="<html>
    Calculates the equilibrium for given, yIn, T, p, considering water gas shift and methane reforming reactions.
</html>",
        revisions="<html>
<ul>
<li><i>24-05-2022</i>
by <a href=\"mailto:rene.lorenz@dlr.de\">Rene Lorenz</a>:<br>
Block form of TEMPEST.Blocks.Functions.simpleWgsMsrEqn</li>
</ul>
</html>"),
    Diagram(graphics={Rectangle(
          extent={{-40,40},{60,-40}},
          lineColor={102,44,145},
          fillColor={217,67,180},
          fillPattern=FillPattern.Backward)}),
    Icon(graphics={Rectangle(
          extent={{-40,40},{60,-40}},
          lineColor={102,44,145},
          fillColor={255,255,255},
          fillPattern=FillPattern.CrossDiag)}));
end WGS_MSR_Equilibrium;

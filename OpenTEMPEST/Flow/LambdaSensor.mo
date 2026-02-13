within OpenTEMPEST.Flow;
model LambdaSensor "\"Enter Desc\""
  import SI = Modelica.SIunits;
  replaceable package Fuel = Medium.Fuel_CH4;
  SI.MoleFraction x_out[Fuel.nX] "molar fractions in outgoing stream / -";
  SI.MoleFraction x_in[Fuel.nX] "molar fractions in ingoing stream / -";
  Real lambda_in "lambda-value of incoming flowstream / -";
  Real lambda_out "lambda-value of outgoing flowstream / -";
  SI.MolarFlowRate n_in "molar flow entering stack / (mole/s)";
  Real eps=ModelicaServices.Machine.small "smallest number, use instead of zero for comparison";
  Real Ke "electron gas coefficient / -";
  Real FUsys "Fuel utilization or Reactant conversion / -";
  Modelica.Blocks.Interfaces.RealInput xmass_out[Fuel.nX]
    annotation (Placement(transformation(extent={{140,40},{100,80}}),
        iconTransformation(extent={{120,32},{80,72}})));
  Modelica.Blocks.Interfaces.RealInput xmass_in[Fuel.nX]
    annotation (Placement(transformation(extent={{140,-60},{100,-20}}),
        iconTransformation(extent={{120,-62},{80,-22}})));
  Modelica.Blocks.Interfaces.IntegerInput N_cells annotation (Placement(
        transformation(extent={{-140,40},{-100,80}}),
                                                   iconTransformation(extent={{-120,28},
            {-80,68}})));
  Modelica.Blocks.Interfaces.RealInput I
    annotation (Placement(transformation(extent={{-140,-62},{-100,-22}}),
        iconTransformation(extent={{-120,-62},{-80,-22}})));
equation
  // convert input mass to mole fractions
  x_out = Fuel.massToMoleFractions(xmass_out, Fuel.MMX);
  x_in = Fuel.massToMoleFractions(xmass_in, Fuel.MMX);
  // lambda value as output of wide range oxygen sensor
  lambda_out=(x_out[4]+x_out[5]+2*x_out[3])/(2*x_out[4]+x_out[5]+x_out[1]+2*x_out[3]+4*x_out[2]);
  // lambda in the entering stream can be assumed known
  lambda_in =(x_in[4]+x_in[5]+2*x_in[3])/(2*x_in[4]+x_in[5]+x_in[1]+2*x_in[3]+4*x_in[2]);
  // calculate the gas electron coefficient for fuel cell or elctrolysis mode
  if I>=0 then Ke=8*x_in[2]+2*x_in[4]+2*x_in[1] "gas electron coefficient in fuel cell mode /-";
  else Ke=2*x_in[3]+2*x_in[5] "gas electron coefficient in electrolysis mode /-";
  end if;
  if abs(I)>eps and abs(lambda_out-lambda_in)>eps then // only defined when there is net current and lambda out is different than lambda in
    n_in=I*N_cells/(2*Modelica.Constants.F*(2*x_in[4]+x_in[5]+x_in[1]+2*x_in[3]+4*x_in[2])*(lambda_out-lambda_in));
    // calculate also fuel utilization / reactant conversion for fuel cell or electrolysis
    FUsys=2*(2*x_in[4]+x_in[5]+x_in[1]+2*x_in[3]+4*x_in[2])/Ke*abs(lambda_out-lambda_in);
  else
    n_in=0;
    FUsys=0;
  end if;
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-120,
            -100},{120,100}}),                                  graphics={
        Rectangle(
          extent={{-80,80},{80,-80}},
          lineColor={28,108,200},
          fillColor={102,44,145},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-78,68},{-20,30}},
          lineColor={244,125,35},
          fillColor={102,44,145},
          fillPattern=FillPattern.None,
          textString="Nc"),
        Text(
          extent={{-94,-22},{-32,-64}},
          lineColor={244,125,35},
          fillColor={102,44,145},
          fillPattern=FillPattern.None,
          textString="I"),
        Text(
          extent={{-2,76},{72,26}},
          lineColor={244,125,35},
          fillColor={102,44,145},
          fillPattern=FillPattern.None,
          textString="Xout"),
        Text(
          extent={{4,-14},{76,-68}},
          lineColor={244,125,35},
          fillColor={102,44,145},
          fillPattern=FillPattern.None,
          textString="Xin"),
        Text(
          extent={{-44,22},{42,-16}},
          lineColor={244,125,35},
          fillColor={102,44,145},
          fillPattern=FillPattern.None,
          textString="lambda")}),                                Diagram(
        coordinateSystem(preserveAspectRatio=false, extent={{-120,-100},{120,
            100}})),
    experiment(
      StopTime=2000000,
      __Dymola_NumberOfIntervals=200,
      Tolerance=2e-05,
      __Dymola_Algorithm="Cvode"));
end LambdaSensor;

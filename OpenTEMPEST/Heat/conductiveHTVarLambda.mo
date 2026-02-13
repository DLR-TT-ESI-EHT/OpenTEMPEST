within OpenTEMPEST.Heat;
model conductiveHTVarLambda
  extends ThermoPower.Icons.HeatFlow;
  import SI = Modelica.SIunits;
  SI.ThermalConductivity lambda "thermalConductivity";
  SI.Temperature T;
  parameter SI.Length L "length";
  parameter SI.Area A "cross sectional area";
  parameter Enumerations.InsulationOptions matOpt
    "Options for insulation materials";
  parameter Boolean matLin = false "c and lambda = linear(f(T))";
  parameter Boolean matConst = false "c and lambda = f(TStart) = const.";
  parameter Real Fac=1 "Correction factor for thermal conductivity";

  ThermoPower.Thermal.HT side1 annotation (Placement(transformation(extent={{-40,
            20},{40,40}}, rotation=0)));
  ThermoPower.Thermal.HT side2 annotation (Placement(transformation(extent={{-40,
            -20},{40,-42}}, rotation=0)));

 // parameter Modelica.SIunits.ThermalConductivity lambda=
 //   if matOpt == InsulationOptions.Ultra then Fac*(2.9352e-8*T^2 - 1.175e-5*T + 0.01881) else 5;

equation
  T = 0.5*(side1.T + side2.T);

  if not matConst then
    if matOpt == Enumerations.InsulationOptions.Ultra then
      if matLin then
        lambda = Fac*(3.989E-05*T + 7.666E-03);
      else
        lambda = Fac*(2.9352e-8*T^2 - 1.175e-5*T + 0.01881);
      end if;
    elseif matOpt == Enumerations.InsulationOptions.Shape then
      lambda = Fac*(2.884e-5*T + 0.01141);
    elseif matOpt == Enumerations.InsulationOptions.Duratec then
      lambda = Fac*(5e-05*T + 0.3363);
    else
      assert(false, "Unsupported material properties option");
      lambda = 0;
    end if;
  else
    lambda = 0.03;
  end if;

  side1.Q_flow = lambda*A*(side1.T - side2.T)/L "Conductive heat transfer";
  side1.Q_flow = -side2.Q_flow "Energy balance";
end conductiveHTVarLambda;

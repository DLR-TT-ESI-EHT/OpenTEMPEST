within OpenTEMPEST.Heat;
model ForcedConvection "Heat transfer correlation for forced convection"
  extends ThermoPower.Thermal.BaseClasses.DistributedHeatTransferFV;
  //! User MUST average temperature to be true for this model!
  parameter Boolean pipe=true
    "= true if flow through a pipe, false if between two parallel plates";
  final parameter Real Re_max=1e6 "Maximum Reynolds number";
  Real Pr_min[Nf] "Minimum Prandtl number";
  //pipe: 0.1 for turbulent flow, 0.6 for transition zone; parallel plates: 0.5 for turbulent, 0.1 for laminar flow
  final parameter Real Pr_max=1000 "Maximum Prandtl number";
  //1000 for all cases but the turbulent flow between plates: 100
  Modelica.SIunits.CoefficientOfHeatTransfer gamma[Nf]
    "Heat transfer coefficients at the nodes";
  Modelica.SIunits.CoefficientOfHeatTransfer gamma_vol[Nw]
    "Heat transfer coefficients at the volumes";
  Medium.Temperature Tvol[Nw] "Fluid temperature in the volumes";
  Medium.DynamicViscosity mu[Nf] "Dynamic viscosity";
  Medium.ThermalConductivity k[Nf] "Thermal conductivity";
  Medium.SpecificHeatCapacity cp[Nf] "Heat capacity at constant pressure";
  Modelica.SIunits.PerUnit Re[Nf] "Reynolds number";
  Modelica.SIunits.PerUnit Pr[Nf] "Prandtl numbers";
  Modelica.SIunits.PerUnit Nu[Nf] "Nusselt numbers";
  Real tmp1[Nf];
  Real tmp2[Nf];

equation
  assert(Nw == Nf - 1, "Number of volumes Nw on wall side should be equal to number of volumes fluid side Nf - 1");
  // Fluid properties at the nodes

  mu[:] = Medium.dynamicViscosity(fluidState[:]);
  k[:] = Medium.thermalConductivity(fluidState[:]);
  cp[:] = Medium.heatCapacity_cp(fluidState[:]);
  //Re[:] = abs(w[:].*Dhyd./(A*mu[:]));
  //Pr[:] = cp[:].*mu[:]./k[:];

  assert(useAverageTemperature == true, "Calculation must be done with average temperature");

  for j in 1:Nf loop
    // mu[j] = Medium.dynamicViscosity(fluidState[j]);
    //k[j] = Medium.thermalConductivity(fluidState[j]);
    //cp[j] = Medium.heatCapacity_cp(fluidState[j]);
    Re[j] = abs(w[j]*Dhyd/(A*mu[j]));
    Pr[j] = cp[j]*mu[j]/k[j];

    assert(
      Pr[j] <= Pr_max,
      "Prandtl number exceeds full accuracy range",
      AssertionLevel.warning);
    assert(
      Pr[j] >= Pr_min[j],
      "Prandtl number underruns full accuracy range",
      AssertionLevel.warning);
    assert(
      Re[j] <= Re_max,
      "Reynolds number exceeds full accuracy range",
      AssertionLevel.warning);
    assert(
      Dhyd/L <= 1,
      "Ratio of hydraulic diameter to lenght exceeds limit",
      AssertionLevel.warning);
    if pipe then

      Pr_min[j] = smooth(0, if Re[j] < 2301 then 0.1 elseif Re[j] > 10000 then 0.1
         else 0.6);
      Nu[j] = (smooth(0, if Re[j] < 2302 then (49.371 + (1.6152*(Re[j]*Pr[j]*
        Dhyd/L)^(1/3) - 0.7)^3 + ((2/(1 + 22*Pr[j]))^(1/6)*sqrt(Re[j]*Pr[j]*
        Dhyd/L))^3)^(1/3) elseif Re[j] > 10000 then tmp2[j] else tmp1[j]));
      gamma[j] = Nu[j]*k[j]/Dhyd;

      tmp1[j] = ((1 - (Re[j] - 2300)/(1e4 - 2300))*(49.371 + (1.615*(2300*Pr[j]*
        Dhyd/L)^(1/3) - 0.7)^3 + ((2/(1 + 22*Pr[j]))^(1/6)*sqrt(2300*Pr[j]*Dhyd/
        L))^3)^(1/3) + ((Re[j] - 2300)/(1e4 - 2300))*(0.0308/8*1e4*Pr[j]/(1 + 12.7
        *sqrt(0.0308/8)*(Pr[j]^(2/3) - 1))*(1 + (Dhyd/L)^(2/3))));
      tmp2[j] = (((1.8*log10(Re[j]+0.001) - 1.5)^(-2))/8*Re[j]*Pr[j]/(1 + 12.7*sqrt(((
        1.8*log10(Re[j]+0.001) - 1.5)^(-2))/8)*(Pr[j]^(2/3) - 1))*(1 + (Dhyd/L)^(2/3)));

    else
      //Mind that Dhyd=2*width of parallel plates!
      if Re[j] < 2200 then
        Pr_min[j] = 0.1;
        gamma[j] = k[j]/Dhyd*(7.541^3 + (1.841*(Re[j]*Pr[j]*Dhyd/L)^(1/3))^3)^(1
          /3);
        //for constant and equal wall temperature conditions on both sides, hydrodynamically developed flow
      elseif Re[j] > 40000 then
        Pr_min[j] = 0.5;
        gamma[j] = k[j]/Dhyd*(((1.8*log10(Re[j]) - 1.5)^(-2))/8*Re[j]*Pr[j]/(1 +
          12.7*sqrt(((1.8*log10(Re[j]) - 1.5)^(-2))/8)*(Pr[j]^(2/3) - 1))*(1 + (
          Dhyd/L)^(2/3)));
        //for constant wall temperature condition, Dhyd/L<=1
      else
        Pr_min[j] = 0.5;
        //assumption as no records shown
        gamma[j] = k[j]/Dhyd*((1 - (Re[j] - 2200)/(4e4 - 2200))*((3.66 + (4 - 0.102
          /1.2))^3 + (1.615*1.14*(Re[j]*Pr[j]*Dhyd/L)^(1/3))^3 + ((2/(1 + 22*Pr[
          j]))^(1/6)*sqrt(2300*Pr[j]*Dhyd/L))^3)^(1/3) + ((Re[j] - 2200)/(4e4 -
          2200))*(0.0217/8*4e4*Pr[j]/(1 + 12.7*sqrt(0.0217/8)*(Pr[j]^(2/3) - 1))
          *(1 + (Dhyd/L)^(2/3))));
        //for both constant wall temperature or constant heat flux, Dhyd/L<=1
      end if;
    end if;
  end for;
  for j in 1:Nw loop
    Tvol[j] = if useAverageTemperature then (T[j] + T[j + 1])/2 else T[j + 1];
    gamma_vol[j] = if useAverageTemperature then (gamma[j] + gamma[j + 1])/2
       else gamma[j + 1];
    Qw[j] = (Tw[j] - Tvol[j])*omega*l*gamma_vol[j]*Nt;
  end for;
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics), Icon(graphics={Text(extent={{-100,-52},
              {100,-80}}, textString="%name")}));
end ForcedConvection;

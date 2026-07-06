within OpenTEMPEST.Medium;
package SimpleExhaust "TEMPEST: Exhaust gas [O2, CO2, H2O, N2]"
  extends Modelica.Media.IdealGases.Common.MixtureGasNasa(
    mediumName = "NaturalGasMixture",
    data = {Modelica.Media.IdealGases.Common.SingleGasesData.O2,
            Modelica.Media.IdealGases.Common.SingleGasesData.CO2,
            Modelica.Media.IdealGases.Common.SingleGasesData.H2O,
            Modelica.Media.IdealGases.Common.SingleGasesData.N2},
    fluidConstants = {Modelica.Media.IdealGases.Common.FluidData.O2,
                      Modelica.Media.IdealGases.Common.FluidData.CO2,
                      Modelica.Media.IdealGases.Common.FluidData.H2O,
                      Modelica.Media.IdealGases.Common.FluidData.N2},
    substanceNames = {"Oxygen", "Carbondioxide", "Water", "Nitrogen"},
    reference_X = {0.05, 0.45, 0.3, 0.2},
    referenceChoice= Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
    excludeEnthalpyOfFormation = false);
end SimpleExhaust;

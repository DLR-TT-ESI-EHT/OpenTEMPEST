within OpenTEMPEST.Medium;
package Fuel_C3 "TEMPEST: FT off gas with C3"
  extends Modelica.Media.IdealGases.Common.MixtureGasNasa(
    mediumName = "NaturalGasMixture",
    data = {Modelica.Media.IdealGases.Common.SingleGasesData.H2,
            Modelica.Media.IdealGases.Common.SingleGasesData.CH4,
            Modelica.Media.IdealGases.Common.SingleGasesData.CO2,
            Modelica.Media.IdealGases.Common.SingleGasesData.CO,
            Modelica.Media.IdealGases.Common.SingleGasesData.H2O,
            Modelica.Media.IdealGases.Common.SingleGasesData.N2,
            Modelica.Media.IdealGases.Common.SingleGasesData.C3H8},
    fluidConstants = {Modelica.Media.IdealGases.Common.FluidData.H2,
                      Modelica.Media.IdealGases.Common.FluidData.CH4,
                      Modelica.Media.IdealGases.Common.FluidData.CO2,
                      Modelica.Media.IdealGases.Common.FluidData.CO,
                      Modelica.Media.IdealGases.Common.FluidData.H2O,
                      Modelica.Media.IdealGases.Common.FluidData.N2,
                      Modelica.Media.IdealGases.Common.FluidData.C3H8},
    substanceNames = {"Hydrogen", "Methane", "Carbondioxide", "Carbonmonoxide", "Water", "Nitrogen", "Propane"},
    reference_X = {0.2, 0.49, 0.05, 0.05, 0.2, 0.0, 0.01},
    referenceChoice= Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
    excludeEnthalpyOfFormation = false);
end Fuel_C3;

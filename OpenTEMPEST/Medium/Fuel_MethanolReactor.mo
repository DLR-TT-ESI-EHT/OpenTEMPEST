within OpenTEMPEST.Medium;
package Fuel_MethanolReactor "TEMPEST: Methane and methanol gas mixture"
  extends Modelica.Media.IdealGases.Common.MixtureGasNasa(
    mediumName = "NaturalGasMixture",
    data = {Modelica.Media.IdealGases.Common.SingleGasesData.H2,
            Modelica.Media.IdealGases.Common.SingleGasesData.CH4,
            Modelica.Media.IdealGases.Common.SingleGasesData.CO2,
            Modelica.Media.IdealGases.Common.SingleGasesData.CO,
            Modelica.Media.IdealGases.Common.SingleGasesData.H2O,
            Modelica.Media.IdealGases.Common.SingleGasesData.N2,
            Modelica.Media.IdealGases.Common.SingleGasesData.CH3OH},
    fluidConstants = {Modelica.Media.IdealGases.Common.FluidData.H2,
                      Modelica.Media.IdealGases.Common.FluidData.CH4,
                      Modelica.Media.IdealGases.Common.FluidData.CO2,
                      Modelica.Media.IdealGases.Common.FluidData.CO,
                      Modelica.Media.IdealGases.Common.FluidData.H2O,
                      Modelica.Media.IdealGases.Common.FluidData.N2,
                      Modelica.Media.IdealGases.Common.FluidData.CH3OH},
    substanceNames = {"Hydrogen", "Methane", "Carbondioxide", "Carbonmonoxide", "Water", "Nitrogen", "Methanol"},
    reference_X = {0.126,1.00E-05,0.297, 0.57696,1.00E-05, 1.00E-05, 1.00E-05},
    referenceChoice= Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
    excludeEnthalpyOfFormation = false);                                        //{0.5, 0.04, 0.2, 0.2, 0.05, 0.005, 0.005},

end Fuel_MethanolReactor;

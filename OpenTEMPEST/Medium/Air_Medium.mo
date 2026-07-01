within OpenTEMPEST.Medium;
package Air_Medium "TEMPEST: Ideal air model (O2/N2)"
  extends Modelica.Media.IdealGases.Common.MixtureGasNasa(
    mediumName = "AirMixture",
    data = {Modelica.Media.IdealGases.Common.SingleGasesData.O2,
            Modelica.Media.IdealGases.Common.SingleGasesData.N2},
    fluidConstants = {
            Modelica.Media.IdealGases.Common.FluidData.O2,
            Modelica.Media.IdealGases.Common.FluidData.N2},
    substanceNames = {"Oxygen", "Nitrogen"},
    reference_X = {0.2329, 0.7671},
    referenceChoice= Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
    excludeEnthalpyOfFormation = false);
    //referenceChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.UserDefined,
end Air_Medium;

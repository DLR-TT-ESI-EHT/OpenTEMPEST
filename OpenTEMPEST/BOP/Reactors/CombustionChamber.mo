within OpenTEMPEST.BOP.Reactors;
model CombustionChamber "Combustion Chamber"
  extends ThermoPower.Gas.BaseClasses.CombustionChamberBase(
    redeclare package Air = Medium.Air_Medium "O2, N2",
    redeclare package Fuel = Medium.Fuel_CH4 "H2, CH4, CO2, CO, H2O, N2",
    redeclare package Exhaust = Medium.FlueGas "O2, Ar, H2O, CO2, N2");
  Real nfFuel(final quantity="MolarFlowRate", unit="mol/s")
    "Molar Combustion rate";
  Modelica.SIunits.PerUnit lambda
    "Stoichiometric ratio (>1 if air flow is greater than stoichiometric)";
protected
  Air.MassFraction ina_X[Air.nXi]=inStream(ina.Xi_outflow);
  Fuel.MassFraction inf_X[Fuel.nXi]=inStream(inf.Xi_outflow);
  Fuel.MoleFraction inf_Y[Fuel.nXi]=Fuel.massToMoleFractions(inf_X, Fuel.MMX);
equation
  //wcomb = inf.m_flow*inf_X[3]/Fuel.data[3].MM "Combustion molar flow rate";
  nfFuel = inf.m_flow/sum(inf_X*Fuel.MMX);
  lambda = (ina.m_flow*ina_X[1]/Air.data[1].MM)/(0.5*nfFuel*inf_Y[1]+2*nfFuel*inf_Y[2]+0.5*nfFuel*inf_Y[4]);
  assert(lambda >= 1, "Not enough oxygen flow");
  der(MX[1]) = ina.m_flow*ina_X[1] + out.m_flow*fluegas.X[1] - (0.5*nfFuel*inf_Y[1]+2*nfFuel*inf_Y[2]+0.5*nfFuel*inf_Y[4])*
    Exhaust.data[1].MM "oxygen";
  der(MX[2]) = out.m_flow*fluegas.X[2] "argon";
  der(MX[3]) = out.m_flow*fluegas.X[3] + (nfFuel*inf_Y[1]+2*nfFuel*inf_Y[2]+nfFuel*inf_Y[5])*
    Exhaust.data[3].MM "water";
  der(MX[4]) = out.m_flow*fluegas.X[4] + (nfFuel*inf_Y[3]+nfFuel*inf_Y[4])*Exhaust.data[
    4].MM "carbondioxide";
  der(MX[5]) = ina.m_flow*ina_X[2] + out.m_flow*fluegas.X[5] + inf.m_flow*
    inf_X[6] "nitrogen";
  annotation (Icon(graphics), Documentation(info="<html>
This model extends the CombustionChamber Base model, with the definition of the gases.
<p>In particular, the air inlet uses the <tt>Media.Air</tt> medium model, the fuel input uses the <tt>Media.NaturalGas</tt> medium model, and the flue gas outlet uses the <tt>Medium.FlueGas</tt> medium model.
<p>The composition of the outgoing gas is determined by the mass balance of every component, taking into account the combustion reaction CH4+2O2--->2H2O+CO2.</p>
<p>The model assumes complete combustion, so that it is only valid if the oxygen flow at the air inlet is greater than the stoichiometric flow corresponding to the flow at the fuel inlet.</p>

</html>",
        revisions="<html>
<ul>
<li><i>31 Jan 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
 Combustion Chamber model restructured using inheritance.
     <p>  First release.
 </li>
</ul>
</html>"));
end CombustionChamber;

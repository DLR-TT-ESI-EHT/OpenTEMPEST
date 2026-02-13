within OpenTEMPEST.Heat.BaseClasses;
partial model Element1D
  "Partial heat transfer element with two 1D thermal ports that does not store energy"

  parameter Integer N "Number of control volumes";
  Modelica.SIunits.HeatFlowRate Q[N]
    "Heat flow rate from dHT0 -> dHT1";
  Modelica.SIunits.TemperatureDifference dT[N] "dHT0.T - dHT1.T";
public
  ThermoPower.Thermal.DHTVolumes dHT0(N=N) annotation (Placement(
        transformation(extent={{-102,-6},{-90,6}}), iconTransformation(extent={{-102,-6},
            {-90,6}})));
  ThermoPower.Thermal.DHTVolumes dHT1(N=N) annotation (Placement(
        transformation(extent={{90,-6},{102,6}}), iconTransformation(extent={{90,-6},
            {102,6}})));
equation

    dT = dHT0.T .- dHT1.T;
    dHT0.Q = Q;
    dHT1.Q = -Q;

  annotation (Documentation(info="<html>
<p>
This partial model contains the basic connectors and variables to
allow heat transfer models to be created that <strong>do not store energy</strong>,
This model defines and includes equations for the temperature
drop across the element, <strong>dT</strong>, and the heat flow rate
through the element from dHT0 to dHT1, <strong>Q</strong>.
</p>
<p>
By extending this model, it is possible to write simple
constitutive equations for many types of heat transfer components.
</p>
</html>"));
end Element1D;

within OpenTEMPEST.Heat;
model PrescribedTemperature1D
  "Variable temperature boundary condition in Kelvin"
  extends ThermoPower.Icons.HeatFlow;
  parameter Integer N "Number of control volumes";
  ThermoPower.Thermal.DHTVolumes port(N=N) annotation (Placement(transformation(
          extent={{-40,-40},{40,-20}}),
                                      iconTransformation(extent={{-40,-40},{40,-20}})));
  Modelica.Blocks.Interfaces.RealInput T[N](unit="K") annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
        rotation=270,
        origin={0,40})));
equation
    port.T = T;
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
            100,100}}), graphics={
                   Text(
          extent={{-102,-44},{98,-68}},
          lineColor={191,95,0},
          textString="%name")}),
    Documentation(info="<html>
<p>
This model represents a variable temperature boundary condition.
The temperature in [K] is given as input signal <strong>T</strong>
to the model. The effect is that an instance of this model acts as
an infinite reservoir able to absorb or generate as much energy
as required to keep the temperature at the specified value.
</p>
</html>"),
       Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
            100,100}}), graphics={
                   Text(
          extent={{-100,-44},{100,-68}},
          lineColor={191,95,0},
          textString="%name")}));
end PrescribedTemperature1D;

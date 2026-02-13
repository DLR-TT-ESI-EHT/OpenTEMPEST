within OpenTEMPEST.Heat;
model ThermalR_HeatC_FV
  "Modeling thermal resistance between an FV object and an object with lumped T and lumped thermal element storing heat"
  import SI = Modelica.SIunits;
  parameter Integer N(min=1) "Number of nodes of DHT";
  parameter SI.ThermalResistance R "Constant thermal resistance";
  parameter SI.HeatCapacity C "Heat capacity of element (= cp*m)";
  parameter SI.Temperature Tstart
    "Start value of the temperature of element";
  SI.Temperature T(start=Tstart, displayUnit="degC")
    "Temperature of element";
  ThermoPower.Thermal.DHTVolumes dHT(final N = N) annotation (Placement(transformation(extent={{-10,
            -40},{10,-20}}), iconTransformation(extent={{-10,-38},{10,-18}})));
  ThermoPower.Thermal.HT hT annotation (Placement(transformation(extent={{-10,20},
            {10,40}}), iconTransformation(extent={{-10,22},{10,42}})));
equation
  T = hT.T;
  C * der(T) = sum(dHT.Q)  + hT.Q_flow "Energy balance";
  for i in 1:N loop
    dHT.T[i] - hT.T = R * dHT.Q[i];
  end for;
  annotation (Icon(graphics={                                                                            Rectangle(extent={{
              -80,22},{80,-18}},                                                                                                    lineColor = {0, 0, 0}, fillColor = {135, 135, 135},
            fillPattern =                                                                                                   FillPattern.CrossDiag)}));
end ThermalR_HeatC_FV;

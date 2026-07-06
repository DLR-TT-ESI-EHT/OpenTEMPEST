within OpenTEMPEST.Heat;
model ConductiveFV "Conductive thermal resistance without mass"
  import SI = Modelica.SIunits;
  parameter Integer N(min=1) = 2 "Number of nodes";
 // parameter SI.SpecificHeatCapacity cap
 //   "Specific heat capacity of metal";
  parameter SI.ThermalConductivity lambda "thermalConductivity";
 // Modelica.SIunits.Temperature Tmea;
  parameter SI.Length L "length between heat ports";
  parameter SI.Area A "cross sectional area";
  //parameter SI.Density rho "Density of element";

 //outer ThermoPower.System system "System wide properties";
  //Units.AbsoluteTemperature T[N](start=Tstart) "Node temperatures";
  SI.HeatFlowRate Q_m "mean heat flux of side 1";
  //SI.HeatFlux phi_2 "mean heat flux of side 2";
  ThermoPower.Thermal.DHTVolumes side1(N=N) "Internal surface"
    annotation (Placement(transformation(extent={{-40,20},{40,40}}, rotation=
            0), iconTransformation(extent={{-40,20},{40,40}})));
  ThermoPower.Thermal.DHTVolumes side2(N=N) "External surface"
    annotation (Placement(transformation(extent={{-40,-42},{40,-20}},
          rotation=0), iconTransformation(extent={{-40,-42},{40,-20}})));
equation

    side1.Q = lambda * A/N * (side1.T - side2.T) / L
    "Conductive heat transfer";

    side2.Q =  -side1.Q
    "Conductive heat transfer";

    Q_m = sum(side1.Q[1:N])/N;

  annotation (Icon(graphics={
          Rectangle(
          extent={{-92,20},{92,-20}},
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid,
          pattern=LinePattern.None)}));
end ConductiveFV;

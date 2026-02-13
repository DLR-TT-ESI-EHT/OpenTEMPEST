within OpenTEMPEST.Blocks.Sources;
model PulseVar
  Real frequency = 2*Modelica.Constants.pi/duration;
  Real out;

  Modelica.Blocks.Interfaces.RealOutput y annotation (Placement(transformation(extent={{100,-10},{120,10}}), iconTransformation(extent={{100,-10},{120,10}})));

  Modelica.Blocks.Interfaces.RealInput duration annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={40,80}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={40,80})));
  Modelica.Blocks.Interfaces.RealInput duty annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={-40,80}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={-40,80})));
  Modelica.Blocks.Interfaces.RealInput Amplitude annotation (Placement(transformation(extent={{-124,-20},{-84,20}}), iconTransformation(extent={{-124,-20},{-84,20}})));
  Modelica.Blocks.Interfaces.BooleanInput pulseOn annotation (Placement(transformation(
        extent={{20,-20},{-20,20}},
        rotation=-90,
        origin={0,-100}), iconTransformation(
        extent={{20,-20},{-20,20}},
        rotation=-90,
        origin={0,-80})));
equation

  if not pulseOn then
    out=Amplitude;
  else
    out =Amplitude*OpenTEMPEST.Blocks.Functions.pulseWave(
      time,
      duration,
      duty);
  end if;

  y = out;

  annotation (Icon(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}}), graphics={
                                Rectangle(
        extent={{-100,-100},{100,100}},
        lineColor={0,0,127},
        fillColor={255,255,255},
        fillPattern=FillPattern.Solid),
        Line(points={{-80,68},{-80,-80}}, color={192,192,192}),
        Polygon(
          points={{-80,90},{-88,68},{-72,68},{-80,90}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(points={{-90,-70},{82,-70}}, color={192,192,192}),
        Polygon(
          points={{90,-70},{68,-62},{68,-78},{90,-70}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(points={{-80,-70},{-40,-70},{-40,44},{0,44},{0,-70},{40,-70},{40,
              44},{79,44}}),
        Text(
          extent={{-147,-152},{153,-112}},
          textString="period=%period")}), Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=10,
      Tolerance=1e-06,
      __Dymola_fixedstepsize=0.0001,
      __Dymola_Algorithm="Dassl"));
end PulseVar;

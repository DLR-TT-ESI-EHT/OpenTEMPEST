within OpenTEMPEST.SOC.Electrochem.Interfaces;
connector SimpFlowInfoPort

  Modelica.SIunits.Temperature TAir;
  Modelica.SIunits.Temperature TFuel;
  Modelica.SIunits.SpecificHeatCapacityAtConstantPressure cpAir;
  Modelica.SIunits.SpecificHeatCapacityAtConstantPressure cpFuel;
  Modelica.SIunits.MassFlowRate mfAir;
  Modelica.SIunits.MassFlowRate mfFuel;

  Modelica.SIunits.Temperature TAirRef;
  Modelica.SIunits.Temperature TFuelRef;
  Modelica.SIunits.EnergyFlowRate HfAirRef;
  Modelica.SIunits.EnergyFlowRate HfFuelRef;

  annotation (Icon(graphics={                                                                                                              Polygon(
          origin={9.74,0},
          fillColor={230,230,230},
          fillPattern=FillPattern.Solid,
          points={{-89.7391,98},{-89.7391,-98},{90.2609,0},{90.2609,0},{-89.7391,
              98}})}));
end SimpFlowInfoPort;

within OpenTEMPEST.BOP.GasConditioning;
model SimpleGasHeater
  extends Flow.BaseClasses.simpleGasHeaterBase;
  import SI = Modelica.SIunits;

   medium.BaseProperties gas(p(start = pStart, fixed=false), T(start = TStart, fixed=false), Xi(start = xStart[1:medium.nXi], fixed=false));

  parameter Integer nHT = 1 "number of heat port in flow direction";

  ThermoPower.Thermal.HT ht[nHT] annotation (Placement(transformation(extent={{-20,
            20},{0,40}}), iconTransformation(extent={{-20,20},{20,60}})));

protected
  SI.MassFraction xIn[medium.nXi] " MAss flow rates entering the reactors";
  SI.MassFraction xOut[medium.nXi] "Mass flow rates from reactor to mixing gas node";
  SI.SpecificEnthalpy hIn;
  SI.SpecificEnthalpy hOut;

equation

  if not useAlphaIn then

    eta = medium.dynamicViscosity(gas.state); // 18.2e-6;//
    cp = medium.specificHeatCapacityCp(gas.state); //  1;//
    lambda = medium.thermalConductivity(gas.state); // 0.0262; //

  else
    eta = -1;
    cp = -1;
    lambda = -1;
  end if;

  // heat flow calculation
  ht[:].Q_flow = alpha .* A .* (ht[:].T-fill(gas.T,nHT));

  // Balances
  0 = inlet.m_flow + outlet.m_flow; //mass Balamce
  0 = inlet.m_flow * hIn + outlet.m_flow * hOut + sum(ht[:].Q_flow);
  xIn = gas.Xi;

  //connecting variables with connectors
  inlet.p = gas.p;
  outlet.p = gas.p;
  hIn = homotopy(if not allowFlowReversal then inStream(inlet.h_outflow) else
    actualStream(inlet.h_outflow), inStream(inlet.h_outflow));
  hOut = homotopy(if not allowFlowReversal then gas.h else actualStream(outlet.h_outflow),
    gas.h);
  xIn = inStream(inlet.Xi_outflow);
  xOut = outlet.Xi_outflow;
  outlet.h_outflow = gas.h;
  xOut = gas.Xi;

  inlet.h_outflow = gas.h;
  inlet.Xi_outflow = medium.reference_X[1:medium.nXi];

  annotation (Icon(graphics={Rectangle(
          extent={{-60,20},{60,-20}},
          lineColor={28,108,200},
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid),
        Text(extent={{-100,-30},{100,-54}},  textString="%name")}));
end SimpleGasHeater;

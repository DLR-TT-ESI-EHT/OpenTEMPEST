within OpenTEMPEST.BOP;
model PumpSimple "simple pump model for incompressible fluids"
    extends ThermoPower.Icons.Water.Pump;
    import SI = Modelica.SIunits;

    replaceable package Medium = Modelica.Media.Water.StandardWater
    constrainedby Modelica.Media.Interfaces.PartialMedium "Medium model"  annotation(choicesAllMatching = true);
    Medium.BaseProperties inletFluid(h(start=hStart))
      "Fluid properties at the inlet";

    Medium.MassFlowRate mf=in_mf "Mass flow rate (total)";
    Medium.SpecificEnthalpy h(start=hStart) "Fluid specific enthalpy";
    Medium.SpecificEnthalpy hIn(start=hStart) "Enthalpy of entering fluid";
    Medium.SpecificEnthalpy hOut(start=hStart) "Enthalpy of outgoing fluid";
    SI.Power P_el "Power Consumption";
    parameter Real eta_isen = 0.8 "isentropic efficency";
    parameter Real eta_el = 0.95 "electrical efficency";
    parameter Real eta_mech = 0.9 "mechanical efficency";
    parameter Medium.SpecificEnthalpy hStart=1e5  "Fluid Specific Enthalpy Start Value"  annotation (Dialog(tab="Initialisation"));
    parameter Boolean allowFlowReversal=system.allowFlowReversal
      "= true to allow flow reversal, false restricts to design direction"  annotation(Evaluate=true);
    outer ThermoPower.System system "System wide properties";

    ThermoPower.Water.FlangeA fl_inlet(redeclare package Medium = Medium,
      h_outflow(start=hStart),
      m_flow(min=if allowFlowReversal then -Modelica.Constants.inf else 0))
      annotation (Placement(transformation(extent={{-100,2},{-60,42}}, rotation=0)));
    ThermoPower.Water.FlangeB fl_outlet(redeclare package Medium = Medium,
      h_outflow(start=hStart),
      m_flow(max=if allowFlowReversal then +Modelica.Constants.inf else 0))
      annotation (Placement(transformation(extent={{40,52},{80,92}}, rotation=0)));

    Modelica.Blocks.Interfaces.RealInput in_mf annotation (Placement(
          transformation(
          origin={-40,76},
          extent={{-10,-10},{10,10}},
          rotation=270)));
    ThermoPower.Electrical.PowerConnection powerConnection annotation (Placement(
          transformation(extent={{-10,80},{10,100}}), iconTransformation(extent={{
              -10,80},{10,100}})));
    ThermoPower.Electrical.Load load(Pnom=1, usePowerInput=true)
      annotation (Placement(transformation(extent={{-10,46},{10,66}})));
    Modelica.Blocks.Sources.RealExpression realExpression(y=P_el)
      annotation (Placement(transformation(extent={{-28,46},{-8,66}})));
equation
    // electrical power requiremed
    P_el = mf/inletFluid.d* (fl_outlet.p - fl_inlet.p) /(eta_isen*eta_el*eta_mech);

     // Mass balance
    fl_inlet.m_flow + fl_outlet.m_flow = 0 "Mass balance";

    // Energy balance
    fl_outlet.h_outflow = inStream(fl_inlet.h_outflow) + P_el/(mf+Modelica.Constants.small)  "Energy balance for mf > 0";
    fl_inlet.h_outflow = inStream(fl_outlet.h_outflow) + P_el/(mf+Modelica.Constants.small)  "Energy balance for mf < 0";
    h = homotopy(if not allowFlowReversal then fl_outlet.h_outflow else if mf >= 0
          then fl_outlet.h_outflow else fl_inlet.h_outflow, fl_outlet.h_outflow)
         "Definition of h";

    // Fluid properties (always uses the properties upstream of the inlet flange)
    inletFluid.p = fl_inlet.p;
    inletFluid.h = inStream(fl_inlet.h_outflow);
    inletFluid.Xi = inStream(fl_inlet.Xi_outflow);

    // Boundary conditions
    mf = fl_inlet.m_flow "Pump total flow rate";
    hIn = homotopy(if not allowFlowReversal then inStream(fl_inlet.h_outflow) else
      if mf >= 0 then inStream(fl_inlet.h_outflow) else h, inStream(fl_inlet.h_outflow));
    hOut = homotopy(if not allowFlowReversal then h else if mf >= 0 then h else
      inStream(fl_outlet.h_outflow), h);
    fl_inlet.Xi_outflow = inStream(fl_outlet.Xi_outflow);
    inStream(fl_inlet.Xi_outflow) = fl_outlet.Xi_outflow;

    connect(powerConnection, load.port) annotation (Line(
        points={{0,90},{0,64.6}},
        color={0,0,255},
        thickness=0.5));
    connect(load.referencePower, realExpression.y) annotation (Line(points={{-3.3,56},{-7,56}}, color={0,0,127}));

    annotation (
      Documentation(revisions="<html>
<ul>
<li><i>10 Jan 2022</i>
by <a href=\"mailto:rene.lorenz@dlr.de\">René Lorenz</a>:<br>
      Adapted/Simplified from <tt>TEMPEST.Flow.FanSimple</tt> model.</li>
</ul>
</html>", info="<html>
<p>Model for a simple pump, based on assumption: </p>
<ol>
<li> incompressible fluid (d<code>rho<\\code>/d<code>p<c\\ode> = 0) </li>
<li>no internal Volume in pump</li>
<li>adiabatic pump</li>
</ol>
<p> 
<ul>
 P = mf (h2 - h1) <br>
 eta_isen= dh_isen/dh_poly, with h_isen = vdp <br>
 P= mf * dh_poly = mf * dh_isen/eta_isen = mf *vdp/eta_isen<br>
 --> resource image with equations to be inserted
</ul>


</html>"));
end PumpSimple;

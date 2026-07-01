within OpenTEMPEST.BOP.Reactors.MeOH;
model LiquidSeparator
  "Model that allows for representing plant components that separate specific components of a mixture flow"

  import SI = Modelica.SIunits;
  replaceable package GasMix = Medium.Fuel_MethanolReactor;
//   GasMix.BaseProperties Gas(
//     p(start=pStart, stateSelect=StateSelect.prefer),
//     T(start=Tstart, stateSelect=StateSelect.prefer),
//     Xi(start=xStart[1:GasMix.nXi], stateSelect=StateSelect.prefer));
  parameter Integer iLiq = 5 "number i of component in Medium.Xi to be separated from mixture";
  parameter SI.Temperature TOut = 50+273.15 "start temperature for iteration";
  parameter SI.SpecificEnthalpy hf0Liq = -15874780 "std specific enthalpy of formation of liquid";
  parameter SI.SpecificHeatCapacity cpLiq = 4.229 "specific heat cpacity of liquid @TOut";
  parameter SI.SpecificEnthalpy dhvLiq = 2256000 " [J/kg] latent heat of evaporation liquid at respective T,p";
  SI.MassFlowRate mf "total mass flow rate through separator";
   SI.SpecificEnthalpy hIn "specific enthalpy of inflowing gas mixture mass flow";
   SI.SpecificEnthalpy hGas "specific enthalpy of gas mixture outflow";
   SI.SpecificEnthalpy hLiq "specific enthalpy of separated liquid outflow";
  SI.MassFraction XIn[GasMix.nXi];
  SI.MassFraction XGas[GasMix.nXi](start= GasMix.reference_X);
  SI.MassFlowRate mfLiq "mass flow rate of separated pure component";
  SI.MassFlowRate miGas[GasMix.nXi] "mass outflow rate of mixture gas";
  SI.MassFlowRate mfGas;
//   SI.MassFlowRate mfGasAbs;

  inner ThermoPower.System system annotation (Placement(transformation(extent={{-100,80},{-80,100}})));
  ThermoPower.Gas.FlangeA infl(redeclare package Medium = GasMix) "Inlet port connector for gas entering the reactor"  annotation (Placement(transformation(extent={{-100,-20},{-60,20}}), iconTransformation(extent={{-100,
            -20},{-60,20}})));
  ThermoPower.Gas.FlangeB outflGas(redeclare package Medium = GasMix) "Outlet port connector of gas mixture stream leaving separator"  annotation (Placement(transformation(
          extent={{60,10},{100,50}}), iconTransformation(extent={{60,10},{
            100,50}})));
  ThermoPower.Water.FlangeB outflLiquid(redeclare package Medium =
        ThermoPower.Water.StandardWater)                                                            "Outlet port connector of single component gas stream leaving separator"  annotation (Placement(transformation(extent={{60,-50},{100,-10}}),
        iconTransformation(extent={{60,-50},{100,-10}})));

    ThermoPower.Thermal.HT hT annotation (Placement(transformation(extent={{-16,60},
             {16,92}}), iconTransformation(extent={{-12,50},{12,74}})));
equation

  // mass balance
   0 = mf + mfLiq + mfGas;
   0 = XIn[iLiq]*mf + mfLiq;
   //    mfGas = sum(miGas[:]);

// component balance:
//   mfGasAbs = if mfGas >= 0 then mfGas else -mfGas; // else (if noEvent(mfGas >=0) then mfGas else -mfGas);
  for i in 1:GasMix.nXi loop
    XGas[i] = miGas[i]/(abs(mfGas)+1e-10);
    if i <> iLiq then
      miGas[i] = XIn[i]*mf;
    elseif i == iLiq then
      miGas[i] = 0;
    end if;
  end for;

//   for i in 1:GasMix.nXi loop
//     XGas[i] = miGas[i]/mfGas;
//   end for;

// energy balance
  hLiq = hf0Liq + cpLiq*(TOut-298.15);  // H0 + cp*dT
  hGas = GasMix.h_TX(T=TOut,X=XGas);
  -hT.Q_flow= (mf*hIn + mfGas*hGas + mfLiq*(hLiq + dhvLiq));

//connecting variables with connectors
  outflGas.p= infl.p;
// outflLiquid.p= infl.p;

  infl.m_flow = mf;
  outflGas.m_flow = mfGas;
  outflLiquid.m_flow = mfLiq;

  hIn = inStream(infl.h_outflow);
  hIn = infl.h_outflow;
  hGas = outflGas.h_outflow;
  hLiq = outflLiquid.h_outflow;

  XIn = inStream(infl.Xi_outflow);
  XIn = infl.Xi_outflow;
  XGas = outflGas.Xi_outflow;
//   XLiq = outflLiquid.Xi_outflow;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
          extent={{-60,50},{60,-50}},
          fillColor={140,236,255},
          fillPattern=FillPattern.HorizontalCylinder,
          pattern=LinePattern.None,
          lineColor={0,0,0}),        Text(
          extent={{-50,36},{50,-32}},
          lineColor={0,0,0},
          textString="liquid
separator")}),                                                   Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end LiquidSeparator;

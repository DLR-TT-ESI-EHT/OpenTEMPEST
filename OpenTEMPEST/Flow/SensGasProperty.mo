within OpenTEMPEST.Flow;
model SensGasProperty
  extends OpenTEMPEST.Flow.SensGasFlow;

  Medium.BaseProperties gas(p(start = pstart, fixed=false), T(start = Tstart, fixed=false), Xi(start = Xstart[1:Medium.nXi], each fixed=false));
  parameter Medium.AbsolutePressure pstart=1e5 "Pressure start value"
    annotation (Dialog(tab="Initialisation"));
  parameter Medium.Temperature Tstart=300 "Temperature start value"
    annotation (Dialog(tab="Initialisation"));
  parameter Medium.MassFraction Xstart[Medium.nX]=Medium.reference_X
    "Start gas composition" annotation (Dialog(tab="Initialisation"));

  parameter Boolean cpOutput = true "Sensor has output for the heat capacity";
  parameter Boolean lambdaOutput = true "Sensor has output for the thermal conductivity";
  parameter Boolean etaOutput = true "Sensor has output for the dynamic viscosity";
  parameter Boolean TOutput = true "Sensor has output for the temperature";
  parameter Boolean rhoOutput = true "Sensor has output for the density";

  Modelica.Blocks.Interfaces.RealOutput cp if cpOutput "specificHeatCapacityCp"
                                                     annotation (Placement(
        transformation(extent={{80,68},{100,88}},rotation=0),
        iconTransformation(extent={{80,90},{100,110}})));
  Modelica.Blocks.Interfaces.RealOutput lambda
                                           if lambdaOutput
    "thermalConductivity"                                  annotation (Placement(
        transformation(extent={{80,50},{100,70}}, rotation=0),
        iconTransformation(extent={{80,70},{100,90}})));
  Modelica.Blocks.Interfaces.RealOutput eta             "dynamicViscosity"
    annotation (Placement(transformation(extent={{80,30},{100,50}}, rotation=0),
        iconTransformation(extent={{80,10},{100,30}})));
  Modelica.Blocks.Interfaces.RealOutput T(unit="K") if TOutput       "Temperature"
    annotation (Placement(transformation(extent={{80,10},{100,30}}, rotation=0),
        iconTransformation(extent={{80,50},{100,70}})));
  Modelica.Blocks.Interfaces.RealOutput rho if rhoOutput   "density"
    annotation (Placement(transformation(extent={{82,-10},{102,10}}, rotation=0),
        iconTransformation(extent={{80,30},{100,50}})));
protected
  Modelica.Blocks.Sources.RealExpression realExpression_cp(y=if cpOutput then
        Medium.specificHeatCapacityCp(gas.state) else -1)
    annotation (Placement(transformation(extent={{16,82},{36,102}})));
  Modelica.Blocks.Sources.RealExpression realExpression_lambda(y=if
        lambdaOutput then Medium.thermalConductivity(gas.state) else -1)
    annotation (Placement(transformation(extent={{20,58},{40,78}})));
  Modelica.Blocks.Sources.RealExpression realExpression_eta(y=if
        etaOutput then Medium.dynamicViscosity(gas.state) else -1)
    annotation (Placement(transformation(extent={{18,38},{38,58}})));
  Modelica.Blocks.Sources.RealExpression realExpression_T(y=if TOutput then gas.T
         else -1)
    annotation (Placement(transformation(extent={{18,20},{38,40}})));
  Modelica.Blocks.Sources.RealExpression realExpression_rho(y=if rhoOutput
         then gas.d else -1)
    annotation (Placement(transformation(extent={{18,2},{38,22}})));
equation

  // Set gas properties
  inlet.p = gas.p;
  gas.h = homotopy(if not allowFlowReversal then inStream(inlet.h_outflow)
     else inStream(inlet.h_outflow), inStream(inlet.h_outflow));
  gas.Xi = homotopy(if not allowFlowReversal then inStream(inlet.Xi_outflow)
 else
     inStream(inlet.Xi_outflow), inStream(inlet.Xi_outflow));

//

  connect(realExpression_cp.y, cp)
    annotation (Line(points={{37,92},{90,92},{90,78}}, color={0,0,127}));
  connect(realExpression_lambda.y, lambda)
    annotation (Line(points={{41,68},{90,68},{90,60}}, color={0,0,127}));
  connect(realExpression_eta.y, eta)
    annotation (Line(points={{39,48},{90,48},{90,40}}, color={0,0,127}));
  connect(rho, rho) annotation (Line(points={{92,0},{92,0}}, color={0,0,127}));
  connect(realExpression_T.y, T)
    annotation (Line(points={{39,30},{90,30},{90,20}}, color={0,0,127}));
  connect(realExpression_rho.y, rho)
    annotation (Line(points={{39,12},{92,12},{92,0}}, color={0,0,127}));
end SensGasProperty;

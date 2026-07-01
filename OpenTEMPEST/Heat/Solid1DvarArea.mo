within OpenTEMPEST.Heat;
model Solid1DvarArea
  "solid 1D volume with area variation along the x-axis (truncated pyramid with base area in y-z-plane)"

  import SI = Modelica.SIunits;
  replaceable package SolidMat = OpenTEMPEST.Solid.SolidMatBase annotation(choicesAllMatching = true);

  SolidMat.BaseProperties Solid[N](T=T, each kCustom_trans=kCustom_trans, each kCustom_long=kCustom_long, each rhoCustom=rhoCustom, each cpCustom=cpCustom);

  parameter Integer N(min=3) = 5 "Number of CVs in the solid";

  // Initial Values
  parameter SI.Temperature TstartIn=773.15 annotation(Dialog(group="Initialisation"));
  parameter SI.Temperature TstartOut=773.15 annotation(Dialog(group="Initialisation"));
  parameter SI.Temperature Tstart[N] = linspace(TstartIn, TstartOut, N) annotation(Dialog(group="Initialisation"));

  // Dimensions
  parameter SI.Length lX = 1 "Total length of solid" annotation(Dialog(group="Dimensions"));

  SI.Length dx =  lX /N "Length of a CV";
  parameter SI.Area A0 = 1 "CS area at x/L=0" annotation(Dialog(group="Dimensions"));
  parameter SI.Area A1 = 1 "CS area at x/L=1" annotation(Dialog(group="Dimensions"));
  SI.Length dA_dx =  (A1-A0)/lX "dA/dx = gradient of area change in x-direction";

  SI.Length delta_x[N+1] "Total length from origin";
  SI.Area Ax[N+1] "Cross sectional areas at CV interfaces";
  SI.Volume dV[N] "Volume of CV";

  SI.Temperature T[N] "Average Temperature in CV";
  SI.ThermalConductivity kv_trans[N-1] "harmonic mean of conductivites between CVs";
  // SI.ThermalConductivity kv_long[N-1] "harmonic mean of conductivites between CVs";

  // Solid Properties
  SI.ThermalConductivity k_trans[N] = Solid[:].k_trans;
  SI.ThermalConductivity k_long[N] = Solid[:].k_long;
  SI.Density rho[N] = Solid[:].rho;
  SI.SpecificHeatCapacity cp[N] = Solid[:].cp;

  parameter SI.ThermalConductivity kCustom_trans = 1 annotation(Dialog(group="Custom Material Only"));
  parameter SI.ThermalConductivity kCustom_long = kCustom_trans annotation(Dialog(group="Custom Material Only"));
  parameter SI.Density rhoCustom = 1 annotation(Dialog(group="Custom Material Only"));
  parameter SI.SpecificHeatCapacity cpCustom = 1 annotation(Dialog(group="Custom Material Only"));

  // Thermal ports
  ThermoPower.Thermal.HT hT_x0 annotation (Placement(transformation(extent={{-120,
            -10},{-100,10}}),
                         iconTransformation(extent={{-120,-10},{-100,10}})));
  ThermoPower.Thermal.HT hT_xN annotation (Placement(transformation(extent={{100,-10},
            {120,10}}), iconTransformation(extent={{100,-10},{120,10}})));
  ThermoPower.Thermal.DHTVolumes dhT_int(N=N)   annotation (Placement(transformation(extent={{-10,-10},
            {10,10}}),  iconTransformation(extent={{-10,-10},{10,10}})));

initial equation
  T[:] = Tstart[:];

equation
  //cross sectional area in x-plane
  delta_x[:] = linspace(0.0,N,N+1).*dx;
  Ax = A0*ones(N+1) .+ dA_dx .* delta_x;
  dV[:] = ( Ax[1:N]+Ax[2:N+1] + sqrt(Ax[1:N].*Ax[2:N+1])) .*dx/3;

  kv_trans[:] = 2.*k_trans[1:N-1].*k_trans[2:N]./(k_trans[1:N-1] .+ k_trans[2:N]);
  // kv_long[:] = 2.*k_long[1:N-1].*k_long[2:N]./(k_long[1:N-1] .+ k_long[2:N]);

  rho[1]*cp[1]*dV[1]*der(T[1])                    = hT_x0.Q_flow + kv_trans[1]*Ax[2]*(T[2] - T[1])/dx + dhT_int.Q[1];
  rho[2:N-1].*cp[2:N-1].*dV[2:N-1].*der(T[2:N-1]) = kv_trans[1:N-2].*Ax[2:N-1].*(T[1:N-2] .- T[2:N-1])./dx .+ kv_trans[2:N-1].*Ax[3:N].*(T[3:N] .- T[2:N-1])./dx .+ dhT_int.Q[2:N-1];
  rho[N]*cp[N]*dV[N]*der(T[N])                    = hT_xN.Q_flow + kv_trans[N-1]*Ax[N]*(T[N-1] - T[N])/dx + dhT_int.Q[N];

  //Boundary Conditions
  hT_x0.Q_flow = (k_trans[1]/(0.5*dx))*A0 * (hT_x0.T - T[1]); // Conductive Flow In
  hT_xN.Q_flow = (k_trans[N]/(0.5*dx))*A1 * (hT_xN.T - T[N]); // Conductive Flow Out
  dhT_int.T[:] = T[:];

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Polygon(
          points={{32,66},{100,80},{100,-80},{32,-66},{32,8},{32,66}},
          lineColor={28,108,200},
          fillColor={116,116,116},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-100,40},{-34,54},{-34,-52},{-100,-40},{-100,12},{-100,40}},
          lineColor={28,108,200},
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-98,8},{-82,-16}},
          lineColor={238,46,47},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="x0"),
        Text(
          extent={{82,12},{98,-10}},
          lineColor={238,46,47},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="xN"),
        Polygon(
          points={{-34,54},{32,66},{32,-66},{-34,-52},{-34,26},{-34,54}},
          lineColor={162,29,33},
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
          Text(
            extent={{-139,23},{139,-23}},
            lineColor={255,255,255},
          textString="%name",
          origin={43,-1},
          rotation=90)}),                                           Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    Documentation(revisions="<html>
</html>", info="<html>
<p>truncated pyramid volume with heat transfer along x-direction. Pyramid base area is y0*z0.</p>
</html>"));
end Solid1DvarArea;

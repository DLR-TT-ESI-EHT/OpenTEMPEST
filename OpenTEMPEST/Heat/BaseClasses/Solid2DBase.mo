within OpenTEMPEST.Heat.BaseClasses;
partial model Solid2DBase

  import SI = Modelica.SIunits;
  replaceable package SolidMat = TEMPEST.Solid.SolidMatBase annotation(choicesAllMatching = true);

  parameter Integer nX(min=3) = 5 "Number of CVs in the solid in x-direction";
  parameter Integer nY(min=3) = 5 "Number of CVs in the solid in y-direction";

  // Initial Values
  parameter Modelica.SIunits.Temperature Tstartbar=1073.15 "uniform start temperature" annotation(Dialog(group="Initialisation"));

  // Dimensions
  parameter SI.Length lX = 1  "Total length of solid" annotation(Dialog(group="Dimensions"));
  parameter SI.Length lY = 1 "width of solid" annotation(Dialog(group="Dimensions"));
  parameter SI.Length lZ = 1 "thickness of solid" annotation(Dialog(group="Dimensions"));

  Real dx = lX/nX "x-Length of a CV";
  Real dy = lY/nY "y-Length of a CV";
  Real dV = dx*dy*lZ "Volume of a CV";
  Real Ax = dy*lZ "cross sectional area in x plane";
  Real Ay = dx*lZ "cross sectional area in y-plane";
  Real Az = dx*dy "cross sectional area in z-plane";

  SI.Temperature T[nX, nY];

  SolidMat.BaseProperties Solid[nX,nY](
    T=T,
    each kCustom_trans=kCustom_trans,
    each kCustom_long=kCustom_long,
    each rhoCustom=rhoCustom,
    each cpCustom=cpCustom);

  // Solid Properties
  SI.ThermalConductivity k_trans[nX, nY] = Solid.k_trans "Thermal conductivities in CV";
  SI.ThermalConductivity k_long[nX, nY] = Solid.k_long "Thermal conductivities in CV";
  SI.Density rho[nX, nY] = Solid.rho;
  SI.SpecificHeatCapacity cp[nX, nY] = Solid.cp;

  SI.ThermalConductivity kx[nX-1, nY]  "Thermal conductivities CV";
  SI.ThermalConductivity ky[nX, nY-1]  "Thermal conductivities between CV";

  SI.HeatFlowRate Qext[nX, nY];

  // Custom material properties
  parameter SI.ThermalConductivity kCustom_trans = 1 "Thermal Conductivity across layers" annotation(Dialog(group="Custom Material Only"));
  //parameter SI.ThermalConductivity kCustom_long = 1 "Thermal Conductivity in plane of layers (=k_trans for homogeneous materials)" annotation(Dialog(group="Custom Material Only"));
  parameter SI.ThermalConductivity kCustom_long = kCustom_trans annotation(Dialog(group="Custom Material Only"));
  parameter SI.Density rhoCustom = 1 annotation(Dialog(group="Custom Material Only"));
  parameter SI.SpecificHeatCapacity cpCustom = 1 annotation(Dialog(group="Custom Material Only"));

  // Thermal ports
  ThermoPower.Thermal.DHTVolumes dhT_x0(N=nY) annotation (Placement(transformation(extent={{-120,
            -10},{-100,10}}),
                       iconTransformation(extent={{-120,-10},{-100,10}})));
  ThermoPower.Thermal.DHTVolumes dhT_xN(N=nY) annotation (Placement(transformation(extent={{90,-10},
            {110,10}}), iconTransformation(extent={{100,-10},{120,10}})));
  OpenTEMPEST.Heat.DHTVolumes2D dhT2_z1(i=nX, j=nY) annotation (Placement(
        transformation(extent={{-80,60},{-60,80}}), iconTransformation(extent={
            {-80,60},{-60,80}})));
  OpenTEMPEST.Heat.DHTVolumes2D dhT2_z0(i=nX, j=nY) annotation (Placement(
        transformation(extent={{-80,-80},{-60,-60}}), iconTransformation(extent
          ={{-80,-80},{-60,-60}})));
  ThermoPower.Thermal.DHTVolumes dhT_y0(N=nX) annotation (Placement(transformation(extent={{60,-80},
            {80,-60}}),iconTransformation(extent={{60,-80},{80,-60}})));
  ThermoPower.Thermal.DHTVolumes dhT_yN(N=nX) annotation (Placement(transformation(extent={{60,60},
            {80,80}}),  iconTransformation(extent={{60,60},{80,80}})));
  OpenTEMPEST.Heat.DHTVolumes2D dhT2_int(i=nX, j=nY) annotation (Placement(
        transformation(extent={{-10,-10},{10,10}}), iconTransformation(extent={
            {-10,-10},{10,10}})));

initial equation
  T = fill(Tstartbar, nX, nY);

equation
  for i in 1:nX-1 loop
    kx[i, :]  = 2.*k_long[i,:] .*k_long[i+1,:] ./(k_long[i,:]  .+ k_long[i+1,:]);
  end for;
  for i in 1:nY-1 loop
    ky[:, i]  = 2.*k_long[:,i] .*k_long[:,i+1] ./(k_long[:,i]  .+ k_long[:,i+1]);
  end for;

  // Energy Balance

  // Internal Nodes
  rho[2:nX-1, 2:nY-1].*cp[2:nX-1, 2:nY-1].*dV.*der(T[2:nX-1, 2:nY-1]) = Ax./dx.*(kx[1:nX-2, 2:nY-1].*(T[1:nX-2, 2:nY-1] .- T[2:nX-1, 2:nY-1]) .+ kx[2:nX-1, 2:nY-1].*(T[3:nX, 2:nY-1] .- T[2:nX-1, 2:nY-1]))  .+
                                                                        Ay./dy.*(ky[2:nX-1, 1:nY-2].*(T[2:nX-1, 1:nY-2] .- T[2:nX-1, 2:nY-1]) .+ ky[2:nX-1, 2:nY-1].*(T[2:nX-1, 3:nY] .- T[2:nX-1, 2:nY-1]))  .+
                                                                        dhT2_z0.Q[2:nX-1, 2:nY-1] .+ dhT2_z1.Q[2:nX-1, 2:nY-1] .+ dhT2_int.Q[2:nX-1, 2:nY-1] .+ Qext[2:nX-1, 2:nY-1];

    // Corner Nodes
  rho[1, 1]*cp[1,1]*dV*der(T[1,1])       = dhT_x0.Q[1]  + kx[1,1]    *Ax/dx*(T[2,1] - T[1,1])        + dhT_y0.Q[1]    + ky[1,1]    *Ay/dy*(T[1,2] - T[1,1])    + dhT2_z0.Q[1,1]   + dhT2_z1.Q[1,1]   + dhT2_int.Q[1,1] + Qext[1,1];
  rho[1, nY]*cp[1,nY]*dV*der(T[1,nY])    = dhT_x0.Q[nY] + kx[1,nY]   *Ax/dx*(T[2,nY] - T[1,nY])      + dhT_yN.Q[1]  + ky[1,nY-1] *Ay/dy*(T[1,nY-1] - T[1,nY])  + dhT2_z0.Q[1,nY]  + dhT2_z1.Q[1,nY]  + dhT2_int.Q[1,nY] + Qext[1,nY];
  rho[nX, 1]*cp[nX,1]*dV*der(T[nX,1])    = dhT_xN.Q[1]  + kx[nX-1,1] *Ax/dx*(T[nX-1,1] - T[nX,1])    + dhT_y0.Q[nX]  + ky[nX,1]   *Ay/dy*(T[nX,2] - T[nX,1])   + dhT2_z0.Q[nX,1]  + dhT2_z1.Q[nX,1]  + dhT2_int.Q[nX,1] + Qext[nX,1];
  rho[nX, nY]*cp[nX,nY]*dV*der(T[nX,nY]) = dhT_xN.Q[nY] + kx[nX-1,nY]*Ax/dx*(T[nX-1, nY] - T[nX,nY]) + dhT_yN.Q[nX] + ky[nX,nY-1]*Ay/dy*(T[nX,nY-1] - T[nX,nY])+ dhT2_z0.Q[nX,nY] + dhT2_z1.Q[nX,nY] + dhT2_int.Q[nX,nY] + Qext[nX,nY];

  // Edge Nodes

  //front (left) face [1,:]
  rho[1, 2:nY-1].*cp[1, 2:nY-1].*dV.*der(T[1, 2:nY-1])    = Ay./dy.*(ky[1, 1:nY-2] .*(T[1, 1:nY-2] .- T[1, 2:nY-1]) .+ ky[1, 2:nY-1] .*(T[1, 3:nY].- T[1, 2:nY-1]))  .+ dhT_x0.Q[2:nY-1] .+ kx[1, 2:nY-1].*Ax./dx   .*(T[2,  2:nY-1] .- T[1, 2:nY-1])    .+ dhT2_z0.Q[1, 2:nY-1]   .+ dhT2_z1.Q[1, 2:nY-1]  .+ dhT2_int.Q[1, 2:nY-1] .+ Qext[1, 2:nY-1];

  //back (right) face [nX,:]
  rho[nX, 2:nY-1].*cp[nX, 2:nY-1].*dV.*der(T[nX, 2:nY-1]) = Ay./dy.*(ky[nX, 1:nY-2].*(T[nX, 1:nY-2].- T[nX, 2:nY-1]).+ ky[nX, 2:nY-1].*(T[nX, 3:nY].- T[nX, 2:nY-1])).+ dhT_xN.Q[2:nY-1] .+ kx[nX-1, 2:nY-1].*Ax./dx.*(T[nX-1,  2:nY-1] .- T[nX, 2:nY-1]).+ dhT2_z0.Q[nX, 2:nY-1]  .+ dhT2_z1.Q[nX, 2:nY-1] .+ dhT2_int.Q[nX, 2:nY-1] .+ Qext[nX, 2:nY-1];

  //front (right)face [:,1]
  rho[2:nX-1, 1].*cp[2:nX-1, 1].*dV.*der(T[2:nX-1, 1])    = Ax./dx.*(kx[1:nX-2, 1] .*(T[1:nX-2, 1] .- T[2:nX-1, 1]) .+ kx[2:nX-1, 1] .*(T[3:nX, 1].- T[2:nX-1, 1]))  .+ dhT_y0.Q[2:nX-1] .+ ky[2:nX-1, 1].*Ay./dy   .*(T[2:nX-1, 2] .- T[2:nX-1, 1])   .+ dhT2_z0.Q[2:nX-1, 1]   .+ dhT2_z1.Q[2:nX-1, 1]  .+ dhT2_int.Q[2:nX-1, 1] .+ Qext[2:nX-1, 1];

  //back face [:,nY]
  rho[2:nX-1, nY].*cp[2:nX-1, nY].*dV.*der(T[2:nX-1, nY]) = Ax./dx.*(kx[1:nX-2, nY].*(T[1:nX-2, nY].- T[2:nX-1, nY]).+ kx[2:nX-1, nY].*(T[3:nX, nY].- T[2:nX-1, nY])).+ dhT_yN.Q[2:nX-1] .+ ky[2:nX-1, nY-1].*Ay./dy.*(T[2:nX-1,  nY-1] .- T[2:nX-1, nY]).+ dhT2_z0.Q[2:nX-1, nY]  .+ dhT2_z1.Q[2:nX-1, nY] .+ dhT2_int.Q[2:nX-1, nY] .+ Qext[2:nX-1, nY];

  // BCs - known heat flux
  //z-plane
  dhT2_z0.Q[:,:] = (k_trans[:,:] ./(0.5*lZ))*Az  .* (dhT2_z0.T[:,:] .- T[:,:]);  // Conductive Flow bottom
  dhT2_z1.Q[:,:] = (k_trans[:,:] ./(0.5*lZ))*Az  .* (dhT2_z1.T[:,:] .- T[:,:]); // Conductive Flow top
  // y-plane
  dhT_y0.Q[:]   = (k_long[:,1] ./(0.5*dy))*Ay .* (dhT_y0.T[:]   .- T[:,1]); // Conductive Flow Left
  dhT_yN.Q[:]   = (k_long[:,nY]./(0.5*dy))*Ay .* (dhT_yN.T[:]   .- T[:,nY]); // Conductive Flow Right
  //x-plane
  dhT_x0.Q[:]   = (k_long[1,:] ./(0.5*dx))*Ax .* (dhT_x0.T[:]   .- T[1,:]);
  dhT_xN.Q[:]   = (k_long[nX,:]./(0.5*dx))*Ax .* (dhT_xN.T[:]   .- T[nX,:]);

  // Internal
  dhT2_int.T = T;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Rectangle(
          extent={{-100,60},{100,20}},
          lineColor={28,108,200},
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-100,-20},{100,-60}},
          lineColor={28,108,200},
          fillColor={116,116,116},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-100,20},{100,-20}},
          lineColor={162,29,33},
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-76,64},{-62,40}},
          lineColor={255,240,0},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="z1"),
        Text(
          extent={{-78,-38},{-60,-62}},
          lineColor={255,240,0},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="z0"),
        Text(
          extent={{62,-34},{80,-62}},
          lineColor={244,125,35},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="y0"),
        Text(
          extent={{62,62},{80,40}},
          lineColor={244,125,35},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="yN"),
          Text(
            extent={{-138,-18},{142,-58}},
            lineColor={255,255,255},
          textString="%name"),
        Text(
          extent={{-96,12},{-78,-10}},
          lineColor={244,125,35},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="x0"),
        Text(
          extent={{82,10},{98,-10}},
          lineColor={244,125,35},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="xN")}),                                         Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    Documentation(revisions="<html>
<ul>
<li><i>27 Aug 2021</i> by <a href=\"faisal.sedeqi@dlr.de\">Faisal Sedeqi</a>:<br>Modified Conduction Between CVs including Direction Dependent Conductivity. </li>
<li><i>03 Aug 2021</i> by <a href=\"faisal.sedeqi@dlr.de\">Faisal Sedeqi</a>:<br>Changed port names. Will change conduction between CVs next.</li>
<li><i>23 Jul 2021</i> by <a href=\"hans.wiggenhauser@dlr.de\">Hans Wiggenhauser</a>:<br>finished and tested.</li>
<li><i>20 Jul 2021</i> by <a href=\"faisal.sedeqi@dlr.de\">Faisal Sedeqi</a>:<br>created, edge energy balance missinge. </li>
</ul>
</html>"));
end Solid2DBase;

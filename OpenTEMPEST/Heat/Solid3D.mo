within OpenTEMPEST.Heat;
model Solid3D

  import SI = Modelica.SIunits;
  replaceable package SolidMat = TEMPEST.Solid.Material.Custom  constrainedby
    TEMPEST.Solid.SolidMatBase                                                                            annotation(choicesAllMatching = true);

  // Initial Values
  parameter Modelica.SIunits.Temperature Tstartbar=1073.15 "uniform start temperature" annotation(Dialog(group="Initialisation"));

  // Dimensions
  parameter SI.Length lX=10e-2    "Total length of solid" annotation(Dialog(group="Dimensions"));
  parameter SI.Length lY=10e-2  "width of solid" annotation(Dialog(group="Dimensions"));
  parameter SI.Length lZ=10e-2                               "thickness of solid" annotation(Dialog(group="Dimensions"));
  parameter Integer nX=3 "number of discretisation elements in x-disrection" annotation(Dialog(group="Dimensions"));
  parameter Integer nY=3 "number of discretisation elements in y-disrection" annotation(Dialog(group="Dimensions"));
  parameter Integer nZ=3 "number of discretisation elements in z-disrection" annotation(Dialog(group="Dimensions"));

  Real dx =  lX/nX "Length of a CV";
  Real dy =  lY/nY "Length of a CV";
  Real dz =  lZ/nZ "Length of a CV";
  Real dV = dx*dy*dz "Volume of a CV";
  Real Ax = dy*dz "coss sectional area perpendicular to x";
  Real Ay = dx*dz "coss sectional area perpendicular to y";
  Real Az = dx*dy "coss sectional area perpendicular to z";
  SI.Temperature T[nX,nY,nZ];

  // Solid Properties
  SolidMat.BaseProperties Solid[nX,nY,nZ](
    T=T,
    each kCustom_trans=kCustom_trans,
    each kCustom_long=kCustom_long,
    each rhoCustom=rhoCustom,
    each cpCustom=cpCustom);
  SI.ThermalConductivity k_trans[nX, nY, nZ] = Solid.k_trans "Thermal conductivities in CV";
  SI.ThermalConductivity k_long[nX, nY, nZ] = Solid.k_long "Thermal conductivities in CV";
  SI.Density rho[nX,nY,nZ]=Solid.rho;
  SI.SpecificHeatCapacity cp[nX,nY,nZ]=Solid.cp;

  SI.ThermalConductivity kx[nX-1, nY, nZ]  "Thermal conductivities between CV";
  SI.ThermalConductivity ky[nX, nY-1, nZ]  "Thermal conductivities between CV";
  SI.ThermalConductivity kz[nX, nY, nZ-1]  "Thermal conductivities between CV";

  parameter SI.ThermalConductivity kCustom_trans = 1 "Thermal Conductivity across layers" annotation(Dialog(group="Custom Material Only"));
  parameter SI.ThermalConductivity kCustom_long = 1 "Thermal Conductivity in plane of layers (=k_trans for homogeneous materials)" annotation(Dialog(group="Custom Material Only"));
  parameter SI.Density rhoCustom=-1 "only used if SolidMat=Custom" annotation(Dialog(group="Custom Material Only"));
  parameter SI.SpecificHeatCapacity cpCustom = -1 "only used if SolidMat=Custom" annotation(Dialog(group="Custom Material Only"));

  // Thermal ports
  OpenTEMPEST.Heat.DHTVolumes2D dhT2_x0(i=nY, j=nZ) annotation (Placement(
        transformation(extent={{-100,-10},{-80,10}}), iconTransformation(extent
          ={{-120,-10},{-100,10}})));
  OpenTEMPEST.Heat.DHTVolumes2D dhT2_xN(i=nY, j=nZ) annotation (Placement(
        transformation(extent={{90,-10},{110,10}}), iconTransformation(extent={
            {100,-10},{120,10}})));
  OpenTEMPEST.Heat.DHTVolumes2D dhT2_zN(i=nX, j=nY) annotation (Placement(
        transformation(extent={{-80,60},{-60,80}}), iconTransformation(extent={
            {-80,60},{-60,80}})));
  OpenTEMPEST.Heat.DHTVolumes2D dhT2_z0(i=nX, j=nY) annotation (Placement(
        transformation(extent={{-80,-80},{-60,-60}}), iconTransformation(extent
          ={{-80,-80},{-60,-60}})));
  OpenTEMPEST.Heat.DHTVolumes2D dhT2_y0(i=nX, j=nZ) annotation (Placement(
        transformation(extent={{60,-80},{80,-60}}), iconTransformation(extent={
            {60,-80},{80,-60}})));
  OpenTEMPEST.Heat.DHTVolumes2D dhT2_yN(i=nX, j=nZ) annotation (Placement(
        transformation(extent={{60,60},{80,80}}), iconTransformation(extent={{
            60,60},{80,80}})));

initial equation
  T[:] = fill(Tstartbar, nX,nY,nZ); //Tstart[:];

equation
  for i in 1:nX - 1 loop
    kx[i, :, :] = 2.*k_long[i, :, :] .* k_long[i + 1, :, :] ./ (k_long[i, :, :] .+ k_long[i + 1, :, :]);
  end for;
  for i in 1:nY - 1 loop
    ky[:, i, :] = 2.*k_long[:, i, :] .* k_long[:, i + 1, :] ./ (k_long[:, i, :] .+ k_long[:, i + 1, :]);
  end for;
  for i in 1:nZ - 1 loop
    kz[:, :, i] = 2.*k_trans[:, :, i] .* k_trans[:, :, i + 1] ./ (k_trans[:, :, i] .+ k_trans[:, :, i + 1]);
  end for;

//Energy Balance
   //internal CVs:
  rho[2:nX-1, 2:nY-1, 2:nZ-1].*cp[2:nX-1, 2:nY-1, 2:nZ-1].*dV.*der(T[2:nX-1, 2:nY-1, 2:nZ-1]) =
                Ax./dx.*(
                kx[1:nX-2, 2:nY-1, 2:nZ-1].*(T[1:nX-2, 2:nY-1, 2:nZ-1] .- T[2:nX-1, 2:nY-1, 2:nZ-1]) .+
                kx[2:nX-1, 2:nY-1, 2:nZ-1].*(T[3:nX, 2:nY-1, 2:nZ-1] .- T[2:nX-1, 2:nY-1, 2:nZ-1]))
             .+ Ay./dy.*(
                ky[2:nX-1, 1:nY-2, 2:nZ-1].*(T[2:nX-1, 1:nY-2, 2:nZ-1] .- T[2:nX-1, 2:nY-1, 2:nZ-1]) .+
                ky[2:nX-1, 2:nY-1, 2:nZ-1].*(T[2:nX-1, 3:nY, 2:nZ-1] .- T[2:nX-1, 2:nY-1, 2:nZ-1]))
             .+ Az./dz.*(
                kz[2:nX-1, 2:nY-1, 1:nZ-2].*(T[2:nX-1, 2:nY-1, 1:nZ-2] .- T[2:nX-1, 2:nY-1, 2:nZ-1]) .+
                kz[2:nX-1, 2:nY-1, 2:nZ-1].*(T[2:nX-1, 2:nY-1, 3:nZ] .- T[2:nX-1, 2:nY-1, 2:nZ-1]));

    //corner CVs:
       //front face: x/Lx = 0
  rho[1,1,1].*cp[1,1,1].*dV.*der(T[1,1,1])= dhT2_x0.Q[1,1] .+ Ax./dx.*kx[1,1,1].*(T[2,1,1] .- T[1,1,1])
                                         .+ dhT2_y0.Q[1,1] .+ Ay./dy.*ky[1,1,1].*(T[1,2,1] .- T[1,1,1])
                                         .+ dhT2_z0.Q[1,1] .+ Az./dz.*kz[1,1,1].*(T[1,1,2] .- T[1,1,1]);

  rho[1,nY,1].*cp[1,nY,1].*dV.*der(T[1,nY,1])= dhT2_x0.Q[nY,1] .+ Ax./dx.*kx[1,nY,1]  .*(T[2,nY,1]  .- T[1,nY,1])
                                            .+ dhT2_yN.Q[1,1]  .+ Ay./dy.*ky[1,nY-1,1].*(T[1,nY-1,1].- T[1,nY,1])
                                            .+ dhT2_z0.Q[1,nY] .+ Az./dz.*kz[1,nY,1]  .*(T[1,nY,2]  .- T[1,nY,1]);

  rho[1,1,nZ].*cp[1,1,nZ].*dV.*der(T[1,1,nZ])= dhT2_x0.Q[1,nZ] .+ Ax./dx.*kx[1,1,nZ]  .*(T[2,1,nZ]  .- T[1,1,nZ])
                                            .+ dhT2_y0.Q[1,nZ] .+ Ay./dy.*ky[1,1,nZ]  .*(T[1,2,nZ]  .- T[1,1,nZ])
                                            .+ dhT2_zN.Q[1,1]  .+ Az./dz.*kz[1,1,nZ-1].*(T[1,1,nZ-1].- T[1,1,nZ]);

  rho[1,nY,nZ].*cp[1,nY,nZ].*dV.*der(T[1,nY,nZ])= dhT2_x0.Q[nY,nZ] .+ Ax./dx.*kx[1,nY,nZ]  .*(T[2,nY,nZ]  .- T[1,nY,nZ])
                                               .+ dhT2_yN.Q[1,nZ]  .+ Ay./dy.*ky[1,nY-1,nZ].*(T[1,nY-1,nZ].- T[1,nY,nZ])
                                               .+ dhT2_zN.Q[1,nY]  .+ Az./dz.*kz[1,nY,nZ-1].*(T[1,nY,nZ-1].- T[1,nY,nZ]);

        //back face: x/Lx = 1

  rho[nX,1,1].*cp[nX,1,1].*dV.*der(T[nX,1,1])= dhT2_xN.Q[1,1]  .+ Ax./dx.*kx[nX-1,1,1].*(T[nX-1,1,1].- T[nX,1,1])
                                            .+ dhT2_y0.Q[nX,1] .+ Ay./dy.*ky[nX,1,1]  .*(T[nX,2,1]  .- T[nX,1,1])
                                            .+ dhT2_z0.Q[nX,1] .+ Az./dz.*kz[nX,1,1]  .*(T[nX,1,2]  .- T[nX,1,1]);

  rho[nX,nY,1].*cp[nX,nY,1].*dV.*der(T[nX,nY,1])= dhT2_xN.Q[nY,1]  .+ Ax./dx.*kx[nX-1,nY,1].*(T[nX-1,nY,1] .- T[nX,nY,1])
                                               .+ dhT2_yN.Q[nX,1]  .+ Ay./dy.*ky[nX,nY-1,1].*(T[nX,nY-1,1] .- T[nX,nY,1])
                                               .+ dhT2_z0.Q[nX,nY] .+ Az./dz.*kz[nX,nY,1]  .*(T[nX,nY,2]   .- T[nX,nY,1]);

  rho[nX,1,nZ].*cp[nX,1,nZ].*dV.*der(T[nX,1,nZ])= dhT2_xN.Q[1,nZ]  .+ Ax./dx.*kx[nX-1,1,nZ] .*(T[nX-1,1,nZ] .- T[nX,1,nZ])
                                               .+ dhT2_y0.Q[nX,nZ] .+ Ay./dy.*ky[nX,1,nZ]   .*(T[nX,2,nZ]   .- T[nX,1,nZ])
                                               .+ dhT2_zN.Q[nX,1]  .+ Az./dz.*kz[nX,1,nZ-1].*(T[nX,1,nZ-1] .- T[nX,1,nZ]);

  rho[nX,nY,nZ].*cp[nX,nY,nZ].*dV.*der(T[nX,nY,nZ])= dhT2_xN.Q[nY,nZ] .+ Ax./dx.*kx[nX-1,nY,nZ] .*(T[nX-1,nY,nZ] .- T[nX,nY,nZ])
                                                  .+ dhT2_yN.Q[nX,nZ] .+ Ay./dy.*ky[nX,nY-1,nZ] .*(T[nX,nY-1,nZ] .- T[nX,nY,nZ])
                                                  .+ dhT2_zN.Q[nX,nY] .+ Az./dz.*kz[nX,nY,nZ-1] .*(T[nX,nY,nZ-1] .- T[nX,nY,nZ]);

   //edge CVs:
       //front face [1,:,:]: bottom, top, left, right
    rho[1,2:nY-1,1].*cp[1,2:nY-1,1].*dV.*der(T[1,2:nY-1,1])= dhT2_x0.Q[2:nY-1,1].+ Ax./dx.*kx[1,2:nY-1,1].*(T[2,2:nY-1,1].-T[1,2:nY-1,1])
                                                          .+ dhT2_z0.Q[1,2:nY-1].+ Az./dz.*kz[1,2:nY-1,1].*(T[1,2:nY-1,2].-T[1,2:nY-1,1])
                                                          .+ Ay./dy.*(
                                                              ky[1,1:nY-2,1].*(T[1,1:nY-2,1].- T[1,2:nY-1,1])
                                                           .+ ky[1,2:nY-1,1].*(T[1,3:nY,1]  .- T[1,2:nY-1,1]));

   rho[1,2:nY-1,nZ].*cp[1,2:nY-1,nZ].*dV.*der(T[1,2:nY-1,nZ])= dhT2_x0.Q[2:nY-1,nZ].+ Ax./dx.*kx[1,2:nY-1,nZ]  .*(T[2,2:nY-1,nZ]  .- T[1,2:nY-1,nZ])
                                                            .+ dhT2_zN.Q[1,2:nY-1] .+ Az./dz.*kz[1,2:nY-1,nZ-1].*(T[1,2:nY-1,nZ-1].- T[1,2:nY-1,nZ])
                                                            .+Ay./dy.*(
                                                                  ky[1,1:nY-2,nZ].*(T[1,1:nY-2,nZ].- T[1,2:nY-1,nZ])
                                                               .+ ky[1,2:nY-1,nZ].*(T[1,3:nY,nZ]  .- T[1,2:nY-1,nZ]));

   rho[1,1,2:nZ-1].*cp[1,1,2:nZ-1].*dV.*der(T[1,1,2:nZ-1])= dhT2_x0.Q[1,2:nZ-1] .+ Ax./dx.*kx[1,1,2:nZ-1].*(T[2,1,2:nZ-1] .- T[1,1,2:nZ-1])
                                                         .+ dhT2_y0.Q[1,2:nZ-1] .+ Ay./dy.*ky[1,1,2:nZ-1].*(T[1,2,2:nZ-1] .- T[1,1,2:nZ-1])
                                                         .+ Az./dz.*(
                                                         kz[1,1,1:nZ-2].*(T[1,1,1:nZ-2] .- T[1,1,2:nZ-1])
                                                      .+ kz[1,1,2:nZ-1].*(T[1,1,3:nZ]   .- T[1,1,2:nZ-1]));

   rho[1,nY,2:nZ-1].*cp[1,nY,2:nZ-1].*dV.*der(T[1,nY,2:nZ-1])= dhT2_x0.Q[nY,2:nZ-1].+ Ax./dx.*kx[1,nY,2:nZ-1].*(T[2,nY,2:nZ-1].-T[1,nY,2:nZ-1])
                                                            .+ dhT2_yN.Q[1,2:nZ-1] .+ Ay./dy.*ky[1,nY-1,2:nZ-1].*(T[1,nY-1,2:nZ-1].-T[1,nY,2:nZ-1])
                                                            .+ Az./dz.*(
                                                            kz[1,nY,1:nZ-2].*(T[1,nY,1:nZ-2] .- T[1,nY,2:nZ-1])
                                                         .+ kz[1,nY,2:nZ-1].*(T[1,nY,3:nZ]   .- T[1,nY,2:nZ-1]));
   //back face [nX,:,:]: bottom, top, left, right
   rho[nX,2:nY-1,1].*cp[nX,2:nY-1,1].*dV.*der(T[nX,2:nY-1,1])= dhT2_xN.Q[2:nY-1,1] .+ Ax./dx.*kx[nX-1,2:nY-1,1].*(T[nX-1,2:nY-1,1].-T[nX,2:nY-1,1])
                                                            .+ dhT2_z0.Q[nX,2:nY-1].+ Az./dz.*kz[nX,2:nY-1,1]  .*(T[nX,2:nY-1,2] .- T[nX,2:nY-1,1])
                                                            .+ Ay./dy.*(
                                                            ky[nX,1:nY-2,1].*(T[nX,1:nY-2,1] .- T[nX,2:nY-1,1])
                                                         .+ ky[nX,2:nY-1,1].*(T[nX,3:nY,1] .- T[nX,2:nY-1,1]));

   rho[nX,2:nY-1,nZ].*cp[nX,2:nY-1,nZ].*dV.*der(T[nX,2:nY-1,nZ])= dhT2_xN.Q[2:nY-1,nZ].+ Ax./dx.*kx[nX-1,2:nY-1,nZ].*(T[nX-1,2:nY-1,nZ].-T[nX,2:nY-1,nZ])
                                                               .+ dhT2_zN.Q[nX,2:nY-1].+ Az./dz.*kz[nX,2:nY-1,nZ-1].*(T[nX,2:nY-1,nZ-1].-T[nX,2:nY-1,nZ])
                                                               .+Ay./dy.*(
                                                               ky[nX,1:nY-2,nZ].*(T[nX,1:nY-2,nZ] .- T[nX,2:nY-1,nZ])
                                                            .+ ky[nX,2:nY-1,nZ].*(T[nX,3:nY,nZ] .- T[nX,2:nY-1,nZ]));

   rho[nX,1,2:nZ-1].*cp[nX,1,2:nZ-1].*dV.*der(T[nX,1,2:nZ-1])= dhT2_xN.Q[1,2:nZ-1] .+ Ax./dx.*kx[nX-1,1,2:nZ-1].*(T[nX-1,1,2:nZ-1].-T[nX,1,2:nZ-1])
                                                            .+ dhT2_y0.Q[nX,2:nZ-1].+ Ay./dy.*ky[nX,1,2:nZ-1] .* (T[nX,2,2:nZ-1]  .-T[nX,1,2:nZ-1])
                                                            .+Az./dz.*(
                                                            kz[nX,1,1:nZ-2].*(T[nX,1,1:nZ-2] .- T[nX,1,2:nZ-1])
                                                         .+ kz[nX,1,2:nZ-1].*(T[nX,1,3:nZ] .- T[nX,1,2:nZ-1]));

   rho[nX,nY,2:nZ-1].*cp[nX,nY,2:nZ-1].*dV.*der(T[nX,nY,2:nZ-1])= dhT2_xN.Q[nY,2:nZ-1].+ Ax./dx.*kx[nX-1,nY,2:nZ-1].*(T[nX-1,nY,2:nZ-1].-T[nX,nY,2:nZ-1])
                                                               .+ dhT2_yN.Q[nX,2:nZ-1].+ Ay./dy.*ky[nX,nY-1,2:nZ-1].*(T[nX,nY-1,2:nZ-1].-T[nX,nY,2:nZ-1])
                                                               .+Az./dz.*(
                                                               kz[nX,nY,1:nZ-2].*(T[nX,nY,1:nZ-2] .- T[nX,nY,2:nZ-1])
                                                            .+ kz[nX,nY,2:nZ-1].*(T[nX,nY,3:nZ] .- T[nX,nY,2:nZ-1]));
   //left face: bottom, top
   rho[2:nX-1,1,1].*cp[2:nX-1,1,1].*dV.*der(T[2:nX-1,1,1])= dhT2_y0.Q[2:nX-1,1].+ Ay./dy.*ky[2:nX-1,1,1].*(T[2:nX-1,2,1].-T[2:nX-1,1,1])
                                                         .+ dhT2_z0.Q[2:nX-1,1].+ Az./dz.*kz[2:nX-1,1,1].*(T[2:nX-1,1,2].-T[2:nX-1,1,1])
                                                         .+Ax./dx.*(
                                                         kx[1:nX-2,1,1].*(T[1:nX-2,1,1] .- T[2:nX-1,1,1])
                                                      .+ kx[2:nX-1,1,1].*(T[3:nX,1,1] .- T[2:nX-1,1,1]));

   rho[2:nX-1,1,nZ].*cp[2:nX-1,1,nZ].*dV.*der(T[2:nX-1,1,nZ])= dhT2_y0.Q[2:nX-1,nZ].+ Ay./dy.*ky[2:nX-1,1,nZ].*(T[2:nX-1,2,nZ].-T[2:nX-1,1,nZ])
                                                            .+ dhT2_zN.Q[2:nX-1,1].+ Az./dz.*kz[2:nX-1,1,nZ-1].*(T[2:nX-1,1,nZ-1].-T[2:nX-1,1,nZ])
                                                            .+ Ax./dx.*(
                                                            kx[1:nX-2,1,nZ].*(T[1:nX-2,1,nZ] .- T[2:nX-1,1,nZ])
                                                         .+ kx[2:nX-1,1,nZ].*(T[3:nX,1,nZ] .- T[2:nX-1,1,nZ]));
   //right face: bottom, top
   rho[2:nX-1,nY,1].*cp[2:nX-1,nY,1].*dV.*der(T[2:nX-1,nY,1])= dhT2_yN.Q[2:nX-1,1].+Ay./dy.*ky[2:nX-1,nY-1,1].*(T[2:nX-1,nY-1,1].-T[2:nX-1,nY,1])
                                                            .+ dhT2_z0.Q[2:nX-1,nY].+ Az./dz.*kz[2:nX-1,nY,1].*(T[2:nX-1,nY,2].-T[2:nX-1,nY,1])
                                                            .+ Ax./dx.*(
                                                            kx[1:nX-2,nY,1].*(T[1:nX-2,nY,1] .- T[2:nX-1,nY,1])
                                                         .+ kx[2:nX-1,nY,1].*(T[3:nX,nY,1] .- T[2:nX-1,nY,1]));

   rho[2:nX-1,nY,nZ].*cp[2:nX-1,nY,nZ].*dV.*der(T[2:nX-1,nY,nZ])= dhT2_yN.Q[2:nX-1,nZ].+ Ay./dy*ky[2:nX-1,nY-1,nZ].*(T[2:nX-1,nY-1,nZ].-T[2:nX-1,nY,nZ])
                                                               .+ dhT2_zN.Q[2:nX-1,nY].+ Az./dz.*kz[2:nX-1,nY,nZ-1].*(T[2:nX-1,nY,nZ-1].-T[2:nX-1,nY,nZ])
                                                               .+ Ax./dx.*(
                                                               kx[1:nX-2,nY,nZ].*(T[1:nX-2,nY,nZ] .- T[2:nX-1,nY,nZ])
                                                            .+ kx[2:nX-1,nY,nZ].*(T[3:nX,nY,nZ] .- T[2:nX-1,nY,nZ]));

   //center CVs: front,back, bottom,top, left,right
   rho[1,2:nY-1,2:nZ-1].*cp[1,2:nY-1,2:nZ-1].*dV.*der(T[1,2:nY-1,2:nZ-1])= dhT2_x0.Q[2:nY-1,2:nZ-1].+ Ax./dx.*kx[1,2:nY-1,2:nZ-1].*(T[2,2:nY-1,2:nZ-1].-T[1,2:nY-1,2:nZ-1])
                                                                        .+ Az./dz.*(kz[1,2:nY-1,1:nZ-2].*(T[1,2:nY-1,1:nZ-2] .- T[1,2:nY-1,2:nZ-1])
                                                                                 .+ kz[1,2:nY-1,2:nZ-1].*(T[1,2:nY-1,3:nZ] .- T[1,2:nY-1,2:nZ-1]))
                                                                        .+ Ay./dy.*(ky[1,1:nY-2,2:nZ-1].*(T[1,1:nY-2,2:nZ-1] .- T[1,2:nY-1,2:nZ-1])
                                                                                 .+ ky[1,2:nY-1,2:nZ-1].*(T[1,3:nY,2:nZ-1] .- T[1,2:nY-1,2:nZ-1]));

   rho[nX,2:nY-1,2:nZ-1].*cp[nX,2:nY-1,2:nZ-1].*dV.*der(T[nX,2:nY-1,2:nZ-1])= dhT2_xN.Q[2:nY-1,2:nZ-1].+ Ax./dx.*kx[nX-1,2:nY-1,2:nZ-1].*(T[nX-1,2:nY-1,2:nZ-1].-T[nX,2:nY-1,2:nZ-1])
                                                                           .+Az./dz.*(kz[nX,2:nY-1,1:nZ-2].*(T[nX,2:nY-1,1:nZ-2] .- T[nX,2:nY-1,2:nZ-1])
                                                                                   .+ kz[nX,2:nY-1,2:nZ-1].*(T[nX,2:nY-1,3:nZ] .- T[nX,2:nY-1,2:nZ-1]))
                                                                           .+Ay./dy.*(ky[nX,1:nY-2,2:nZ-1].*(T[nX,1:nY-2,2:nZ-1] .- T[nX,2:nY-1,2:nZ-1])
                                                                                   .+ ky[nX,2:nY-1,2:nZ-1].*(T[nX,3:nY,2:nZ-1] .- T[nX,2:nY-1,2:nZ-1]));

   rho[2:nX-1,2:nY-1,1].*cp[2:nX-1,2:nY-1,1].*dV.*der(T[2:nX-1,2:nY-1,1])= dhT2_z0.Q[2:nX-1,2:nY-1].+ Az./dz.*kz[2:nX-1,2:nY-1,1].*(T[2:nX-1,2:nY-1,2].-T[2:nX-1,2:nY-1,1])
                                                                        .+ Ax./dx.*(kx[1:nX-2,2:nY-1,1].*(T[1:nX-2,2:nY-1,1] .- T[2:nX-1,2:nY-1,1])
                                                                                 .+ kx[2:nX-1,2:nY-1,1].*(T[3:nX,2:nY-1,1].-T[2:nX-1,2:nY-1,1]))
                                                                        .+ Ay./dy.*(ky[2:nX-1,1:nY-2,1].*(T[2:nX-1,1:nY-2,1] .- T[2:nX-1,2:nY-1,1])
                                                                                 .+ ky[2:nX-1,2:nY-1,1].*(T[2:nX-1,3:nY,1] .- T[2:nX-1,2:nY-1,1]));

   rho[2:nX-1,2:nY-1,nZ].*cp[2:nX-1,2:nY-1,nZ].*dV.*der(T[2:nX-1,2:nY-1,nZ])= dhT2_zN.Q[2:nX-1,2:nY-1].+ Az./dz.*kz[2:nX-1,2:nY-1,nZ-1].*(T[2:nX-1,2:nY-1,nZ-1].-T[2:nX-1,2:nY-1,nZ])
                                                                           .+ Ax./dx.*(kx[1:nX-2,2:nY-1,nZ].*(T[1:nX-2,2:nY-1,nZ] .- T[2:nX-1,2:nY-1,nZ])
                                                                                    .+ kx[2:nX-1,2:nY-1,nZ].*(T[3:nX,2:nY-1,nZ] .- T[2:nX-1,2:nY-1,nZ]))
                                                                           .+ Ay./dy.*(ky[2:nX-1,1:nY-2,nZ].*(T[2:nX-1,1:nY-2,nZ] .- T[2:nX-1,2:nY-1,nZ])
                                                                                    .+ ky[2:nX-1,2:nY-1,nZ].*(T[2:nX-1,3:nY,nZ] .- T[2:nX-1,2:nY-1,nZ]));

   rho[2:nX-1,1,2:nZ-1].*cp[2:nX-1,1,2:nZ-1].*dV.*der(T[2:nX-1,1,2:nZ-1])= dhT2_y0.Q[2:nX-1,2:nZ-1].+ Ay./dy.*ky[2:nX-1,1,2:nZ-1].*(T[2:nX-1,2,2:nZ-1].-T[2:nX-1,1,2:nZ-1])
                                                                         .+ Ax./dx.*(kx[1:nX-2,1,2:nZ-1].*(T[1:nX-2,1,2:nZ-1] .- T[2:nX-1,1,2:nZ-1])
                                                                                  .+ kx[2:nX-1,1,2:nZ-1].*(T[3:nX,1,2:nZ-1] .- T[2:nX-1,1,2:nZ-1]))
                                                                         .+ Az./dz.*(kz[2:nX-1,1,1:nZ-2].*(T[2:nX-1,1,1:nZ-2] .- T[2:nX-1,1,2:nZ-1])
                                                                                  .+ kz[2:nX-1,1,2:nZ-1].*(T[2:nX-1,1,3:nZ] .- T[2:nX-1,1,2:nZ-1]));

   rho[2:nX-1,nY,2:nZ-1].*cp[2:nX-1,nY,2:nZ-1].*dV.*der(T[2:nX-1,nY,2:nZ-1])= dhT2_yN.Q[2:nX-1,2:nZ-1].+Ay./dy.*ky[2:nX-1,nY-1,2:nZ-1].*(T[2:nX-1,nY-1,2:nZ-1].-T[2:nX-1,nY,2:nZ-1])
                                                                            .+ Ax./dx.*(kx[1:nX-2,nY,2:nZ-1].*(T[1:nX-2,nY,2:nZ-1] .- T[2:nX-1,nY,2:nZ-1])
                                                                                     .+ kx[2:nX-1,nY,2:nZ-1].*(T[3:nX,nY,2:nZ-1] .- T[2:nX-1,nY,2:nZ-1]))
                                                                            .+ Az./dz.*(kz[2:nX-1,nY,1:nZ-2].*(T[2:nX-1,nY,1:nZ-2] .- T[2:nX-1,nY,2:nZ-1])
                                                                                     .+ kz[2:nX-1,nY,2:nZ-1].*(T[2:nX-1,nY,3:nZ] .- T[2:nX-1,nY,2:nZ-1]));

  dhT2_z0.Q = (k_trans[:,:,1] ./(0.5*dz))*Az .* (dhT2_z0.T .- T[:,:,1]);  // Conductive Flow bottom
  dhT2_zN.Q = (k_trans[:,:,nZ]./(0.5*dz))*Az .* (dhT2_zN.T .- T[:,:, nZ]); // Conductive Flow top
  dhT2_y0.Q = (k_long[:,1,:]  ./(0.5*dy))*Ay .* (dhT2_y0.T .- T[:,1,:]); // Conductive Flow Left
  dhT2_yN.Q = (k_long[:,nY,:] ./(0.5*dy))*Ay .* (dhT2_yN.T .- T[:,nY,:]); // Conductive Flow Right
  dhT2_x0.Q = (k_long[1,:,:]  ./(0.5*dx))*Ax .* (dhT2_x0.T .- T[1,:,:]); // Conductive Flow back
  dhT2_xN.Q = (k_long[nX,:,:] ./(0.5*dx))*Ax .* (dhT2_xN.T .- T[nX,:,:]); // Conductive Flow front

 annotation(Dialog(group="Custom Material Only"),
                                                Dialog(group="Custom Material Only"),
                                                             Dialog(group="Custom Material Only"),
              Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Rectangle(
          extent={{-100,60},{100,20}},
          lineColor={28,108,200},
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-80,60},{-60,38}},
          lineColor={244,125,35},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="zN"),
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
          extent={{-78,-40},{-58,-62}},
          lineColor={244,125,35},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="z0"),
        Text(
          extent={{62,60},{82,38}},
          lineColor={244,125,35},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="yN"),
        Text(
          extent={{62,-38},{82,-60}},
          lineColor={244,125,35},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="y0"),
        Text(
          extent={{-96,12},{-76,-10}},
          lineColor={244,125,35},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="x0"),
        Text(
          extent={{76,10},{96,-12}},
          lineColor={244,125,35},
          fillColor={95,95,95},
          fillPattern=FillPattern.None,
          textString="xN"),
          Text(
            extent={{-140,-18},{140,-58}},
            lineColor={255,255,255},
          textString="%name")}),                                    Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    Documentation(revisions="<html>
<ul>
<li><i>31 Aug 2021</i> by <a href=\"faisal.sedeqi@dlr.de\">Faisal Sedeqi</a>:<br>Modified to consider conduction in different planes</li>
<li><i>03 Aug 2021</i> by <a href=\"faisal.sedeqi@dlr.de\">Faisal Sedeqi</a>:<br>Removed Radiation ports. Will change conduction between CVs in next change.</li>
<li><i>19 Jul 2021</i> by <a href=\"Rene.Lorenz@dlr.de\">Rene Lorenz</a>:<br>First release. </li>
</ul>
</html>"));
end Solid3D;

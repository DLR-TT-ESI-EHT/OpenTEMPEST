within OpenTEMPEST.Heat;
model HeatTransferNodeUA
  extends OpenTEMPEST.Heat.BaseClasses.HeatTransferNodeBase;

    import SI = Modelica.SIunits;
    parameter Real UAnom=1;
    parameter SI.MassFlowRate mfNom=1;

    Real UA(start = UAnom);

  Modelica.Blocks.Interfaces.RealInput mf annotation (Placement(transformation(
        extent={{20,-20},{-20,20}},
        rotation=90,
        origin={-78,66}), iconTransformation(
        extent={{20,-20},{-20,20}},
        rotation=90,
        origin={-84,40})));
equation

  UA = UAnom*(mf/mfNom);

  hotWall.Q = UA/N*(hotWall.T .- coldWall.T);

end HeatTransferNodeUA;

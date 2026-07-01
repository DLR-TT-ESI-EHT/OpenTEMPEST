within OpenTEMPEST.Flow;
model Flow1D2w
  extends OpenTEMPEST.Flow.Flow1D(
                        Qext = wall.Q .+ wall2.Q);

  ThermoPower.Thermal.DHTVolumes wall2(N=N)
                                       annotation (Placement(
        transformation(extent={{-80,-46},{80,-26}}),  iconTransformation(extent=
           {{-80,-50},{80,-32}})));

equation

  wall2.T[:] = Gas[:].T;

end Flow1D2w;

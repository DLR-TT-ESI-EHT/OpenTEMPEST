within OpenTEMPEST.SOC.Stack;
model StackCFSimpCell
  "Stack based on the crossflow cell with 1D fuel and air channels"
  import SI = Modelica.SIunits;

  // general parameters
  parameter Integer nX=5 "Number of control volumes in first direction" annotation(Dialog(group="General Parameters"));
  parameter Integer nY=5 "Number of control volumes in second direction" annotation(Dialog(group="General Parameters"));
  parameter Integer Ncell(min=1)=5 "Number of cells in the stack" annotation(Dialog(group="General Parameters"));
  parameter Boolean LUDS=false "Set true if Linear Upwind Difference Scheme wanted (more accuray), false for Upwind Difference Scheme (more speed and stability)" annotation(Dialog(group="General Parameters"));
  parameter Boolean heatTransferCorrelationFormDuct=true "true for Nusselt correlation duct geometry with characteristic length=2*lZ (default), false for plate geometry with characteristic length=lX" annotation(Dialog(group="General Parameters"));

  // Dimensions from http://dx.doi.org/10.2139/ssrn.3987808
  parameter SI.Length lX = 0.1 "Length of solid" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lY = 0.1 "Width of solid" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lXpen = lX "Length of pen" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lYpen = lY "Width of pen" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lZpen = 3.425e-4 "Thickness of pen" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lXac = lX "Length of air channel" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lYac = lY "Width of air channel" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lZac = 1e-3 "Height of air channel" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lXfc = lX "Length of fuel channel" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lYfc = lY "Width of fuel channel" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lZfc = 1e-3 "Height of fuel channel" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length lZsolid = 0.2e-3 "Height of interconnector" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length Bac = porAC*lYac "Width of air channel without ribs" annotation(Dialog(tab="Dimensions"));
  parameter SI.Length Bfc = porFC*lYfc "WIdth of fuel channel without ribs" annotation(Dialog(tab="Dimensions"));

  // Initialization
  parameter SI.Temperature TStart=1023.15   "Uniform start temperature" annotation(Dialog(tab="Initialization"));
  parameter SI.AbsolutePressure pStart=101325   "Starting pressure" annotation(Dialog(tab="Initialization"));
  parameter SI.CurrentDensity Jstart = 0 "Starting current density in PEN" annotation(Dialog(tab="Initialization"));
  parameter SI.MassFraction xStartFC[FCMedium.nXi] = FCMedium.reference_X "Starting mass fraction in fuel channel" annotation(Dialog(tab="Initialization"));
  parameter SI.MassFraction xStartAC[ACMedium.nX] = ACMedium.reference_X "Starting mass fraction in air channel" annotation(Dialog(tab="Initialization"));

  // PEN parameters
  replaceable model Electrochem =
      OpenTEMPEST.SOC.Electrochem.Components.Crossflow_Electrochem
  constrainedby OpenTEMPEST.SOC.Electrochem.Components.ElectrochemBase         "Electrochemical model used" annotation (
  Dialog(tab="PEN"),
  Placement(transformation(extent={{78,50},{98,70}})),
  choicesAllMatching=true);
  parameter SI.ThermalConductivity kCustom_trans=2.16 "Thermal Conductivity across layers of PEN in W/mK" annotation(Dialog(tab="PEN"));
  parameter SI.ThermalConductivity kCustom_long = 2.16 "Thermal Conductivity in plane of layers of PEN in W/mK (=k_trans for homogeneous materials)" annotation(Dialog(tab="PEN"));
  parameter SI.Density rhoPEN=5900 "Density of PEN in kg/m3" annotation(Dialog(tab="PEN"));
  parameter SI.SpecificHeatCapacity cpPEN= 500 "Specific heat capacity of PEN in J/kgK" annotation(Dialog(tab="PEN"));
  parameter SI.SpectralEmissivity epsilonPEN = 0.8 "Emissivity of Anode-Electrolyte-Cathode unit" annotation(Dialog(tab="PEN"));

  // FC parameters
  replaceable package FCMedium = Medium.Fuel_CH4         annotation(Dialog(tab="Fuel Channel"));
  parameter Real porFC = 0.4 "Porosity in fuel channel" annotation(Dialog(tab="Fuel Channel"));
  parameter Real NuFCPEN = 12 "Nusselt number fuel channel on PEN side" annotation(Dialog(tab="Fuel Channel"));
  parameter Real NuFCIC = 10 "Nusselt number fuel channel on IC side" annotation(Dialog(tab="Fuel Channel"));
  parameter Real alfafc = 0.015 "Weight for 2D temperature in z-direction convection for fuel channel" annotation(Dialog(tab="Fuel Channel"));

  // AC parameters
  replaceable package ACMedium = Medium.Air_Medium         annotation(Dialog(tab="Air Channel"));
  parameter Real porAC = 0.4 "Porosity in air channel" annotation(Dialog(tab="Air Channel"));
  parameter Real pDrop(max=0.99) = 0.04 "Pressure drop as a factor of inlet pressure (between 0 and 0.99)" annotation(Dialog(tab="Air Channel"));
  parameter Real NuACPEN = 8 "Nusselt number air channel on PEN side" annotation(Dialog(tab="Air Channel"));
  parameter Real NuACIC = 7.5 "Nusselt number air channel on IC side" annotation(Dialog(tab="Air Channel"));
  parameter Real alfaac = 0.12 "Weight for 2D temperature in z-direction convection for air channel" annotation(Dialog(tab="Air Channel"));

  // IC parameters
  parameter SI.ThermalConductivity kIC = 27 "Thermal conductivity of ribs in W/mK" annotation(Dialog(tab="Interconnects"));
  parameter SI.SpecificHeatCapacity cpIC = 500 "Ribs heat capacity in W/mK" annotation(Dialog(tab="Interconnects"));
  parameter SI.Density rhoIC=8000   "Ribs denisty in kg/m3" annotation(Dialog(tab="Interconnects"));
  parameter SI.SpectralEmissivity epsilonIC = 0.1 "Emissivity of interconnects" annotation(Dialog(tab="Interconnects"));

  OpenTEMPEST.SOC.Cell.CrossFlow.CellCFSimp cell[Ncell](
    each nX=nX,
    each nY=nY,
    each lX=lX,
    each lY=lY,
    each lXpen=lXpen,
    each lYpen=lYpen,
    each lZpen=lZpen,
    each lXac=lXac,
    each lYac=lYac,
    each lZac=lZac,
    each lXfc=lXfc,
    each lYfc=lYfc,
    each lZfc=lZfc,
    each lZsolid=lZsolid,
    each TStart=TStart,
    each pStart=pStart,
    each Jstart=Jstart,
    each xStartFC=xStartFC,
    each xStartAC=xStartAC,
    redeclare model Electrochem = Electrochem,
    each kCustom_trans=kCustom_trans,
    each kCustom_long=kCustom_long,
    each rhoPEN=rhoPEN,
    each cpPEN=cpPEN,
    each epsilonPEN=epsilonPEN,
    each porFC=porFC,
    each NuFCPEN=NuFCPEN,
    each NuFCIC=NuFCIC,
    each porAC=porAC,
    each pDrop=pDrop,
    each NuACPEN=NuACPEN,
    each NuACIC=NuACIC,
    each kIC=kIC,
    each cpIC=cpIC,
    each rhoIC=rhoIC,
    each epsilonIC=epsilonIC,
    each LUDS=LUDS,
    each heatTransferCorrelationFormDuct=heatTransferCorrelationFormDuct,
    each alfafc=alfafc,
    each alfaac=alfaac)
    annotation (Placement(transformation(extent={{-28,-24},{32,32}})));
  OpenTEMPEST.Heat.Solid2D topPlate(
    redeclare package SolidMat = OpenTEMPEST.Solid.Material.Crofer22APU,
    nX=nX,
    nY=nY,
    Tstartbar=TStart,
    lX=lX,
    lY=lY,
    lZ=lZsolid)
    annotation (Placement(transformation(extent={{-16,76},{16,102}})));
  OpenTEMPEST.Heat.Solid2D bottomPlate(
    redeclare package SolidMat = OpenTEMPEST.Solid.Material.Crofer22APU,
    nX=nX,
    nY=nY,
    Tstartbar=TStart,
    lX=lX,
    lY=lY,
    lZ=lZsolid)
    annotation (Placement(transformation(extent={{-16,-102},{16,-76}})));
  Flow.Manifold manifoldFuel(redeclare package Medium = Medium.Fuel_CH4,
      nPorts_b=Ncell)
    annotation (Placement(transformation(extent={{-64,42},{-54,68}})));
  Flow.Manifold_out manifoldFuel_out(redeclare package Medium = Medium.Fuel_CH4,
      nPorts_a=Ncell)
    annotation (Placement(transformation(extent={{36,40},{64,70}})));
  Flow.Manifold manifoldAir(redeclare package Medium = Medium.Air_Medium,
      nPorts_b=Ncell)
    annotation (Placement(transformation(extent={{-64,-62},{-54,-36}})));
  Flow.Manifold_out manifoldAir_out(redeclare package Medium =
        Medium.Air_Medium, nPorts_a=Ncell)
    annotation (Placement(transformation(extent={{36,-62},{64,-32}})));
  Modelica.Electrical.Analog.Interfaces.PositivePin pin_p annotation (Placement(transformation(extent={{-110,
            -28},{-90,-8}})));
  Modelica.Electrical.Analog.Interfaces.NegativePin pin_n annotation (Placement(transformation(extent={{-110,20},
            {-90,40}})));
  ThermoPower.Gas.FlangeA fuelIn(redeclare package Medium = Medium.Fuel_CH4)
    annotation (Placement(transformation(extent={{-110,48},{-90,68}})));
  ThermoPower.Gas.FlangeA airIn(redeclare package Medium = Medium.Air_Medium)
    annotation (Placement(transformation(extent={{-110,-62},{-90,-42}})));
  ThermoPower.Gas.FlangeB fuelOut(redeclare package Medium = Medium.Fuel_CH4)
    annotation (Placement(transformation(extent={{90,30},{110,50}})));
  ThermoPower.Gas.FlangeB airout(redeclare package Medium = Medium.Air_Medium)
    annotation (Placement(transformation(extent={{90,-44},{110,-24}})));

  Heat.DHT_DHTVolumes2D dHT_DHT2_yN(i=Ncell, j=nX)
    annotation (Placement(transformation(extent={{58,80},{72,94}})));
  Heat.DHT_DHTVolumes2D dHT_DHT2_y0(i=Ncell, j=nX)
    annotation (Placement(transformation(extent={{38,-82},{52,-68}})));
  Heat.DHT_DHTVolumes2D dHT_DHT2_xN(i=Ncell, j=nY)
    annotation (Placement(transformation(extent={{52,0},{66,14}})));
  Heat.DHT_DHTVolumes2D dHT_DHT2_x0(i=Ncell, j=nY) annotation (Placement(
        transformation(
        extent={{7,-7},{-7,7}},
        rotation=0,
        origin={-63,7})));
  Heat.DHTVolumes2D dHT2_top(i=nX, j=nY)
    annotation (Placement(transformation(extent={{-84,94},{-72,106}})));
  Heat.DHTVolumes2D dHT2_bottom(i=nX, j=nY)
    annotation (Placement(transformation(extent={{-84,-104},{-72,-92}})));
  Heat.DHTVolumes2D dHT2_x0(i=Ncell, j=nY)
    annotation (Placement(transformation(extent={{-106,-6},{-94,6}})));
  Heat.DHTVolumes2D dHT2_xN(i=Ncell, j=nY)
    annotation (Placement(transformation(extent={{94,-8},{106,4}})));
  Heat.DHTVolumes2D dHT2_yN(i=Ncell, j=nX)
    annotation (Placement(transformation(extent={{94,62},{106,74}})));
  Heat.DHTVolumes2D dHT2_y0(i=Ncell, j=nX)
    annotation (Placement(transformation(extent={{94,-70},{106,-58}})));
  ThermoPower.Thermal.DHTVolumes dHTtopPlate_x0(N=nY)
    annotation (Placement(transformation(extent={{-106,76},{-94,88}})));
  ThermoPower.Thermal.DHTVolumes dHTtopPlate_xN(N=nY)
    annotation (Placement(transformation(extent={{94,82},{106,94}})));
  ThermoPower.Thermal.DHTVolumes dHTbottomPlate_xN(N=nY)
    annotation (Placement(transformation(extent={{94,-90},{106,-78}})));
  ThermoPower.Thermal.DHTVolumes dHTbottomPlate_x0(N=nY)
    annotation (Placement(transformation(extent={{-106,-86},{-94,-74}})));

  ThermoPower.Thermal.DHTVolumes dHTtopPlate_yN(N=nX)
    annotation (Placement(transformation(extent={{30,94},{42,106}})));
  ThermoPower.Thermal.DHTVolumes dHTbottomPlate_yN(N=nX)
    annotation (Placement(transformation(extent={{30,-106},{42,-94}})));
  ThermoPower.Thermal.DHTVolumes dHTtopPlate_y0(N=nX)
    annotation (Placement(transformation(extent={{-44,94},{-32,106}})));
  ThermoPower.Thermal.DHTVolumes dHTbottomPlate_y0(N=nX)
    annotation (Placement(transformation(extent={{-44,-106},{-32,-94}})));
equation

  // Electrical connections - External
  connect(pin_n, cell[1].pin_n);
  connect(pin_p, cell[Ncell].pin_p);

  // Gas flow connections
  // Air
  connect(airIn, manifoldAir.port_a);
  connect(manifoldAir.ports_b[1:Ncell], cell.airIn);
  connect(cell.airout, manifoldAir_out.ports_a[1:Ncell]);
  connect(manifoldAir_out.port_b, airout);
  // Fuel
  connect(fuelIn, manifoldFuel.port_a);
  connect(manifoldFuel.ports_b[1:Ncell], cell.fuelIn);
  connect(cell.fuelOut, manifoldFuel_out.ports_a[1:Ncell]);
  connect(manifoldFuel_out.port_b, fuelOut);

  // z-direction
  // External
  connect(bottomPlate.dhT2_z0, dHT2_bottom);
  connect(topPlate.dhT2_z1, dHT2_top);
  // Plates - Cells
  connect(topPlate.dhT2_z0, cell[Ncell].dHT2_z1);
  connect(cell[1].dHT2_z0, bottomPlate.dhT2_z1);
  // Inter cell connections
  for i in 1:Ncell-1 loop
    connect(cell[i].pin_p, cell[i + 1].pin_n);
    connect(cell[i].dHT2_z1, cell[i + 1].dHT2_z0);
  end for;

  // y-direction
  // Top and bottom plates
  connect(dHTtopPlate_y0, topPlate.dhT_y0);
  connect(dHTbottomPlate_y0, bottomPlate.dhT_y0);
  connect(dHTtopPlate_yN, topPlate.dhT_yN);
  connect(dHTbottomPlate_yN, bottomPlate.dhT_yN);
  // Cells
  for i in 1:Ncell loop
    connect(dHT_DHT2_y0.DHT_port[i], cell[i].dHT_y0);
    connect(dHT_DHT2_yN.DHT_port[i], cell[i].dHT_yN);
  end for;
  // External connections
  connect(dHT_DHT2_y0.DHT2D_port, dHT2_y0);
  connect(dHT_DHT2_yN.DHT2D_port, dHT2_yN);

  // x-direction
  // Top and bottom plates
  connect(topPlate.dhT_x0, dHTtopPlate_x0);
  connect(topPlate.dhT_xN, dHTtopPlate_xN);
  connect(bottomPlate.dhT_x0, dHTbottomPlate_x0);
  connect(bottomPlate.dhT_xN, dHTbottomPlate_xN);
  // Cells
  for i in 1:Ncell loop
    connect(dHT_DHT2_x0.DHT_port[i], cell[i].dHT_x0);
    connect(dHT_DHT2_xN.DHT_port[i], cell[i].dHT_xN);
  end for;
  // External connections
  connect(dHT_DHT2_x0.DHT2D_port, dHT2_x0);
  connect(dHT_DHT2_xN.DHT2D_port, dHT2_xN);

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={28,108,200},
          fillColor={0,127,127},
          fillPattern=FillPattern.CrossDiag,
          lineThickness=0.5)}),                                  Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    experiment(__Dymola_Algorithm="Cvode"),
    Documentation(info="<html>
<html>
<body>
<!--StartFragment-->&lt;h2&gt;SimplifiedStack1D&lt;/h2&gt;
&lt;p&gt;
This model represents a 1D stack of solid oxide fuel cells using a combination of detailed and simplified cells. 
It consists of &lt;code&gt;Ncell&lt;/code&gt; total cells, where &lt;code&gt;NdetailedCell&lt;/code&gt; are calculated in detail and the remaining are grouped into &lt;code&gt;nSimpMult&lt;/code&gt; simplified vertical blocks. 
The model includes thermal, electrical, and flow interactions for each cell, fuel and air manifolds, interconnects, and stack plates. 
Simplified cells receive averaged information from corresponding detailed cells to reduce computational cost while preserving stack-level accuracy.
&lt;/p&gt;<!--EndFragment-->
</body>
</html>
</html>", revisions="<html>
</html>"));
end StackCFSimpCell;

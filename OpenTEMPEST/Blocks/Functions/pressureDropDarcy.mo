within OpenTEMPEST.Blocks.Functions;
function pressureDropDarcy
  "Calculates the pressure drop in a porous medium using Darcy's law"
  extends Modelica.Icons.Function;

  input Modelica.SIunits.MassFlowRate mf "Mass flow rate of the fluid";
  input Modelica.SIunits.DynamicViscosity eta "Dynamic viscosity of the fluid";
  input Modelica.SIunits.Density rho "Density of the fluid";
  input Modelica.SIunits.Length L "length of flow channel";
  input Modelica.SIunits.Area A "Total cross-sectional area of flow channel";
  input Modelica.SIunits.Area k "Permeability of the porous medium in the flow channel";

  output Modelica.SIunits.PressureDifference dp "Pressure loss";

algorithm

  dp := mf*eta*L/A/rho/k;

end pressureDropDarcy;

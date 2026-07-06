within OpenTEMPEST.Blocks.ComputationBlocks;
model PlossConverter "Is needed to take the power loss in PowerElectronics into account."
  extends Modelica.Blocks.Interfaces.SISO;
  extends Modelica.Icons.UnderConstruction; // comparison on if statment needs to be chacked

  parameter Modelica.SIunits.Power Pdesign=2377 "Design power in watts (only used in old modeltype)";
  parameter Boolean AC "True if the power source is AC";
  parameter Real etaAC=0.95 "AC losses";

  parameter Enumerations.SwitchingMethod converterMode
    "Select the type of converter";

protected
  parameter Real K0_GaNFets=0.00759;
  parameter Real K1_GaNFets=-0.001153;
  parameter Real K2_GaNFets=0.004358;
  // Corrected scientific notation
  parameter Real K0_MOSFET=0.01418;
  parameter Real K1_MOSFET=-0.03306;
  parameter Real K2_MOSFET=0.045428;
  // Corrected scientific notation
  Modelica.SIunits.Power P;
  Modelica.SIunits.Power PDClossRel;
  Modelica.SIunits.Power PACloss;
  Modelica.SIunits.Power Ploss;

  //old modeltype parameters
  parameter Real a=0.048212553240941;
  parameter Real b=0.009846461934665;
  parameter Real c=0.011514564183595;
  Real x;
  //Real PlossRel;

equation
  u = P;
  if P==0 then
    y=0;
  else
    y = Ploss;
  end if;

  if AC then
        PACloss = (1 - etaAC)*P;
        x = (P - PACloss)/Pdesign;

  else  PACloss =0;
        x = P/Pdesign;
  end if;

  if converterMode == Enumerations.SwitchingMethod.old then
    PDClossRel = a*x^2 + b*x + c;
  elseif converterMode == Enumerations.SwitchingMethod.GaNFets then
          PDClossRel = K0_GaNFets*x^2 + K1_GaNFets*x + K2_GaNFets;
  else    PDClossRel = K0_MOSFET*x^2 + K1_MOSFET*x + K2_MOSFET;
  end if;

    Ploss = PACloss + PDClossRel*Pdesign;

  annotation (Documentation(info="<html>
<p>This model gives us the original power output from a DC/DC converter. we obtained a polynomial equation by implementing the curve fitting technique for the following curve obtained from the following literature: <a href=\"literature:https://doi.org/10.1016/j.ijhydene.2015.12.186\">https://doi.org/10.1016/j.ijhydene.2015.12.186</a> </p>
<p>DC/DC power conversion using GaNFets <a href=\"https://www.e3s-conferences.org/articles/e3sconf/pdf/2017/04/e3sconf_espc2017_18003.pdf\">https://www.e3s-conferences.org/articles/e3sconf/pdf/2017/04/e3sconf_espc2017_18003.pdf</a></p>
<p>DC/DC powe conversion using MoSFets <a href=\"https://dspace.nwu.ac.za/handle/10394/34177\">https://dspace.nwu.ac.za/handle/10394/34177</a></p>
<p><br>The <span style=\"font-family: monospace;\">Ploss</span> model calculates the power losses in a power converter based on the input power, design parameters, and the type of switching technology used (GaNFets or MOSFET).</p>
<table cellspacing=\"0\" cellpadding=\"2\" border=\"1\" width=\"100%\"><tr>
<td><p align=\"center\">Variable </p></td>
<td><p align=\"center\">Description</p></td>
</tr>
<tr>
<td><p align=\"center\">u</p></td>
<td><p align=\"center\">Input power to the system</p></td>
</tr>
<tr>
<td><p align=\"center\">y</p></td>
<td><p align=\"center\">caclculated total power loss</p></td>
</tr>
<tr>
<td><p align=\"center\">Pdesign </p></td>
<td><p align=\"center\">Design power</p></td>
</tr>
<tr>
<td><p align=\"center\">AC</p></td>
<td><p align=\"center\">Boolean: true if the power source is AC, false if DC</p></td>
</tr>
<tr>
<td><p align=\"center\">ConverterMode</p></td>
<td><p align=\"center\">Switching technology: GaNFets or MOSFET</p></td>
</tr>
<tr>
<td><p align=\"center\">P</p></td>
<td><p align=\"center\">input power to the model</p></td>
</tr>
<tr>
<td><p align=\"center\">PDCloss</p></td>
<td><p align=\"center\">Power loss due to DC/DC conversion</p></td>
</tr>
<tr>
<td><p align=\"center\">PACloss</p></td>
<td><p align=\"center\">Power loss due to the AC/DC conversion</p></td>
</tr>
<tr>
<td><p align=\"center\">Ploss</p></td>
<td><p align=\"center\">Total Power loss</p></td>
</tr>
</table>
</html>", revisions="<html>
<ul>
<li><i>16 April 2025</i> by <a href=\"mailto:sasikanth.vadde@dlr.de\">Sasikanth Vadde</a> and <a href=\"mailto:hans.wiggenhauser@dlr.de\">Hans Wiggenhauser</a>:<br> New version with conversion type specific power loss and separate AC loss. </li>
<li><i>01 January 2020</i> by Marius Tomberg</a> First release. </li>
</ul>
</html>"));
end PlossConverter;

within OpenTEMPEST.Blocks.ComputationBlocks;
model PlossConverter "Is needed to take the power loss in PowerElectronics into account."
  extends Modelica.Blocks.Interfaces.SISO;
  extends Modelica.Icons.UnderConstruction; // comparison on if statment needs to be chacked

  parameter Modelica.SIunits.Power Pdesign=2377 "Design power in watts (only used in old modeltype)";
  parameter Boolean constACLoss=true
                                "Rectifier/inverter loss constant if True, else pAC parameters used for 2nd order polynomial loss function. No loss for etaAC=1";
  parameter Real etaAC=0.95 "If constACLoss=true. Constant rectifier/inverter loss";
  parameter Real pAC[3] = {4.50434939e-03, 5.73763570e-03, 7.53027504e-05}  "If constACLoss=false. p[1]+x² + p[2]+x +p[3]. Default from fit to TNPC L5 IGBT curve in Fig 10 in DOI: 10.30420/566262171";

  parameter TEMPEST.Blocks.Enumerations.SwitchingMethod converterType "Select the type of converter";

// protected
// GaNFET parameters fitted to data in Fig 7 of DOI: 10.1051/e3sconf/20171618003:
  parameter Real pGan_pos[3] = {0.00929225, 0.00116731, 0.00358863} "If converterType==GaNFets. p[1]+x² + p[2]+x +p[3] for pos. power range";  // {0.00759, -0.001153, 0.004358}
  parameter Real pGan_neg[3] = {-0.01151183, -0.00069356, -0.00336578} "If converterType==GaNFets. p[1]+x² + p[2]+x +p[3] for neg. power range";

// MOSFET parameters fitted to data in https://repository.nwu.ac.za/items/19b3358c-a14d-4a30-856d-ac1472d93448:
  parameter Real pMOS_pos[3] = {0.01418, -0.03306, 0.045428} "If converterType==MOSFET. p[1]+x² + p[2]+x +p[3] for pos. power range";
  parameter Real pMOS_neg[3] = pMOS_pos                      "If converterType==MOSFET. p[1]+x² + p[2]+x +p[3] for neg. power range";

  // Corrected scientific notation
  Modelica.SIunits.Power P;
  Modelica.SIunits.Power PDClossRel;
  Modelica.SIunits.Power PACloss;
  Modelica.SIunits.Power Ploss;

  //old modeltype parameters
//   parameter Real a=0.048212553240941;
//   parameter Real b=0.009846461934665;
//   parameter Real c=0.011514564183595;
  Real x;
  //Real PlossRel;

equation
  u = P;
  if P==0 then
    y=0;
  else
    y = Ploss;
  end if;

  if constACLoss==true then
        PACloss = (1 - etaAC)*P;
        x = (P - PACloss)/Pdesign;

  else
        PACloss =pAC[1]*(P/Pdesign)^2 + pAC[2]*(P/Pdesign) + pAC[3];
        x = (P - PACloss)/Pdesign;
  end if;

  if converterType == TEMPEST.Blocks.Enumerations.SwitchingMethod.old then
      PDClossRel = 0.04821255*x^2 + 0.00984646*x + 0.01151456;
  elseif converterType == TEMPEST.Blocks.Enumerations.SwitchingMethod.GaNFets then
      if P >=0 then
          PDClossRel = pGan_pos[1]*x^2 + pGan_pos[2]*x + pGan_pos[3];
      else
          PDClossRel = pGan_neg[1]*x^2 + pGan_neg[2]*x + pGan_neg[3];
      end if;
  else
      if P >=0 then
          PDClossRel = pMOS_pos[1]*x^2 + pMOS_pos[2]*x + pMOS_pos[3];
      else
          PDClossRel = pMOS_neg[1]*x^2 + pMOS_neg[2]*x + pMOS_neg[3];
      end if;
  end if;

    Ploss = PACloss + PDClossRel*Pdesign;

  annotation (Documentation(info="<html>
<p>This model gives us the original power output from a DC/DC converter. we obtained a polynomial equation by implementing the curve fitting technique for the following curve obtained from the following literature: <a href=\"literature:https://doi.org/10.1016/j.ijhydene.2015.12.186\">https://doi.org/10.1016/j.ijhydene.2015.12.186</a> </p>
<p>DC/DC power conversion using GaNFets <a href=\"https://www.e3s-conferences.org/articles/e3sconf/pdf/2017/04/e3sconf_espc2017_18003.pdf\">https://www.e3s-conferences.org/articles/e3sconf/pdf/2017/04/e3sconf_espc2017_18003.pdf</a></p>
<p>DC/DC power conversion using MoSFets <a href=\"https://dspace.nwu.ac.za/handle/10394/34177\">https://dspace.nwu.ac.za/handle/10394/34177</a></p>
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
<td><p align=\"center\">converterType</p></td>
<td><p align=\"center\">Switching technology: GaNFets or MOSFET</p></td>
</tr>
<tr>
<td><p align=\"center\">P</p></td>
<td><p align=\"center\">input power to the model</p></td>
</tr>
<tr>
<td><p align=\"center\">PDCloss</p></td>
<td><p align=\"center\">Power loss of DC/DC conversion</p></td>
</tr>
<tr>
<td><p align=\"center\">PACloss</p></td>
<td><p align=\"center\">Power loss of the AC/DC conversion</p></td>
</tr>
<tr>
<td><p align=\"center\">Ploss</p></td>
<td><p align=\"center\">Total Power loss</p></td>
</tr>
</table>
</html>", revisions="<html>
</html>"));
end PlossConverter;

within OpenTEMPEST.Solid;
record solidMaterialsRecord
  "record holding sets of coefficients of selected solid materials for calculating material parameters"

  extends Modelica.Icons.Record;
  String name "Name of solid material";
  Real cRho[4] "kg/m³";
  Real cK_trans[4] "W/m*K; effective transversal heat conductivity (in series)";
  Real cK_long[4] "W/m*K; effective longitudinal heat conductivity (in parallel)";
  Real cCp[4] "J/kg*K";

  annotation (Documentation(info="<html>
<p>
This data record contains the coefficients for calculating rho, lambda and cp as a plynomial function of temperature:
</p>
<p>
The equations have the structure of 3rd order poynomial:
<p>
c3*T^3+ c2*T^2 + c1*T +c0
<p>
where c0, c1, c2 and c3 are the coefficients and T the temperature.
</p>

</html>"));

end solidMaterialsRecord;

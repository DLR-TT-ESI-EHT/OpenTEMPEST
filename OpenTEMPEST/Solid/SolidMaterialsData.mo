within OpenTEMPEST.Solid;
package SolidMaterialsData
  "package of solid materials, where each holds coefficients for determining its temperature dependant rho, k and cp"
  extends Modelica.Icons.Package;

   constant OpenTEMPEST.Solid.solidMaterialsRecord Custom(
    name="Custom",
    cRho=zeros(4),
    cK_trans=zeros(4),
    cK_long=zeros(4),
    cCp=zeros(4));

////steel and alloys:
  constant OpenTEMPEST.Solid.solidMaterialsRecord Steel(
    name="Steel",
    cRho={0,0,0,7900},
    cK_trans={0,0,0,21},
    cK_long={0,0,0,21},
    cCp={0,0,0,650});

  constant OpenTEMPEST.Solid.solidMaterialsRecord Crofer22APU(
    name="Crofer22APU",
    cRho={0,0,0,7700},
    cK_trans={-2.9149e-09,2.3268e-05,-2.9350e-02,32.564},
    cK_long={-2.9149e-09,2.3268e-05,-2.9350e-02,32.564},
    cCp={-1.905e-07,1.486e-04,0.351,353});

  annotation (Documentation(info="<html> Material properties for solid objects </html>"));
end SolidMaterialsData;

within OpenTEMPEST.Solid;
package Material "Solid materials ready to use as SolidMat in solid models"
  package Steel "Steel"
    extends SolidMatBase(data=OpenTEMPEST.Solid.SolidMaterialsData.Steel);
  end Steel;

 package Crofer22APU "Crofer 22 APU ferritic interconnect steel"
   extends SolidMatBase(data=OpenTEMPEST.Solid.SolidMaterialsData.Crofer22APU);
 end Crofer22APU;

  package Custom "Custom material with user defined constant properties"
    extends SolidMatBase(data=OpenTEMPEST.Solid.SolidMaterialsData.Custom);
  end Custom;
end Material;

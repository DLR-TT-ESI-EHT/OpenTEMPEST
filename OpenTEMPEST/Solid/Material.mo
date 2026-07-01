within OpenTEMPEST.Solid;
package Material
  package Steel
    extends SolidMatBase(data=OpenTEMPEST.Solid.SolidMaterialsData.Steel);
  end Steel;

 package Crofer22APU
   extends SolidMatBase(data=OpenTEMPEST.Solid.SolidMaterialsData.Crofer22APU);
 end Crofer22APU;

  package Custom
    extends SolidMatBase(data=OpenTEMPEST.Solid.SolidMaterialsData.Custom);
  end Custom;
end Material;

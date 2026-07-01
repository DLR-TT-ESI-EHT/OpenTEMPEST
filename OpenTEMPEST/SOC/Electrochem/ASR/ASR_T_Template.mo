within OpenTEMPEST.SOC.Electrochem.ASR;
model ASR_T_Template "Template: see instructions inside"

  /*
  This template can be used to create new ASR models, the base level only assumes
  a relation with Temperature, but relations with pressure and other variables 
  can also be added 
  
  
  All you need to do in this level is:
  0. Uncomment the extends section
  1. declare your parameters and variables (i.e. parameter Real A = 2)
  2. write the equation for your ASR taking note that it is discretised (i.e. ASR[1:N] = 1/A*Tpen[1:N];)
  3. ensure your parameters are passed down at the PEN/cell level
  
  EXAMPLE: 
  
  extends ASR_Base;
  
  parameter Real A;
  
  equation
      ASR = 1/A*Tpen;
  
  */

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end ASR_T_Template;

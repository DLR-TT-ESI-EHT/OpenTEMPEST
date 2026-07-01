within OpenTEMPEST.SOC.Electrochem.Components;
partial model BV_ElectrochemBase
  "Air and Ohmic Overvoltages calculated here. Voltage Balance equation also here. Fuel side overvoltages calculated at level above as well as q_electroChem"
  extends OpenTEMPEST.SOC.Electrochem.Components.ElectrochemBase;
  import SI = Modelica.SIunits;

  replaceable model ASRobj = ASR.ASR_Ohm_ESC (
      A=35.316e4,
      B=-1.264e-2,
      Yo=27.266e-2) constrainedby ASR.ASR_OhmBase                           annotation(choicesAllMatching=true);

  ASRobj asrobj(Tpen=Tpen);

  // Geometric parameters
  parameter SI.Length tauFE_Measure = 30e-6 "Measured fuel electrode thickness in m",
                      tauAE_Measure = 55e-6 "Measured air electrode thickness in m",
                      tauEl_Measure = 90e-6 "Measured electrolyte thickness in m";

  parameter SI.Length tauFE=30e-6 "Fuel electrode thickness in m",
                      tauAE=55e-6 "Air electrode thickness in m",
                      tauEl=90e-6 "Electrolyte thickness in m";

  parameter Real ASROhmDegrationRate=0 "Degradation rate of the ohmic resitance in Ohm.cm²/kh";
  parameter Real ROhmFac=1 "Constant factor for ohmic resistance (e. g. for sensitivity analysis)";

  parameter Real alphaFE=0.5 "Transfer coefficient for fuel electrode H2/H2O";
  parameter Real alphaAE=0.5 "Transfer coefficient for air electrode";
  parameter SI.Length dpFE=4e-6 "Pore diameter of porous fuel electrode",
                     dpAE=4e-6 "Pore diameter of porous air electrode";
  parameter Real PhsiFE=3.39845 "Fitting parameter for fuel electrode - effective tortousity";
  parameter Real PhsiAE=7.51231 "Fitting parameter for air electrode - effective tortousity";
  parameter Real epsiFE=0.4 "Porosity of fuel electrode";
  parameter Real epsiAE=0.4 "Porosity of air electrode";

  parameter Boolean diffusionActive=true;
  parameter Real Xstart[MediumF.nX] = {0.1, 1e-5, 1e-5, 1e-5, 0.9, 1e-5};

  replaceable Real E_ae= 106810 "Activation energy air electrode";
  replaceable Real d=0.298 "Exponent for air electrode Jo calculation";
  replaceable Real gammaAE = 2.44e6 "Pre-exponential factor";


  Real ROhm "Ohmic resistance";

  SI.Voltage Uideal "Ideal voltage";
  SI.Voltage UactFE(start=-0.02) "Fuel electrode activation overpotential";
  SI.Voltage UDiffFE "Fuel electrode diffusion overpotential";
  SI.Voltage Uohm "Ohmic overpotential";
  SI.Voltage UactAE(start=0.05) "Air electrode activation overpotential";
  SI.Voltage UDiffAE "Air electrode diffusion overpotential";

  //SI.SpecificEnthalpy dHr;
  SI.SpecificGibbsFreeEnergy delGr "Specific gibbs free energy change of reaction (defined in level above)";
  SI.SpecificGibbsFreeEnergy go_T "Gibbs Energy change with temperature (defined in level above)";
  Real yF_tpb[nSpeciesF](start=Xstart) "Initial molar composition at the fuel electrode triple-phase boundary";
  Real yA_tpb[nSpeciesA](start={0.21, 0.79}) "Initial molar composition at the air electrode triple-phase boundary";

  Real sigma_a "Conductivity of fuel electrode",
       sigma_e "Conductivity of electrolyte",
       sigma_c "Conductivity of cathode";

  SI.CurrentDensity Jo_AE "Air electrode exchange current density";

 // Diffusivities - Knudssen
   SI.DiffusionCoefficient Deff_h2_k "Knudsen diff coeff of h2";
   SI.DiffusionCoefficient Deff_h2o_k "Knudsen diff coeff of h2o";
   SI.DiffusionCoefficient Deff_co2_k "Knudsen diff coeff of h2o";
   SI.DiffusionCoefficient Deff_co_k "Knudsen diff coeff of h2o";
   SI.DiffusionCoefficient Deff_o2_k "Knudsen diff coeff of o2";
   // Diffusivities - Binary FE H reaction
   SI.DiffusionCoefficient Deff_h2n2 "Binary diff coeff of h2-n2";
   SI.DiffusionCoefficient Deff_h2co "Binary diff coeff of h2-co";
   SI.DiffusionCoefficient Deff_h2co2 "Binary diff coeff of h2-co2";
   SI.DiffusionCoefficient Deff_h2ch4 "Binary diff coeff of h2-ch4";
   SI.DiffusionCoefficient Deff_h2h2o "Binary diff coeff of h2-h2o";
   SI.DiffusionCoefficient Deff_h2on2 "Binary diff coeff of h2o-n2";
   SI.DiffusionCoefficient Deff_h2oco "Binary diff coeff of h2o-co";
   SI.DiffusionCoefficient Deff_h2oco2 "Binary diff coeff of h2o-co2";
   SI.DiffusionCoefficient Deff_h2och4 "Binary diff coeff of h2o-ch4";
   // Diffusivites - Binary FE C reaction
   SI.DiffusionCoefficient Deff_coco2 "Binary diff coeff of CO-CO2";
   SI.DiffusionCoefficient Deff_con2 "Binary diff coeff of CO-N2";
   SI.DiffusionCoefficient Deff_coch4 "Binary diff coeff of CO-CH4";
   SI.DiffusionCoefficient Deff_co2n2 "Binary diff coeff of CO2-N2";
   SI.DiffusionCoefficient Deff_co2ch4 "Binary diff coeff of CO2-CH4";
   // Diffusivities - Binary AE
   SI.DiffusionCoefficient Deff_o2n2 "Binary diff coeff of o2-n2";

equation

  // Voltage balance
  Uop = Uideal - Uohm - UactFE - UactAE - UDiffAE - UDiffFE;

  // Nernst
  Uideal = -delGr/(2*Modelica.Constants.F);

  //Calculating Ohmic Voltage Drop
  ROhm = ROhmFac*(asrobj.ASR -((tauFE_Measure/sigma_a) + (tauAE_Measure/sigma_c) + (tauEl_Measure/sigma_e)) + ((tauFE/sigma_a) + (tauAE/sigma_c) + (tauEl/sigma_e)))*(1+ASROhmDegrationRate*time/1000/3600); // /100^2

  sigma_e = 5.15e7 / Tpen * Modelica.Math.exp(-10300 / Tpen);
  sigma_c = 42e6 / Tpen * Modelica.Math.exp(-1200 / Tpen);
  sigma_a = 95e6 / Tpen * Modelica.Math.exp(-1150 / Tpen);

  Uohm = J*ROhm;

  r = J/2/Modelica.Constants.F;

  // Air Activation kinetics
  J = Jo_AE*(Modelica.Math.exp(alphaAE * 2 * Modelica.Constants.F / Modelica.Constants.R / Tpen*UactAE) - Modelica.Math.exp(-(1 - alphaAE) * 2 * Modelica.Constants.F / Modelica.Constants.R /Tpen*UactAE));
  // Exchange current density
  Jo_AE = gammaAE * Tpen * ((P_A/p0*abs(yA_tpb[1]))^d)* Modelica.Math.exp(-E_ae / (Modelica.Constants.R * Tpen));


  // Diffusivities - for calculating Ohmic overvoltages
  if diffusionActive then
   Deff_h2_k  =(dpFE/3)*smooth(1, ThermoPower.Functions.sqrtReg((8*Modelica.Constants.R
      *Tpen)/(Modelica.Constants.pi*Medium.Fuel_CH4.MMX[1])))*(epsiFE/PhsiFE);
   Deff_h2o_k =(dpFE/3)*smooth(1, ThermoPower.Functions.sqrtReg((8*Modelica.Constants.R
      *Tpen)/(Modelica.Constants.pi*Medium.Fuel_CH4.MMX[5])))*(epsiFE/PhsiFE);
   Deff_co_k  =(dpFE/3)*smooth(1, ThermoPower.Functions.sqrtReg((8*Modelica.Constants.R
      *Tpen)/(Modelica.Constants.pi*Medium.Fuel_CH4.MMX[4])))*(epsiFE/PhsiFE);
   Deff_co2_k =(dpFE/3)*smooth(1, ThermoPower.Functions.sqrtReg((8*Modelica.Constants.R
      *Tpen)/(Modelica.Constants.pi*Medium.Fuel_CH4.MMX[3])))*(epsiFE/PhsiFE);

       // Effective binary diffusion - FE - H reaction
    Deff_h2n2   = 0.8114e-4*((Tpen/310)^1.75)*(p0/P_F)*(epsiFE/PhsiFE);
    Deff_h2h2o  = 0.8684e-4*((Tpen/310)^1.75)*(p0/P_F)*(epsiFE/PhsiFE);
    Deff_h2co   = 0.8224e-4*((Tpen/310)^1.75)*(p0/P_F)*(epsiFE/PhsiFE);
    Deff_h2co2  = 0.6769e-4*((Tpen/310)^1.75)*(p0/P_F)*(epsiFE/PhsiFE);
    Deff_h2ch4  = 0.7722e-4*((Tpen/310)^1.75)*(p0/P_F)*(epsiFE/PhsiFE);
    Deff_h2on2  = 0.2231e-4*((Tpen/310)^1.75)*(p0/P_F)*(epsiFE/PhsiFE);
    Deff_h2oco  = 0.2233e-4*((Tpen/310)^1.75)*(p0/P_F)*(epsiFE/PhsiFE);
    Deff_h2oco2 = 0.1588e-4*((Tpen/310)^1.75)*(p0/P_F)*(epsiFE/PhsiFE);
    Deff_h2och4 = 0.2233e-4*((Tpen/310)^1.75)*(p0/P_F)*(epsiFE/PhsiFE);
       // Effective binary diffusion - FE - C reaction
    Deff_coco2  = 0.1622e-4*((Tpen/310)^1.75)*(p0/P_F)*(epsiFE/PhsiFE);
    Deff_con2   = 0.2159e-4*((Tpen/310)^1.75)*(p0/P_F)*(epsiFE/PhsiFE);
    Deff_coch4  = 0.2325e-4*((Tpen/310)^1.75)*(p0/P_F)*(epsiFE/PhsiFE);
    Deff_co2n2  = 0.1615e-4*((Tpen/310)^1.75)*(p0/P_F)*(epsiFE/PhsiFE);
    Deff_co2ch4 = 0.1777e-4*((Tpen/310)^1.75)*(p0/P_F)*(epsiFE/PhsiFE);
       // Effective Knudsen Diffusion - AE
    Deff_o2_k =(dpAE/3)*smooth(1, ThermoPower.Functions.sqrtReg((8*Modelica.Constants.R
      *Tpen)/(Modelica.Constants.pi*Medium.Air_Medium.MMX[1])))*(epsiAE/PhsiAE);
       // Effective binary diffusion - AE
    Deff_o2n2 = 0.2176e-4*((Tpen/310)^1.75)*(p0/P_A)*(epsiAE/PhsiAE);

     // Air Electrode TPB Composition
     yA_tpb[1]  = yA[1] - J*Modelica.Constants.R*Tpen / (4*Modelica.Constants.F*P_A) / Deff_o2_k*tauAE - yA[2]*(Modelica.Math.exp(J*Modelica.Constants.R*Tpen / (4*Modelica.Constants.F*P_A) / Deff_o2n2*tauAE) - 1);
     yA_tpb[2] = 1-yA_tpb[1];
     // Diffusion Overvoltage
     UDiffAE = -Modelica.Constants.R * Tpen / (4 * Modelica.Constants.F) * (Modelica.Math.log(abs(yA_tpb[1])) - Modelica.Math.log(abs(yA[1])));
  else
    Deff_h2_k  =-1;
    Deff_h2o_k =-1;
    Deff_co_k  =-1;
    Deff_co2_k =-1;
       // Effective binary diffusion - FE - H reaction
    Deff_h2n2   = -1;
    Deff_h2h2o  = -1;
    Deff_h2co   = -1;
    Deff_h2co2  = -1;
    Deff_h2ch4  = -1;
    Deff_h2on2  = -1;
    Deff_h2oco  = -1;
    Deff_h2oco2 = -1;
    Deff_h2och4 = -1;
       // Effective binary diffusion - FE - C reaction
    Deff_coco2  = -1;
    Deff_con2   = -1;
    Deff_coch4  = -1;
    Deff_co2n2  = -1;
    Deff_co2ch4 = -1;
       // Effective Knudsen Diffusion - AE
    Deff_o2_k = -1;
       // Effective binary diffusion - AE
    Deff_o2n2 = -1;
      // Air Electrode TPB Composition
    yA_tpb = yA;
    // Diffusion Overvoltage
    UDiffAE = 0;
  end if;

  annotation (Documentation(info="<html>
  <body>
    <p>References for parameters:</p>
    <ul>
      <li>DOI: 10.1149/1945-7111/ab6820</li>
      <li>DOI: 10.1149/1945-7111/ad5e01</li>
    </ul>
  </body>
</html>"));
end BV_ElectrochemBase;

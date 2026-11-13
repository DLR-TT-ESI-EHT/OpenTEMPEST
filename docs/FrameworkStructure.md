
# Framework Structure & Modeling Levels

openTEMPEST is built as a **modular, object-oriented framework** designed to simulate electrochemical energy systems across multiple scales. Depending on the specific research question, it enables proven simplification approaches at all modeling levels. 

<div align="center">
  <img src="Overview_Cell_Stack_Module.svg" alt="Overview of the modelling levels: cell, stack and module" width="1000">
</div>

In essence, the modeling levels are defined as follows:
1. Electrochemical objects
2. Cell models
3. Stack models
4. Multi-stack / module models
5. System models

In the following, the framework is explained taking the solid oxide cell technology as an example.




## 1. Electrochemical Objects
 
*Package: [`OpenTEMPEST.SOC.Electrochem`](library-reference/OpenTEMPEST_SOC_Electrochem.md)*

Reusable components that describe the electrochemical characteristics independent of geometry and usually obtained from cell experiments. Depending on cell technology and available measurement data, detailed or lumped approaches are available:

- **Lumped Electrochemistry (ASR)**: simplifies the U,i-characteristics as a temperature dependent area-specific resistance (ASR)  
[`OpenTEMPEST.SOC.Electrochem.Components.ASR_ElectrochemBase`](library-reference/OpenTEMPEST_SOC_Electrochem_Components.md#asr_electrochembase)

- **Detailed Electrochemistry**: breakdown of voltage losses into:
  - Ohmic overpotential (by means of an ohmic ASR (*T*) expression) 
  - Activation overpotentials (by means of Butler Volmer equations)
  - Diffusion overpotentials (by means of Dusty Gas Model)
  
  Models extending from [`OpenTEMPEST.SOC.Electrochem.Components.BV_ElectrochemBase`](library-reference/OpenTEMPEST_SOC_Electrochem_Components.md#bv_electrochembase) consider this.

By taking advantage of the OOP capabilities, we are able to abstract away different mathematical formulations for the different resistances in the electrodes. A user made electrochemistry object could be made to have different functionality and be utilised in our, or their own cell and stack models, as long as they extend from [`ElectrochemBase`](library-reference/OpenTEMPEST_SOC_Electrochem_Components.md#electrochembase). See [Change the electrochemistry](UserWorkflows.md#5-change-the-electrochemistry-replaceable-models) for how to choose or write one.

## 2. Cell Models (Repeating Units (RU))



<div align="center">
  <img src="SRU_Overview.svg" alt="Single repeating unit of a planar cell" width="700">
</div>

*Package: [`OpenTEMPEST.SOC.Cell`](library-reference/OpenTEMPEST_SOC_Cell.md)*

Cell models are composed of the following sub-components, which form a repeating unit for stack modeling:
- **Electrochemically active solid** (positive-electrolyte-negative, PEN), defined by its geometry and thermophysical properties - mainly solves the energy balance in the cell solid, incorporating the electrochemical objects from Modeling Level 1 and thermally interfacing with the channels and interconnects via convection and radiation, respectively.
- **Air and Fuel gas channels** defined by its geometry and thermophysical properties - solves the energy, mass and momentum balance equations. Interacts with the PEN and Interconnects via convection, considers side chemical reactions that may occur on the catalytic surfaces in the channels. Discretised using FVM, replaceable discretisation interpolation schemes found in [`OpenTEMPEST.Flow.FluxInterpolators`](library-reference/OpenTEMPEST_Flow_FluxInterpolators.md) are used to resolve the balance equations.
- **Interconnector**, virtually distributed into a half air electrode interconnect (AIC) and fuel electrode interconnect (FIC) for each repeating unit. Interacts with gas channels via convection and adjacent half interconnect via conduction and PEN via radiation. It extends from the `OpenTEMPEST.Heat.SolidXD` models (i.e. Solid 1D or Solid2D etc.), just like the PEN model, which holds the energy balance equations.


As of now, the following planar cell models are available:

| Model | Discretization               | Use Case                                                                                                                                                                                                                             |
| -------| ------------------------------| --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1D ([`Cell1D.Cell`](library-reference/OpenTEMPEST_SOC_Cell_Cell1D.md#cell))    | Along-the-flow (co-flow)     | Precise quantification of property distribution (temperature, current density etc.) as well as spatial and temporal gradients, observation of superimposed temperature-dependent phenomena (internal reforming/cracking/methanation) |
| 2D ([`CrossFlow.Cell`](library-reference/OpenTEMPEST_SOC_Cell_CrossFlow.md#cell))    | Along-the-flows (cross-flow) | Precise quantification of property distribution as well as spatial and temporal gradients, especially relevant when flow directions are not parallel or lateral plane phenomena becomes significant                                  |


These models can be written in many different ways. Simplification approaches can also be employed at the cell level, for example the model [`OpenTEMPEST.SOC.Cell.CrossFlow.CellCFSimp`](library-reference/OpenTEMPEST_SOC_Cell_CrossFlow.md#cellcfsimp) reduces the flow channels in the crossflow cell from 2D to 1D while keeping the solid surfaces as 2D. This has the benefit of significantly reducing the number of equations needed to solve the balance equations, as the number of computationally heavy fluid states is reduced, while keeping thermal information in the PEN to a high precision - this does come at the cost of some model accuracy. [`OpenTEMPEST.SOC.Cell.CrossFlow.CellGeneral`](library-reference/OpenTEMPEST_SOC_Cell_CrossFlow.md#cellgeneral) gives a clean interface to have functionality of both models based on the parameters set.

These component and cell models can be used to make different cell types, such as the arrayed (segmented) cells, which will be included very soon to the OpenTEMPEST library, as well as different flow configurations. Even cells with different flow geometries such as tubular cells can be modelled in this way.


<div align="center">
  <img src="Cell_FlowConfigurations.svg" alt="Cell flow configurations" width="500">
</div>


## 3. Stack Models

<div align="center">
  <img src="Stack_Overview.svg" alt="Stack model overview" width="500">
</div>

Stack models combine a set of cell models, along with end plates, manifolds, flow distribution and specifications of heat transfer to surrounding.



*Package: [`OpenTEMPEST.SOC.Stack`](library-reference/OpenTEMPEST_SOC_Stack.md)*

| Model | Discretization | Use case |
|-------|----------------|----------|
| [`Stack0D`](library-reference/OpenTEMPEST_SOC_Stack.md#stack0d) | Lumped (0D) stack temperature, ASR electrochemistry, full internal reforming | System studies, control design and long simulations where only the overall stack behaviour matters |
| [`Stack1D`](library-reference/OpenTEMPEST_SOC_Stack.md#stack1d) | `Ncell` 1D co-flow cells with end plates | Temperature and current density distribution along the flow and over the stack height |
| [`SimplifiedStack1D`](library-reference/OpenTEMPEST_SOC_Stack.md#simplifiedstack1d) | A few detailed 1D cells, the other cells lumped into simplified blocks | Large co-flow stacks at a fraction of the computational cost |
| [`StackCF`](library-reference/OpenTEMPEST_SOC_Stack.md#stackcf) | `Ncell` 2D cross-flow cells with end plates | Cross-flow stacks, 2D temperature fields in the cell plane |
| [`StackCFSimpCell`](library-reference/OpenTEMPEST_SOC_Stack.md#stackcfsimpcell) | Cross-flow cells with 1D channels and 2D solids | Cross-flow stacks with fewer equations |
| [`SimplifiedStackCF`](library-reference/OpenTEMPEST_SOC_Stack.md#simplifiedstackcf) | A few detailed cross-flow cells, the other cells lumped into simplified blocks | Large cross-flow stacks |

The simplified stacks follow the approach of M. Tomberg et al. 2022 ([doi:10.1149/1945-7111/ac7009](https://doi.org/10.1149/1945-7111/ac7009)): only a few cells are resolved in detail and the cells in between are represented by simplified blocks that interpolate the flows and keep the thermal behaviour of the stack. See [Work with stack models](UserWorkflows.md#6-work-with-stack-models) for how to use them.



## 4. Multi-Stack (Module) Models

Module models enable the depiction of commercial-like assemblies of multiple stacks inside a shared thermal environment. 

We are still working on providing an example of such a model in OpenTEMPEST - please bear with us :)

<!--- 
STILL MISSING: 
- TABLE OF AVAILABLE MODELS (OR EXAMPLES)
- REFERENCES TO MODELICA FILES

--->

## 5. System Models

**Balance of Plant (BoP)** component models enable the full depiction of real test rigs or large-scale commercial systems. BoP models are implemented to depict their off-design behavior in order to correctly reproduce the thermal system behavior in part-load conditions.

**Available BoP components** (package [`OpenTEMPEST.BOP`](library-reference/OpenTEMPEST_BOP.md)):

| Component | Models |
|---|---|
| Gas-gas heat exchangers | [`GasConditioning.BlackBoxHexGeneric`](library-reference/OpenTEMPEST_BOP_GasConditioning.md#blackboxhexgeneric) with effectiveness models in [`HexEtaTypes`](library-reference/OpenTEMPEST_BOP_GasConditioning_HexEtaTypes.md) |
| Electrical heaters | [`GasConditioning.ElectricHeater0D`](library-reference/OpenTEMPEST_BOP_GasConditioning.md#electricheater0d), [`ElectricHeater1D`](library-reference/OpenTEMPEST_BOP_GasConditioning.md#electricheater1d), [`ElectricalHeater0D`](library-reference/OpenTEMPEST_BOP_GasConditioning.md#electricalheater0d) |
| Afterburners | [`Reactors.Combustor0D`](library-reference/OpenTEMPEST_BOP_Reactors.md#combustor0d) |
| Evaporators | [`Vaporizer`](library-reference/OpenTEMPEST_BOP.md#vaporizer), [`Evaporator`](library-reference/OpenTEMPEST_BOP.md#evaporator) |
| Compressors and blowers | [`Compressor`](library-reference/OpenTEMPEST_BOP.md#compressor), [`FanSimple`](library-reference/OpenTEMPEST_BOP.md#fansimple) |
| Pre-reformers | [`Reactors.Reformer1D`](library-reference/OpenTEMPEST_BOP_Reactors.md#reformer1d), [`Reactors.Reformer1D_C3`](library-reference/OpenTEMPEST_BOP_Reactors.md#reformer1d_c3) |
| Downstream chemical reactors (electrolysis mode) | [`Reactors.MeOH`](library-reference/OpenTEMPEST_BOP_Reactors_MeOH.md) (methanol synthesis) |
| Power electronics | [`PowerElectronics`](library-reference/OpenTEMPEST_BOP.md#powerelectronics) |

Flow components (sources, splitters, manifolds, sensors) are in [`OpenTEMPEST.Flow`](library-reference/OpenTEMPEST_Flow.md) and control blocks in [`OpenTEMPEST.Blocks`](library-reference/OpenTEMPEST_Blocks.md). See [Build a system model](UserWorkflows.md#7-build-a-system-model) for a worked example.


[Return to main page](../README.md)
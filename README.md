<p align="center">
  <img src="./docs/TEMPESTlogo.svg" alt="TEMPEST logo" width="500">
</p>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21771899.svg)](https://doi.org/10.5281/zenodo.21771899)



Thank you for visiting. We're still working on the documentation!
For now feel free to try out the models in your Modelica environment (We have tested it with dymola so far)

**openTEMPEST** (**T**ransient **E**lectrochemical **M**odel for **P**rocess and **E**nergy **S**ys**T**ems) is an open-source Modelica/Dymola simulation framework for modeling and analyzing electrochemical systems at multiple scales from single cells to full integrated systems. The framework enables the detailed investigation of steady-state and transient performance and allows a reliable extrapolation of cell to stack and multi-stack behavior. As of now, openTEMPEST covers high temperature solid oxide cell fuel cell and electrolysis (SOFC/SOEC) as well as alkaline electrolysis (AEL) systems. 

<!---
**STILL MISSING:**
- Contact: good to have a generic DLR mail / additional person 
- Blog post, how to use OpenTEMPEST, some how-tos
- Screenshot-based tutorials how to start with given cell data, most effectively, if we link it to one of our own papers
- how to cite openTEMPEST, probably doing a "V1" updated on ZENODO with a DOI
--->
## Requirements
- Dymola (no guarantee that the models run with OpenModelica or other Modelica environments, but we'll be working on it)
- Modelica standard Library (3.2.3) - we will be updating to MSL 4 later in the year
- Thermopower library (https://github.com/casella/ThermoPower)

## Who Should Use This Framework?

- Researchers scaling up SOC/AEL systems from lab to application
- System and control engineers evaluating dynamic behavior of electrochemical systems, potentially coupled to upstream or downstream processes
- Experimental researchers validating electrochemical performance at system level
- Lecturers teaching electrochemical system design


## Documentation Index

Click a topic below to dive deeper:

- [Framework Structure & Modeling Levels](docs/FrameworkStructure.md)
- [Package Installation and First Steps](docs/Installation.md)
- [Common Workflows & Examples](docs/UserWorkflows.md): running examples, building cell, stack and system models, swapping the electrochemistry, extending the library
- [Library Reference](docs/library-reference/README.md): every package, model and function with its parameters and connectors
- [Writing Documentation](docs/DocumentationStyle.md): how to document models and regenerate the library reference


## List of Publications

| Year | Title                                                                                                                                                      | DOI                            | Description                                                                                                     | Citation                                                      |     |
| ------| ------------------------------------------------------------------------------------------------------------------------------------------------------------| --------------------------------| -----------------------------------------------------------------------------------------------------------------| ---------------------------------------------------------------| -----|
| 2026 | Improving Low-Load Efficiency of Marine SOFC Powertrains by Means of Efficient Module Hot Standby Operation                                                | 10.1149/1945-7111/ae7479       | Part-load operation of an SOFC system for maritime drivetrains                                                  | J. Hollmann et al 2026 J. Electrochem. Soc. 173 124503        |     |
| 2025 | Control concept for solid oxide electrolysis reactors to promote highly transient operation in modular plants                                              | 10.1016/j.ecmx.2025.101236     | Development of a control concept for multi-stack modules with fluctuating power supply                          | R. Lorenz et al 2025 Energy Convers. Manag. X 28 101236       |     |
| 2024 | Solid Oxide Cell Reactor Model for Transient and Stationary Electrochemical H2O and CO2 Conversion Process Studies                                         | 10.1149/1945-7111/ad5e01       | Describes the model equations and the FVM approach utilised as well as the co-electrolyiss modelling approach   | F. Sedeqi et al 2024 J. Electrochem. Soc. 171 074507          |     |
| 2023 | Operation strategies for a flexible megawatt scale electrolysis system for synthesis gas and hydrogen production with direct air capture of carbon dioxide | 10.1039/d2se01473d             | Development of operating strategies for a SOEC system with upstream DAC                                         | M. Tomberg et al 2023 Sustainable Energy Fuels 7 (2): 471–484 |     |
| 2022 | A New Approach to Modeling Solid Oxide Cell Reactors with Multiple Stacks for Process System Simulation                                                    | 10.1149/1945-7111/ac7009       | Description of the stack simplification approach with experimental validation on a multi-stack module           | M. Tomberg et al 2022 J. Electrochem. Soc. 169 054530         |     |
| 2022 | Operation of a Solid Oxide Fuel Cell Reactor with Multiple Stacks in a Pressured System with Fuel Gas Recirculation                                        | 10.1002/ente.202101075         | Model of a multi-stack SOFC system in a pressurised environment and the challenges associated with such systems | M. Tomberg et al 2022 Energy Technol. 10                      |     |
| 2018 | Transient reversible solid oxide cell reactor operation – Experimentally validated modeling and analysis                                                   | 10.1016/j.apenergy.2018.09.186 | Description of the 1D+1D modelling framework                                                                     | S. Srikanth et al 2018 Appl. Energy 232 473–488               |     |

## Known Limitations
We are actively developing OpenTEMPEST, the work has been going on for several years and has been fruitful for our work, nevertheless there are still some limitations - some intentional due to the scope we work (process system level), and some are a consequence of the limited capacity to work on both framework development and research - these are a few of the ones we currently see.

- Developed for Dymola, we can't guarantee how well it works in other modelica software such as OpenModelica
- No built-in mechanical stress models for cell degradation (yet)
- No detailed microkinetic models - we don't plan to go in this direction because of computational burden
- Currently still based on MSL 3.2.3 and Thermopower - we are still working on a seamless way to migrate to newer libraries, but the current framework works well for our purposes and we hope for yours too

## Feedback or Questions?

We welcome feedback and collaboration. Open an [Issue](https://codeberg.org/DLR-TT-ESI-EHT/OpenTEMPEST/issues) or reach out to us: tempest@dlr.de

---

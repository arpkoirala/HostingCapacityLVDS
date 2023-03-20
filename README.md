# Hosting capacity of LVDS using StochasticPowerModels


StochasticPowerModels.jl is an extension package of PowerModels.jl for 
Stochastic (Optimal) Power Flow. It is designed to enable inclusion of 
uncertainty in Steady-State Power Network Optimization. 

In this project file, the base code of StochasticPowerModels.jl is utilized to obtain Stochastic PV hosting capacity of LVDS. Activating the environment should be enough using:
 ```
    using Pkg
    Pkg.activate(".")
    Pkg.instantiate()
```
## Core Problem Specification

- Stochastic Optimal Power Flow (sOPF) HC
- Deterministic HC


## Varients of HC

- Equal PV
- Unequal PV

## Core Stochastic Specification
For now, we only support Polynomial Chaos Expansion based stochastic HC and deterministic HC in IV formulation.

- Polynomial Chaos Expansion
    - with auxiliary variables/constraints

## Network Data with Stochastic Data Extension

- Matpower ".m" files, extended to include:
    - stochastic germ: `mpc.sdata`,
    - stochastic bus data: `mpc.bus_sdata`, including: `dst_id`, `μ`, `σ`, `λvmin` and `λvmax`,
    - stochastic gen data: `mpc.gen_sdata`, including: `λpmin`, `λpmax`, `λqmin` and `λqmax`, and
    - stochastic branch data: `mpc.branch_sdata`, including: `λcmax`.

For an example of the base code, the user is referred to `/test/data/matpower/case5_spm.m`

## Example in Hosting Capacity Calculation
- For an example in different type of HC calculation:
    - `examples/case_spm_HC.jl` is pointed out where the codes for HC calculation to obtain Fig. 6, 7 and 8 can be found of the following paper:
    
    Koirala, Arpan; Van Aacker, Tom; Hashmi, Md Umar; D'hulst, Reinhilde; Van Hertem, Dirk (2022): Chance-constrained optimization based PV hosting capacity calculation using general Polynomial Chaos. TechRxiv. Preprint. https://doi.org/10.36227/techrxiv.21394950.v1 
    
## Network Data with Stochastic Data Extension
The Folder /text/data/Spanish/All_feeder has Spanish network on folder in JSON format. 

The file `test/data/Spanish/CreatePMDDictionary.jl` is the parser file to convert the JSON file into `PowerModels.jl` format.

The original dataset consists of a full Low voltage network of sub-urban region with 30 transformers, 160 feeders, 10290 nodes and 8087 consumers, with load profiles of 20 days from actual smart-meter.
The paper is available in https://www.sciencedirect.com/science/article/pii/S0142061519318836
and the link for full data-set: https://data.mendeley.com/datasets/685vgp64sm/1

For the original networks, the line impedance is specified 4x4 matrice without mutual impedance and the load from smart meter data for 20 days.

However, in this part the load and irradiance are defined as Beta distribution in folders `beta_lm_2016_8_6.csv` and `beta_pm_2016_8_6.csv` respectively for a high irradiance day in spring. 

For each feeder there are 4 JSON file describing the feeder topology:
	
- *_configuration.json, 
- *_branches.json, 
- *_buses.json, and 
- *_devices.json 
and 1 csv file:
- *.csv- which is the linking file between devices and the load uncertainty. 

## Example in Hosting Capacity Calculation for Non-synthetic European LVDS
- The HC calculation of all feeders:
    - `examples/case_spm_HC_all_feeder.jl` is pointed out where the codes for HC calculation for all feeder and code to obtain Fig. 11 of the following paper:
    
    Koirala, Arpan; Van Aacker, Tom; Hashmi, Md Umar; D'hulst, Reinhilde; Van Hertem, Dirk (2022): Chance-constrained optimization based PV hosting capacity calculation using general Polynomial Chaos. TechRxiv. Preprint. https://doi.org/10.1109/TPWRS.2023.3258550 

    
## Installation
For this application the installation of StochasticPowerModels.jl is not suggested as this application might not be updated along with the changes in the main package and will remain as an standalone application.

## Acknowledgements

The primary developer is StochasticPowerModels.jl Tom Van Acker ([@timmyfaraday](https://github.com/timmyfaraday)), 
with support from the following contributors:
- Arpan Koirala ([@arpkoirala](https://github.com/arpkoirala)), KU Leuven, ACR formulation
- Frederik Geth ([@frederikgeth](https://github.com/frederikgeth)), CSIRO, reduced IVR formulation

The latest stable release of StochasticPowerModels can be obtained at:

```
https://github.com/timmyfaraday/StochasticPowerModels.jl.git
```
## License

This code is provided under a BSD license.

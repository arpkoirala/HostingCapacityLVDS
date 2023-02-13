# Hosting capacity of LVDS using StochasticPowerModels


StochasticPowerModels.jl is an extension package of PowerModels.jl for 
Stochastic (Optimal) Power Flow. It is designed to enable inclusion of 
uncertainty in Steady-State Power Network Optimization. In this branch of StochasticPowerModels.jl, the base code of StochasticPowerModels.jl is utilized to obtain Stochastic PV hosting capacity of LVDS.

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
    
    

## Installation

The latest stable release of StochasticPowerModels can be installed using the 
Julia package manager:

```
] add https://github.com/timmyfaraday/StochasticPowerModels.jl.git
```

In order to test whether the package works, run:

```
] test StochasticPowerModels
```

## Acknowledgements

The primary developer is Tom Van Acker ([@timmyfaraday](https://github.com/timmyfaraday)), 
with support from the following contributors:
- Arpan Koirala ([@arpkoirala](https://github.com/arpkoirala)), KU Leuven, ACR formulation
- Frederik Geth ([@frederikgeth](https://github.com/frederikgeth)), CSIRO, reduced IVR formulation

## License

This code is provided under a BSD license.

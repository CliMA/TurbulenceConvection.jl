# TurbulenceConvection.jl #

|||
|---------------------:|:----------------------------------------------|
| **Docs Build**       | [![docs build][docs-bld-img]][docs-bld-url]   |
| **Documentation**    | [![dev][docs-dev-img]][docs-dev-url]          |
| **GHA CI**           | [![gha ci][gha-ci-img]][gha-ci-url]           |
| **Code Coverage**    | [![codecov][codecov-img]][codecov-url]        |
| **Bors enabled**     | [![bors][bors-img]][bors-url]                 |

[docs-bld-img]: https://github.com/CliMA/TurbulenceConvection.jl/actions/workflows/docs.yml/badge.svg
[docs-bld-url]: https://github.com/CliMA/TurbulenceConvection.jl/actions/workflows/docs.yml

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://CliMA.github.io/TurbulenceConvection.jl/dev/

[gha-ci-img]: https://github.com/CliMA/TurbulenceConvection.jl/actions/workflows/ci.yml/badge.svg
[gha-ci-url]: https://github.com/CliMA/TurbulenceConvection.jl/actions/workflows/ci.yml

[codecov-img]: https://codecov.io/gh/CliMA/TurbulenceConvection.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/CliMA/TurbulenceConvection.jl

[bors-img]: https://bors.tech/images/badge_small.svg
[bors-url]: https://app.bors.tech/repositories/35146


TurbulenceConvection (Single Column Atmospheric Model in Julia) provides a framework for testing parameterizations of clouds and turbulence.
It is particularly designed to support eddy-diffusivity mass-flux modeling frameworks.

Information about the EDMF parameterization implemented in TurbulenceConvection can be found in:

Tan, Z., C. M. Kaul, K. G. Pressel, Y. Cohen, T. Schneider, and J. Teixeira, 2018:
An extended eddy-diffusivity mass-flux scheme for unified representation of
subgrid-scale turbulence and convection. Journal of Advances in Modeling Earth Systems, 2018.
(see [Tan et al., 2018](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1002/2017MS001162)).

The code is written in Julia, and was translated from [SCAMPy](https://github.com/CliMA/SCAMPy) for the commit 496dad0c2438235684823511cacbf5761d6a237c.

Code Contributors:
	Colleen Kaul (PNNL) --initial/primary developer,
	Yair Cohen (Caltech),
        Jia He (Caltech),
	Anna Jaruga (JPL/Caltech),
        Ignacio Lopez-Gomez (Caltech),
	Kyle Pressel (PNNL),
	Charles Kawczynski (Caltech),

Additional Acknowledgements:
	Tapio Schneider (Caltech),
	Joao Teixeira (JPL).

# User Installation #

Installation is easy:

```julia-repl
(@v1.x) pkg> add https://github.com/CliMA/TurbulenceConvection.jl
```

If you plan to develop TurbulenceConvection.jl, you may want to clone it instead:


```
git clone https://github.com/CliMA/TurbulenceConvection.jl
cd TurbulenceConvection.jl
julia --project
julia> ]
pkg> instantiate
```

# Running #

```
$ cd TurbulenceConvection.jl
```

TurbulenceConvection.jl can be run in the same way that SCAMPy used to, given one of the following cases:

 - Bomex
 - life_cycle_Tan2018
 - Soares
 - Rico
 - TRMM_LBA
 - ARM_SGP
 - GATE_III
 - DYCOMS_RF01
 - GABLS
 - SP

with, for example:

```
julia --project integration_tests/utils/generate_namelist.jl Soares
julia --project integration_tests/utils/main.jl Soares

```
or, interactively, with
```julia-repl
julia --project
julia> include(joinpath("integration_tests", "Soares.jl"))
```

# Automated plotting  #

Upon running a particular experiment (described above), comparison plots (against [SCAMPy](https://github.com/CliMA/SCAMPy)) are automatically generated in, for example, `Output.Bomex.01/stats/comparison/`.

# Table of prognostic and diagnostic variables

 - Prognostic [✓]
 - Diagnostic [x]
 - NA [-]

| **Variable** | **Grid-mean** | **Environment** | **Updrafts** |
| ------------ | ------------- | --------------- | ------------ |
| `u`          | [✓]           |  [x]            |  [✓]         |
| `v`          | [✓]           |  [x]            |  [✓]         |
| `w`          | [-]           |  [x]            |  [✓]         |
| `θ_liq_ice`  | [✓]           |  [x]            |  [✓]         |
| `q_tot`      | [✓]           |  [x]            |  [✓]         |
| `a`          | [-]           |  [x]            |  [✓]         |
| `tke`        | [x]           |  [✓]            |  [-]         |



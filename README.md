# TurbulenceConvection.jl #

|||
|---------------------:|:----------------------------------------------|
| **Documentation**    | [![dev][docs-dev-img]][docs-dev-url]          |
| **DOI**              | [![DOI][zenodo-img]][zenodo-latest-url]       |
| **Docs Build**       | [![docs build][docs-bld-img]][docs-bld-url]   |
| **GHA CI**           | [![gha ci][gha-ci-img]][gha-ci-url]           |
| **Code Coverage**    | [![codecov][codecov-img]][codecov-url]        |
| **Bors enabled**     | [![bors][bors-img]][bors-url]                 |

[zenodo-img]: https://zenodo.org/badge/DOI/10.5281/zenodo.6392396.svg
[zenodo-latest-url]: https://doi.org/10.5281/zenodo.6392396

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

Tan, Z., Kaul, C. M., Pressel, K. G., Cohen, Y., Schneider, T., and Teixeira, J. (2018)
**An extended eddy-diffusivity mass-flux scheme for unified representation of
subgrid-scale turbulence and convection.** *Journal of Advances in Modeling Earth Systems*. [doi](https://doi.org/10.1002/2017MS001162)

Cohen, Y., Lopez-Gomez, I., Jaruga, A., He, J., Kaul, C., and Schneider, T. (2020) **Unified entrainment and detrainment closures for extended eddy-diffusivity mass-flux schemes.** *Journal of Advances in Modeling Earth Systems*, 12, e2020MS002162. [doi](https://doi.org/10.1029/2020MS002162)

Lopez-Gomez, I., Cohen, Y., He, J., Jaruga, A., Schneider, T. (2020) **A generalized mixing length closure for eddy-diﬀusivity mass-flux schemes of turbulence and convection.** *Journal of Advances in Modeling Earth Systems*, 12, e2020MS002161. [doi](https://doi.org/10.1029/2020MS002161)

He, J., Cohen, Y., Lopez-Gomez, I., Jaruga, A., Schneider, T. (2021) **An improved perturbation pressure closure for eddy-diffusivity mass-flux schemes**. [preprint](https://doi.org/10.1002/essoar.10505084.2)

The code is written in Julia, and was translated from [SCAMPy](https://github.com/CliMA/SCAMPy) for the commit 496dad0c2438235684823511cacbf5761d6a237c.

Code Contributors (alphabetical):
    Yair Cohen (Caltech),
    Jia He (Caltech),
    Anna Jaruga (Caltech),
    Colleen Kaul (PNNL) --initial/primary developer of SCAMPy,
    Charles Kawczynski (Caltech),
    Ignacio Lopez-Gomez (Caltech),
    Kyle Pressel (PNNL),

Additional Acknowledgements:
    Tapio Schneider (Caltech),
    Joao Teixeira (JPL).

# User Installation #

Installation is easy:

```julia-repl
(@v1.6) pkg> add TurbulenceConvection
```

TurbulenceConvection.jl requires Julia 1.6 or higher. If you plan to develop TurbulenceConvection.jl, you may want to clone it instead:


```
git clone https://github.com/CliMA/TurbulenceConvection.jl
cd TurbulenceConvection.jl
julia --project=integration_tests
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
 - DYCOMS_RF02
 - GABLS

with, for example:

```
julia --project=integration_tests driver/generate_namelist.jl Soares
julia --project=integration_tests driver/main.jl Soares
```
or by calling our integraion test driver,
```
julia --project=integration_tests integration_tests/driver.jl --case Soares
```
or, interactively, with
```julia-repl
julia --project=integration_tests
julia> case_name = "Soares" # default is "Bomex"
julia> include(joinpath("integration_tests", "driver.jl"))
```

# Automated plotting  #

Upon running a particular experiment (described above), comparison plots (against [SCAMPy](https://github.com/CliMA/SCAMPy)) are automatically generated in, for example, `Output.Bomex.01/stats/comparison/`.

# Table of prognostic and diagnostic variables

 - Prognostic [✓]
 - Diagnostic [x]
 - NA [-]

| **Variable**  | **Grid-mean** | **Environment** | **Updrafts** |
| ------------- | ------------- | --------------- | ------------ |
| `u`           | [✓]           |  [x]            |  [x]         |
| `v`           | [✓]           |  [x]            |  [x]         |
| `w`           | [-]           |  [x]            |  [✓]         |
| `θ_liq_ice`   | [✓]           |  [x]            |  [✓]         |
| `q_tot`       | [✓]           |  [x]            |  [✓]         |
| `a`           | [-]           |  [x]            |  [✓]         |
| `tke`         | [x]           |  [✓]            |  [-]         |
| `θ′θ′`        | [x]           |  [✓]            |  [-]         |
| `q_tot′q_tot′`| [x]           |  [✓]            |  [-]         |
| `θ′q_tot′`    | [x]           |  [✓]            |  [-]         |


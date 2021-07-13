# TurbulenceConvection.jl #

TurbulenceConvection (Single Column Atmospheric Model in Python) provides a framework for testing parameterizations of clouds and turbulence.
It is particularly designed to support eddy-diffusivity mass-flux modeling frameworks.

Information about the EDMF parameterization implemented in TurbulenceConvection can be found in:

Tan, Z., C. M. Kaul, K. G. Pressel, Y. Cohen, T. Schneider, and J. Teixeira, 2018:
An extended eddy-diffusivity mass-flux scheme for unified representation of
subgrid-scale turbulence and convection. Journal of Advances in Modeling Earth Systems, 2018.
(see [Tan et al., 2018](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1002/2017MS001162)).

The code is written in Julia, and was translated from [SCAMPy](https://github.com/CliMA/SCAMPy) on 7/13/2021.

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

# Installation #

Installation is easy:

```julia-repl
(@v1.x) pkg> add https://github.com/CliMA/TurbulenceConvection.jl
```

# Running #

```
$ cd TurbulenceConvection.jl
```

Note that, we will soon support running TurbulenceConvection.jl in the same way that SCAMPy used to. For now, run one of the cases:

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

with

```
julia --project integration_tests/Soares.jl
```
or, interactively, with
```julia-repl
julia --project
julia> include(joinpath("integration_tests", "Soares.jl"))
```

# automated plotting  #

All automated plots are located in the `viz/` folder and are generated. This section is under construction.

To run all automatic plots please try:

```
$ cd tests/

$ py.test -s -v plots/

$ cd ../

```
To run an individual plot please try:

```
$ cd tests/

$ py.test -s -v plots/test_plot_Soares.py

$ cd ../

```

The subfolder TurbulenceConvection/tests/les_data contains several netCDF files. These are reduced data files from the stats files of pycles
in which only the relevant data is saved. The code "reduce_pycles_netcdf.py" reduces pycles files given the input location of a pycles file.



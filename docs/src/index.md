# TurbulenceConvection.jl

```@meta
CurrentModule = TurbulenceConvection
```

TurbulenceConvection.jl (TC.jl) is a Julia implementation of the Extended Eddy-Diffusivity Mass-Flux (EDMF) model published in [^1], [^2], and [^3]. The package contains a Single Column Model (SCM) which solvers for a single column of a climate model with vertical fluxes computed by the dynamics of the EDMF's 'updrafts' and 'environment' subdomains.
The performance of TC.jl is routinely monitored by a Continuous Integration (CI) and plots of model variances from a range of simulations are plotted in Buildkite. These plots can be viewed via the colored check (`✓`, `x` or `⦿`) near the commit tag in the 'code' page on github web interface.

The Extended EDMF model can be run with several specifications such as:
 - Case: what model setup is run. typically corresponding to an observation campaign (i.e. DYCOMS_RF01 for a standard stratocumulus case) or an LES driven simulation of the SCM based on CliMA's LES library [^4].
- Closure types: such as entrainment, mixing length etc.
- Model specifications: such as forcing properties, domain and grid etc.
- Parameters: Such as the physical parameters listed in Table 2 in [^2] and Table 1 in [^3] as well as parameters relating to machine learning models (Fourier Neural Operators, Neural Networks etc.)

The package depends on the following clima packages: [ClimaCore.jl](https://github.com/CliMA/ClimaCore.jl), [Thermodyanmics.jl](https://github.com/CliMA/Thermodynamics.jl), [SurfaceFlux.jl](https://github.com/CliMA/SurfaceFluxes.jl), [OperatorFlux.jl](https://github.com/CliMA/OperatorFlux.jl)



## Authors

CalibrateEDMF.jl is being developed by the [Climate Modeling Alliance](https://clima.caltech.edu). The main developers are Charles Kawczynski, Yair Cohen, Anna Jaruga, Ignacio Lopez-Gomez, Haakon Ludvig Langeland Ervik, Costa Christopoulos.

## References

[^1]: [Tan2018](@cite)
[^2]: [Cohen2020](@cite)
[^3]: [Lopez2020](@cite)
[^4]: [Shen2022](@cite)
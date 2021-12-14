# Adding closures to EDMF

Entrainment closure functions are written in a file that can be augmented with additional closure types. The type of entrainment closure is specified in the namelist, which is a collection of model parameters, flags, and other options used to run EDMF. After the new closure is implemented, EDMF can be run with your new closure by running `julia --project integration_tests/LES_driven_SCM.jl` in the top level of the repo (locally), or with the submit script below on the Caltech HPC cluster.

## Adding a general closure function:
*1.)* Add a `entr_detr` function defintion to `src/closures/entr_detr.jl` which implements your closure. The code uses multiple 
dispatch to determine which type of entrainment to used based on the `εδ_model_type` argument. All `entr_detr` 
definitions must take:
		`param_set` (a collection of parameters, many of which are defined in the namelist, as well as parameters needed for your closure)
		`εδ_model_vars` - a structure containing all input variables your closure might depend on. 
		`εδ_model_type` - a type specifying the closure, defined in `src/types.jl`


Currently, the closures we're after are functions of local, non-dimensional groups of the prognostic model variables. These are returned by the `non_dimensional_groups` function. **New ML models should take these non-dimensional groups as inputs and return two continuous, non-negative outputs (one for entrainment `nondim_ε` and 
one for detrainment `nondim_δ`)**. Note we are not modifying the dimensional part of entrainment/detrainment (`dim_scale`), but only the non-dimensional component multiplying it.


```julia
function entr_detr(param_set, εδ_model_vars, εδ_model_type::<YourClosureType>)
	c_gen = ICP.c_gen(param_set)
	dim_scale = dimensional_part(param_set, εδ_model_vars)
	area_limiter = max_area_limiter(param_set, εδ_model_vars)

	foo =  <Your Non-Dimensional Function>

	nondim_groups = non_dimensional_groups(param_set, εδ_model_vars)
	nondim_ε, nondim_δ = <foo(nondim_groups)>

	# dynamic entrainment / detrainment
	ε_dyn = dim_scale * nondim_ε
	δ_dyn = dim_scale * nondim_δ + area_limiter

	# turbulent entrainment
	ε_turb, K_ε = compute_turbulent_entrainment(param_set, εδ_model_vars)

	return EntrDetr(ε_dyn, δ_dyn, ε_turb, K_ε)
end
```
Most of this code is boiler plate - you only need to make changes at the bracketed components. `c_gen` is a vector of parameters needed for your closure, specified in the namelist. 


*2.)* Update `driver/generate_namelist.jl` to specify your closure type and initial parameters.
 Add initial parameters as the `c_gen` namelist argument as a vector.

```julia
namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment"] = <Name of your closure>
```
.
.
.
```julia
namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["general_ent_params"] = <Your initial parameters.>
```

*3.)* Add a type definition for your closure in `src/types.jl`
Create a struct with a unique name for your closure.

```julia
struct <YourClosureType> end
```
.
.
.
Add logic that maps namelist `entrainment` name to the type you defined above. 


```julia
	valid_options = ["moisture_deficit", "NN", ... <Name of your closure>],
```
.
.
.
```julia
elseif entr_type == <Name of your closure>
	<YourClosureType>
```


## An example: Simple Neural Network Closure

Here we've added a simple neural network which maps non-dimensional groups (`nondim_groups`) to the non-dimensional components
of entrainment & detrainment (`nondim_ε` & `nondim_δ`, respectively). 
The entrainment type is `NNEntr` and the namelist name is "NN".


#### src/closures/entr_detr.jl
```julia
function entr_detr(param_set, εδ_model_vars, εδ_model_type::NNEntr)
    c_gen = ICP.c_gen(param_set)
    dim_scale = dimensional_part(param_set, εδ_model_vars)
    area_limiter = max_area_limiter(param_set, εδ_model_vars)

    # Neural network closure
    nn_arc = (4, 2, 2)  # (#inputs, #neurons, #outputs)
    nn_model = Chain(
        Dense(reshape(c_gen[1:8], nn_arc[2], nn_arc[1]), c_gen[9:10], sigmoid),
        Dense(reshape(c_gen[11:14], nn_arc[3], nn_arc[2]), c_gen[15:16], softplus),
    )

    nondim_groups = non_dimensional_groups(param_set, εδ_model_vars)
    nondim_ε, nondim_δ = nn_model(nondim_groups)

    # dynamic entrainment / detrainment
    ε_dyn = dim_scale * nondim_ε
    δ_dyn = dim_scale * nondim_δ + area_limiter

    # turbulent entrainment
    ε_turb, K_ε = compute_turbulent_entrainment(param_set, εδ_model_vars)

    return EntrDetr(ε_dyn, δ_dyn, ε_turb, K_ε)
end
```


#### driver/generate_namelist.jl
```julia
namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment"] = "NN" 
```
. 
.
.
```julia
namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["general_ent_params"] =
	SA.SVector(0.3038, 0.719,-0.910,-0.483,
				0.739, 0.0755, 0.178, 0.521,
				0.0, 0.0, 0.843,-0.340,
				0.655, 0.113, 0.0, 0.0)
```

#### `src/types.jl`
```julia
struct NNEntr end
```
.
.
.
```julia
	valid_options = ["moisture_deficit", "NN"],
```
. 
.
.
```julia
elseif entr_type == "NN"
    NNEntr()
```
## Submitting TC.jl on the Caltech cluster
Use the following sbatch script (or something similar) to submit TC.jl on the Caltech central cluster.

```bash
#!/bin/bash

#Submit this script with: sbatch submit_TC

#SBATCH --time=01:00:00   # walltime
#SBATCH --ntasks=10  # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --reservation=clima
#SBATCH -J "TC_jl"   # job name

module purge
module load julia/1.7.0 hdf5/1.10.1 netcdf-c/4.6.1 openmpi/4.0.1
julia --project integration_tests/LES_driven_SCM.jl
```
# The Extended Eddy-Diffusivity Mass-Flux scheme

The Eddy-Diffusivity Mass-Flux (hereafter EDMF) scheme is a unified parameterization scheme for subgrid scale (SGS) turbulence and convection, including clouds which are imperative for the global radiative balance. The EDMF decomposes horizontally the grid box of a low resolution Earth System Model (ESM) into two or more "subdomains": namely a turbulent environment and one or more coherent updrafts. The EDMF scheme solves equations for first and second moment in these subdomains, from which the first, second and third (skewness) moments in each grid box of a ESM can be computed. The grid box skewness is essential for the representation of extreme events.
Our Extended EDMF model is a time dependent version of the traditional EDMF. The model requires several closures for unknown processes, which are the crux of the problem in simulating turbulence and clouds correctly in a changing climate. 
These closures include:
 * Entrainment and detrainment - mass exchange between the subdomains
 * Non-hydrostatic (i.e. unbalanced by gravity) pressure effects
 * Eddy diffusivity and mixing length to account for turbulent mixing
 * Microphysical cloud process such as the formation and removal of precipitation
 * Radiative effects

These closures that are the foci of our machine learning effort, and specifically the entrainment and detrainment that are notorious effecting the global behavior of climate models.

A conceptual sketch of the EDMF scheme in an Earth Systems Model grid box is shown in the figure below. The projected distribution of the combined environment and updrafts as well as the various physical processes requiring closures are sketched.
![alt text](figures/combined_globe_edmf.pdf)

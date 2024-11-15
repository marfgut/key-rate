# Key rate calculation in QKD with Julia
This repository contains Julia code for calculating the key rate in quantum key distribution (QKD) under realistic noise models, as part of my Bachelor's thesis in Physics. The code is based on a converging semidefinite programming (SDP) hierarchy of lower bounds on the conditional von Neumann entropy with mutually unbiased bases (MUB). It uses JuMP for modeling and the MOSEK solver to perform the necessary computations.

## Files
### `functions.jl`: contains the main code, where all the necessary functions for performing the calculations are defined. 
### `plots.jl`: contains scripts for generating **plots** to visualize the results of the simulations and calculations performed in `functions.jl`. 
### `mub.jld2`: contains the mutually unbiased bases for dimensions 2 through 13. For dimensions 6, 10, and 12, the bases are only approximately unbiased.

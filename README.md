# SocialMediaModel

This repository contains the Julia code accompanying the paper 
"Modelling opinion dynamics under the impact of influencer and media strategies". 


## Components
In particular, the repository contains: 
- the agent-based formulation (`ABMsolve()`) and the PDE approximation (`PDEsolve()`) of an opinion model that models how individuals, influencers and media interact on social media platforms
- plotting functions (in `plotting.jl`)
- code to perform ensemble simulations (`ABMsolveensemble()` and `PDEsolveensemble`)
- methods to study optimal strategies of influencers and media with `PDEsolveopt()`. 


## Usage
To use the package in Julia, it has to be installed with
`Pkg.add https://github.com/LuzieH/SocialMediaModel.git`
and then import the package with `using SocialMediaModel`. 
Most functions are exported and can be run e.g. with either `SocialMediaModel.ABMsolve()` or simply, `ABMsolve()`. 
Tests can be run with `runtests()`. 
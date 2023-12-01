# SocialMediaModel

This repository contains the Julia code accompanying the paper 
["Modelling opinion dynamics under the impact of influencer and media strategies"](https://doi.org/10.1038/s41598-023-46187-9). 

![](https://github.com/LuzieH/SocialMediaModel/blob/master/img/abm_single_4inf.gif)


## Components
In particular, the repository contains: 
- the agent-based formulation (`ABMsolve()`) and the PDE approximation (`PDEsolve()`) of an opinion model that models how individuals, influencers and media interact on social media platforms
- plotting functions (in `plotting.jl`)
- code to perform ensemble simulations (`ABMsolveensemble()` and `PDEsolveensemble`)
- methods to study optimal strategies of influencers and media with `PDEsolveopt()`. 


## Usage
To use the package in Julia, it has to be installed with
`Pkg.add https://github.com/LuzieH/SocialMediaModel`
and then imported with `using SocialMediaModel`. 
Most functions are exported, so instead of calling `SocialMediaModel.ABMsolve()`, they can simply be called with `ABMsolve()`. 
Tests can be run with `runtests()`. 


## Citation
If you use this code in your work, we politely ask you to acknowledge it in your manuscript by citing:
1. Helfmann, L., Djurdjevac Conrad, N., Lorenz-Spreen, P., & Schütte, C. (2023). Supplementary code for the paper Modelling opinion dynamics under the impact of influencer and media strategies. DOI: 10.12752/9267
2. Helfmann, L., Djurdjevac Conrad, N., Lorenz-Spreen, P., & Schütte, C. (2023). Modelling opinion dynamics under the impact of influencer and media strategies. Scientific Reports, 13(1), 19375. DOI: 10.1038/s41598-023-46187-9

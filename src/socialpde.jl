module socialpde
    include("setting.jl")
    include("pde.jl")
    include("abm.jl")
    include("ensemble.jl")
    include("plotting.jl")
    include("controlstrategies.jl")
    include("strategyexperiments.jl")
    include("tests.jl")

    export  ABMsolve,  ABMsolveplot, solve, solveplot, plotsnapshots, ABMplotsnapshots, gifsingle,
     plotfollowernumbers, ABMgifsingle, plotensemblesnapshots, solveensemble, ABMsolveensemble, runensembles,
     parameters, parametersstronginf, ABMconstruct, PDEconstruct, PDEconstructcoarse, PDEconstructmeso, parameterscontrol,
     solvefixedtargets, solveopt, GLmed, GLinf, infrightcorner, infcounteraction, medcounteraction,
     testabm, testpde, testensemble, testcontrol, testoptimization, runtests
end 
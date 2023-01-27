using Revise
using LinearAlgebra

LinearAlgebra.BLAS.set_num_threads(parse(Int,get(ENV, "SLURM_JOB_CPUS_PER_NODE", "20")))


print("Running on ", string(Threads.nthreads())," threads with ",string(LinearAlgebra.BLAS.get_num_threads())," BLAS")

includet("setting.jl")
includet("pde.jl")
includet("abm.jl")
includet("ensemble.jl")
includet("plotting.jl")
includet("controlstrategies.jl")
includet("strategyexperiments.jl")
includet("tests.jl")
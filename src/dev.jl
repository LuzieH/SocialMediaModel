""" This is the entry file, to be included, for
using this not as a package but a simple project """

using Revise
using LinearAlgebra

LinearAlgebra.BLAS.set_num_threads(parse(Int,get(ENV, "SLURM_JOB_CPUS_PER_NODE", "20")))
print("Running on ", string(Threads.nthreads())," threads with ",string(LinearAlgebra.BLAS.get_num_threads())," BLAS threads.")

includet("setting.jl")
includet("pde.jl")
includet("abm.jl")
includet("ensemble.jl")
includet("plotting.jl")
includet("controlstrategies.jl")
includet("strategyexperiments.jl")
includet("tests.jl")
includet("paperfigures.jl")


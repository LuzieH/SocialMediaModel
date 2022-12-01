using Revise

ENV["MATLABPATH"]="/home/htc/bzfhelfm/.julia/packages/Mex/5WarT/mexjulia"
ENV["MATLAB_ROOT"] = "/software/Matlab/current"

includet("pde.jl")
includet("abm.jl")
includet("control.jl")
includet("run.jl")
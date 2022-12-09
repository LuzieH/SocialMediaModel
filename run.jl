quickrun() = solveopt(; p = PDEconstructcoarse(), q= parameters_control(), r=parameters_optcont(ntarg=8), alg=:LN_COBYLA, mtime = -1, meval=100)

function normalrun(;mtime = -1, meval=100) 
    return solveopt(; p = PDEconstructcoarse(), q= parameters_control(), r=parameters_optcont(ntarg=8), alg=:LN_COBYLA, mtime =mtime,meval=meval)
end

function maxcorner(;mtime = -1, meval=100, alg=:LN_COBYLA) 
    return solveopt(; p = PDEconstructmeso(), q= parameters_control(), r=parameters_optcont(ntarg=8,maximize="corner",start="zero"), alg=alg, mtime =mtime,meval=meval)
end

function maxinf(;mtime = -1, meval=100) 
    return solveopt(; p = PDEconstructcoarse(), q= parameters_control(), r=parameters_optcont(ntarg=8,maximize="follower",start="mean"), alg=:LN_COBYLA, mtime =mtime,meval=meval)
end
 

function globalcorner(;globaleval = 100, localeval = 100,localtime = -1)
    return solveopt_global(;globaleval = globaleval, localeval = localeval, localtime = localtime, p = PDEconstructcoarse(), q= parameters_control(), r=parameters_optcont(ntarg=8,maximize="corner",start="zero"), alg=:LN_COBYLA, ftol_rel=1e-4, xtol_rel = 1e-2)
end
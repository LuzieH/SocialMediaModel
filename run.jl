quickrun() = solveopt(; p = PDEconstructcoarse(), q= parameters_control(), r=parameters_optcont(ntarg=3, Tmax=1.), alg=:LN_COBYLA, mtime = -1, meval=100)

fullrun() = solveopt(; p = PDEconstructcoarse(), q= parameters_control(), r=parameters_optcont(ntarg=8, Tmax=5.), alg=:LN_COBYLA, mtime = -1, meval=1000)

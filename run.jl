quickrun() = solveopt(; p = PDEconstructcoarse(), q= parameters_control(), r=parameters_optcont(ntarg=8), alg=:LN_COBYLA, mtime = -1, meval=100)

largebounds() = solveopt(; p = PDEconstructcoarse(), q= parameters_control(), r=parameters_optcont(ntarg=8,speedbound=10), alg=:LN_COBYLA, mtime = -1, meval=100)

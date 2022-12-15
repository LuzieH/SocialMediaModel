quickrun() = solveopt(; p = PDEconstructcoarse(), q= parameters_control(), r=parameters_optcont(ntarg=8), alg=:LN_COBYLA, mtime = -1, meval=100)

function normalrun(;mtime = -1, meval=100) 
    return solveopt(; p = PDEconstructcoarse(), q= parameters_control(), r=parameters_optcont(ntarg=8), alg=:LN_COBYLA, mtime =mtime,meval=meval)
end

function maxcorner(;mtime = -1, meval=100, alg=:LN_COBYLA) 
    return solveopt(; p = PDEconstructmeso(), q= parameters_control(), r=parameters_optcont(ntarg=8,maximize="corner",start="zero"), alg=alg, mtime =mtime,meval=meval, ftol_rel=1e-6, xtol_rel = 1e-3,x0 = rand(Uniform(-2,2),2*r.ntarg))
end

function maxinf(;mtime = -1, meval=100) 
    return solveopt(; p = PDEconstructcoarse(), q= parameters_control(), r=parameters_optcont(ntarg=8,maximize="follower",start="mean"), alg=:LN_COBYLA, mtime =mtime,meval=meval)
end
 

function globalcorner(;globaleval = 100, localeval = 100, localtime = -1)
    return solveopt_global(;globaleval = globaleval, localeval = localeval, localtime = localtime, p = PDEconstructcoarse(), q= parameters_control(), r=parameters_optcont(ntarg=8,maximize="corner",start="zero"), alg=:LN_COBYLA, ftol_rel=1e-4, xtol_rel = 1e-2)
end

function counter_med(;mtime = -1, meval=100, alg=:LN_COBYLA)
    return solveopt(; p = PDEconstructcoarse(), r=parameters_optcont(ntarg=8,maximize="counter_follower",start="zero"), alg=alg, mtime = mtime, meval=meval,countercontrol="med",stubborntarget=[1.5 1.5],x0 = rand(Uniform(-2,2),2*r.ntarg))
end


function globalmultistart(;mtime = 800000, meval=-1, alg=:LN_COBYLA)
    return solveopt(; p = PDEconstructcoarse(), q= parameters_control(),r=parameters_optcont(ntarg=8,maximize="corner",start="zero"), alg=alg, mtime = mtime, meval=meval, ftol_rel=1e-2, xtol_rel = 1e-1, multistart =true)
end


# function maxinf1step(;mtime = -1, meval=100) 
#     return solveopt(; p = PDEconstructcoarse(), q= parameters_control(), r=parameters_optcont(ntarg=1,maximize="follower",start="zero"), alg=:LN_COBYLA, mtime =mtime,meval=meval, ftol_rel=1e-3, xtol_rel = 1e-2, multistart =true)
# end

# function countermed1step(;mtime = -1, meval=100)
#     return solveopt(; p = PDEconstructcoarse(), r=parameters_optcont(ntarg=1,maximize="counter_follower",start="zero",alpha=0.005), alg=:LN_COBYLA, mtime = mtime, meval=meval, ftol_rel=1e-3, xtol_rel = 1e-2, multistart =true,countercontrol="med",stubborntarget=[1.5 1.5])
# end #if speed penalty is too large, counter control will not move at all, since the effect is anyway only small

# function counterinf1step(;mtime = -1, meval=100)
#     return solveopt(; p = PDEconstructcoarse(), r=parameters_optcont(ntarg=1,maximize="counter_follower",start="zero",alpha=0.005), alg=:LN_COBYLA, mtime = mtime, meval=meval, ftol_rel=1e-3, xtol_rel = 1e-2, multistart =true,countercontrol="inf",stubborntarget=[1.5 1.5])
# end



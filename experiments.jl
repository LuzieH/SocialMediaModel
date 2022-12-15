
function counterglobal(;mtime = 800000, meval=-1, alg=:LN_COBYLA, ntarg = 3,countercontrol = "med", p = PDEconstructmeso())
    return solveopt(; p = p, q= parameters_control(),r=parameters_optcont(ntarg=ntarg,maximize="counter_follower",start="zero"), alg=alg, mtime = mtime, meval=meval, ftol_rel=1e-4, xtol_rel = 1e-2, multistart =true,countercontrol=countercontrol,stubborntarget=[1.5 1.5],x0 = rand(Uniform(-2,2),2*ntarg))
end

function counterlocal(x0; mtime = 800000, meval=-1, alg=:LN_COBYLA, ntarg = 3, countercontrol = "med", p = PDEconstructmeso())
    return solveopt(;p=p, q= parameters_control(),r=parameters_optcont(ntarg=ntarg,maximize="counter_follower",start="zero"), alg=alg, mtime = mtime, meval=meval, ftol_rel=1e-10, xtol_rel = 1e-5, multistart =false,countercontrol=countercontrol,stubborntarget=[1.5 1.5],x0 = x0)
end

#global first with meso and 5 days = 432000, then fine and ... days
#started thursday 15.12 around 12 oclock

function Gbob3med(;mtime =432000,meval = -1, p = PDEconstructmeso())
    return counterglobal(;mtime = mtime, meval=meval, alg=:LN_BOBYQA, ntarg = 3,countercontrol = "med",p=p)
end

function Gbob3inf(;mtime = 432000,meval = -1, p = PDEconstructmeso())
    return counterglobal(;mtime = mtime, meval=meval, alg=:LN_BOBYQA, ntarg = 3,countercontrol = "inf",p=p)
end

function Gbob4med(;mtime = 432000,meval = -1, p = PDEconstructmeso())
    return counterglobal(;mtime = mtime, meval=meval, alg=:LN_BOBYQA, ntarg = 4,countercontrol = "med",p=p)
end

function Gbob4inf(;mtime = 432000,meval = -1, p = PDEconstructmeso())
    return counterglobal(;mtime = mtime, meval=meval, alg=:LN_BOBYQA, ntarg = 4,countercontrol = "inf",p=p)
end

function Gcob3med(;mtime = 432000,meval = -1, p = PDEconstructmeso())
    return counterglobal(;mtime = mtime, meval=meval, alg=:LN_COBYLA, ntarg = 3,countercontrol = "med",p=p)
end

function Gcob3inf(;mtime = 432000,meval = -1, p = PDEconstructmeso())
    return counterglobal(;mtime = mtime, meval=meval, alg=:LN_COBYLA, ntarg = 3,countercontrol = "inf",p=p)
end  

function Gcob4med(;mtime = 432000,meval = -1, p = PDEconstructmeso())
    return counterglobal(;mtime = mtime, meval=meval, alg=:LN_COBYLA, ntarg = 4,countercontrol = "med",p=p)
end
     
function Gcob4inf(;mtime = 432000,meval = -1, p = PDEconstructmeso())
    return counterglobal(;mtime = mtime, meval=meval, alg=:LN_COBYLA, ntarg = 4,countercontrol = "inf",p=p)
end

#local first with meso and 1.5 day=130000, then with fine and ... days

function Lbob3med(x0; mtime =130000,meval = -1, p = PDEconstructmeso())
    return counterlocal(x0;mtime = mtime, meval=meval, alg=:LN_BOBYQA, ntarg = 3,countercontrol = "med",p=p)
end

function Lbob3inf(x0; mtime = 130000,meval = -1, p = PDEconstructmeso())
    return counterlocal(x0;mtime = mtime, meval=meval, alg=:LN_BOBYQA, ntarg = 3,countercontrol = "inf",p=p)
end

function Lbob4med(x0; mtime =130000,meval = -1, p = PDEconstructmeso())
    return counterlocal(x0;mtime = mtime, meval=meval, alg=:LN_BOBYQA, ntarg = 4,countercontrol = "med",p=p)
end

function Lbob4inf(x0; mtime =130000,meval = -1, p = PDEconstructmeso())
    return counterlocal(x0;mtime = mtime, meval=meval, alg=:LN_BOBYQA, ntarg = 4,countercontrol = "inf",p=p)
end

function Lcob3med(x0; mtime =130000,meval = -1, p = PDEconstructmeso())
    return counterlocal(x0;mtime = mtime, meval=meval, alg=:LN_COBYLA, ntarg = 3,countercontrol = "med",p=p)
end

function Lcob3inf(x0; mtime =130000,meval = -1, p = PDEconstructmeso())
    return counterlocal(x0;mtime = mtime, meval=meval, alg=:LN_COBYLA, ntarg = 3,countercontrol = "inf",p=p)
end
        
function Lcob4med(x0; mtime =130000,meval = -1, p = PDEconstructmeso())
    return counterlocal(x0;mtime = mtime, meval=meval, alg=:LN_COBYLA, ntarg = 4,countercontrol = "med",p=p)
end
        
function Lcob4inf(x0; mtime = 130000,meval = -1, p = PDEconstructmeso())
    return counterlocal(x0;mtime = mtime, meval=meval, alg=:LN_COBYLA, ntarg = 4,countercontrol = "inf",p=p)
end
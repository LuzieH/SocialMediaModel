function counterglobal(;mtime = 800000, meval=-1, alg=:LN_COBYLA, ntarg = 3,countercontrol = "med", p = PDEconstructmeso())
    return solveopt(; p = p, q= parametersstronginf(),r=parameterscontrol(ntarg=ntarg,maximize="counter_follower",start="zero"), alg=alg, mtime = mtime, meval=meval, ftol_rel=1e-4, xtol_rel = 1e-2, multistart =true,countercontrol=countercontrol,stubborntarget=[1.5 1.5],x0 = rand(Uniform(-2,2),2*ntarg))
end

function counterlocal(x0; mtime = 800000, meval=-1, alg=:LN_COBYLA, ntarg = 3, countercontrol = "med", p = PDEconstructmeso())
    return solveopt(;p=p, q= parametersstronginf(),r=parameterscontrol(ntarg=ntarg,maximize="counter_follower",start="zero"), alg=alg, mtime = mtime, meval=meval, ftol_rel=0., xtol_rel = 0., multistart =false,countercontrol=countercontrol,stubborntarget=[1.5 1.5],x0 = x0)
end

function Gcobmed(;mtime = 432000,meval = -1, p = PDEconstructmeso())
    return counterglobal(;mtime = mtime, meval=meval, alg=:LN_COBYLA, ntarg = 4,countercontrol = "med",p=p)
end

function Lcobmed(x0; mtime =130000,meval = -1, p = PDEconstructmeso())
    return counterlocal(x0;mtime = mtime, meval=meval, alg=:LN_COBYLA, ntarg = 4,countercontrol = "med",p=p)
end
         
function Gcobinf(;mtime = 432000,meval = -1, p = PDEconstructmeso())
    return counterglobal(;mtime = mtime, meval=meval, alg=:LN_COBYLA, ntarg = 4,countercontrol = "inf",p=p)
end

function Lcobinf(x0; mtime = 130000,meval = -1, p = PDEconstructmeso())
    return counterlocal(x0;mtime = mtime, meval=meval, alg=:LN_COBYLA, ntarg = 4,countercontrol = "inf",p=p)
end

""" final experiments to find strategies of media and influencer counteraction """

# global 10 days = 864000, local 3 days=259200
function GLmed(;mtime = 864000, mtime2 = 259200)
    maxf,maxx, ret, numevals, x_list,followersum_list, penalty_list = Gcobmed(;mtime = mtime)
    maxf2,maxx2, ret2, numevals2, x_list2,followersum_list2,penalty_list2 = Lcobmed(maxx; mtime =mtime2)
    return maxf,maxx, ret, numevals, x_list,followersum_list, penalty_list, maxf2, maxx2, ret2, numevals2, x_list2,followersum_list2, penalty_list2
end

function GLinf(;mtime = 864000, mtime2 = 259200)
    maxf,maxx, ret, numevals, x_list,followersum_list, penalty_list = Gcobinf(;mtime = mtime)
    maxf2,maxx2, ret2, numevals2, x_list2,followersum_list2, penalty_list2 = Lcobinf(maxx; mtime =mtime2)
    return maxf,maxx, ret, numevals, x_list,followersum_list, penalty_list, maxf2, maxx2, ret2, numevals2, x_list2,followersum_list2, penalty_list2
end


""" simulations and plots of different strategies """

function infrightcorner()
    followersum, speedpenalty, sols, Ps = solvefixedtargets([1.5 1.5];  p=PDEconstruct(),q=parametersstronginf(), r=parameterscontrol(ntarg = 1 ,start="zero"), scenario="control",countercontrol = "no",stubborntarget=[1.5 1.5])
    plotsnapshots(sols, Ps, [5. 7. 8.5 10.]; scenario="infmax",followercount=true)
    return followersum, speedpenalty
end

function infcounteraction(;targets =  [ 0.5604084588036462
    0.5455196695338128
    0.642358069258012
    0.6527252366140318
    0.7454806015269124
    0.734952127347543
    0.8060867535176255
    0.7718781212141405])

    followersum, speedpenalty, sols, Ps = solvefixedtargets(targets;  p=PDEconstruct(),q=parametersstronginf(), r=parameterscontrol(ntarg = 4 ,start="zero"), scenario="control",countercontrol = "inf",stubborntarget=[1.5 1.5])
    plotsnapshots(sols, Ps, [5. 6. 7.5 10.]; scenario="counterinf",followercount=true)
    return followersum, speedpenalty
end

function medcounteraction(;targets = [  1.8106679226546316
    -0.9805088659711532
     1.7642225907551847
    -1.7710551914852335
     1.6338751208236737
    -1.8358368944346843
     1.5683845474974514
    -1.7962804090565507])

    followersum, speedpenalty, sols, Ps = solvefixedtargets(targets;  p=PDEconstruct(),q=parametersstronginf(), r=parameterscontrol(ntarg = 4 ,start="zero"), scenario="control",countercontrol = "med",stubborntarget=[1.5 1.5])
    plotsnapshots(sols, Ps, [5. 6. 7.5 10.];  scenario="countermed",followercount=true)
    return followersum, speedpenalty
end
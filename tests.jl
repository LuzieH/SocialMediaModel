function testabm()
    xs, xinfs, infs, meds, state, (p,q) = ABMsolveplot(;NT = 10, ts = [1 2],  p = ABMconstruct(), q=parameters(), scenario="4inf")
    plotfollowernumbers(xinfs,state,(p,q);scenario=scenario)
    ABMgifsingle(xs, xinfs, state, infs, meds, (p,q))
    return nothing
end

function testpde()
    sol, (p,q) = solveplot(; tmax=0.1, ts = [0. 0.1], alg=nothing, scenario="4inf", p = PDEconstructcoarse(), q= parameters())
    gifsingle(sol,(p,q))
    return nothing
end

function testensemble()
    runensembles(2; NT=10, tmax=0.1, savepoints = 2, q=parameters(),sigma=0.1)
    return nothing
end

function testcontrol()
    solvefixedtargets([0.3 0.3];  p=PDEconstructcoarse(),q=parametersstronginf(), r=parameterscontrol(ntarg = 1 ,start="zero",Tmax = 0.2,tequil = 0.2), scenario="control",countercontrol = "no",stubborntarget=[1.5 1.5])
    return nothing
end

function testoptimization()
    GLmed(;mtime = 60, mtime2 = 60)
    GLinf(;mtime = 60, mtime2 = 60)
    return nothing
end

function testabm(init = "4inf")
    println("Testing ABM")
    xs, xinfs, infs, meds, state, (p,q) = ABMsolveplot(;NT = 10, ts = [1 2],  p = ABMconstruct(), q=parameters(), init=init,save=false)
    ABMplotfollowernumbers(xinfs,state,(p,q);name=init,save=false)
    ABMgifsingle(xs, xinfs, state, infs, meds, (p,q),save=false)
    return nothing
end

function testpde(init="4inf")
    println("Testing PDE")
    sol, (p,q) = PDEsolveplot(; tmax=0.1, ts = [0. 0.1], alg=nothing, init=init, p = PDEconstructcoarse(), q= parameters(),save=false)
    PDEgifsingle(sol,(p,q),save=false)
    return nothing
end

function testensemble()
    println("Testing Ensemble Simulation")
    runensembles(2; NT=10, tmax=0.1, savepoints = 2, q=parameters(),sigma=0.1,save=false)
    return nothing
end

function testcontrol()
    println("Testing Control Scenario")
    PDEsolvefixedtargets([0.3 0.3];  p=PDEconstructcoarse(),q=parametersstronginf(), r=parameterscontrol(ntarg = 1 ,start="zero",Tmax = 0.2,tequil = 0.2), init="uniform",countercontrol = "no",stubborntarget=[1.5 1.5])
    
    return nothing
end

function testoptimization()
    println("Testing Optimization of Countercontrol")
    PDEsolveopt(; p = PDEconstructcoarse(), q= parametersstronginf(),r=parameterscontrol(ntarg=1,tequil=0.2, Tmax=0.2), mtime = 10, meval=-1,multistart =true,countercontrol="med",stubborntarget=[1.5 1.5],x0 = rand(Uniform(-2,2),2),save=false)
    return nothing
end

function runtests()
    testabm("4inf")
    testpde("4inf")
    testensemble()
    testcontrol()
    testoptimization()
end
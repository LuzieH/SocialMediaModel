
function paperfigures()
    
    # single ABM realization
    ABMsolveplot();

    # for appendix
    ABMsolveplot(init="4infstrongc",q=parameters(a=1,b=1,c=3))
    ABMsolveplot(init="4infstronga",q=parameters(a=3,b=1,c=1))
    #further good experiments for appendix
    #ABMsolveplot(init="4infweakgamma2.5",q=parameters(frictionI=2.5))
    #ABMsolveplot(init="4infweakgamma1",q=parameters(frictionI=1))
    #ABMsolveplot(init="4infweakgamma2.5noise",q=parameters(sigmahat=1,frictionI=2.5),seed=3)

    # ensemble
    @load string("data/pde_ensemble_4inf.jld2") us zs ys p q
    @load string("data/abm_ensemble_4inf.jld2") us2 zs2 ys2 p2 q2 
    plotensemblesnapshots(us, zs, ys, (p,q), us2, zs2, ys2, (p2,q), 2.; name="4inf",save=save)

    # optimizations
    infrightcorner();
    infcounteraction();
    medcounteraction();
end
using NLopt
using Distributions

""" prepare equilibration of system """
function prep(tequil = 5.; p = PDEconstruct(), q= parametersstronginf(), r=parameterscontrol(), scenario = "control",  savedt=0.05, atol = 1e-6, rtol = 1e-3,dtmin = 0.001,countercontrol="no",returnsol=false)
    uzy0,controlled_inf,controlled_med,q = constructinitial(scenario,(p,q))
    q1 = (; q..., controlled_inf,controlled_med)
    start = r.start
    # Solve the ODE
    prob1 = ODEProblem(f,uzy0,(0.0,tequil),(p,q1))
    @time sol1 = DifferentialEquations.solve(prob1, dtmin = dtmin, force_dtmin = true, nothing, saveat = 0:savedt:tequil,save_start=true, abstol = atol, reltol = rtol)
    Jold = q1.J

    #add new influencer for control
    u,z,y = sol2uyz(sol1,tequil)
    if countercontrol == "no"
        q2 = merge(q1, (;J=q1.J+1,controlled_inf = [q1.controlled_inf..., 3]))
        (;N_x, N_y, grid_points,N) = p
        J= q2.J
    elseif countercontrol == "inf" #add one influencer to control and one to stubbornly go into corner
        q2 = merge(q1, (;J=q1.J+2,controlled_inf = [q1.controlled_inf..., 1,3]))
        (;N_x, N_y, grid_points,N) = p
        J= q2.J        
    elseif countercontrol =="med"
        q2 = merge(q1, (;J=q1.J+1,controlled_inf = [q1.controlled_inf..., 3], controlled_med = [0, 1]))
        (;N_x, N_y, grid_points,N) = p
        J= q2.J        
    end

    u2 = zeros(N_x, N_y, 2, J)
    u2[:,:,:,1:Jold] = u
    if start =="influencer"
        startlocation = y[:,1]' #position of influencer in right corner
    elseif start =="zero"
        startlocation =[0 0]
    elseif start=="mean"
        startlocation = 1/sum(u2,dims=(1,2,4))[1,1,2,1] * reshape(sum(u2, dims=4)[:,:,2,:],1,N)*grid_points
    end
    y2 = zeros(2, J)
    if countercontrol == "no"
        y2[:,1:J-1] = y
        #TODO maybe make constant speed such that starts and ends within 5. time steps
        y2[:,J] = startlocation 
        startlocationmed = nothing 
    elseif countercontrol == "inf" #add one influencer to control and one to stubbornly go into corner
        y2[:,1:J-2] = y
        #TODO maybe make constant speed such that starts and ends within 5. time steps
        y2[:,J-1] = startlocation 
        y2[:,J] = startlocation   
        startlocationmed = nothing  
    elseif countercontrol =="med"
        y2[:,1:J-1] = y
        #TODO maybe make constant speed such that starts and ends within 5. time steps
        y2[:,J] = startlocation     
        startlocationmed = z[:,2]'
    end   

    uzy0 = ArrayPartition(u2,z,y2)
    if returnsol==false
        return uzy0, startlocation, (p,q2), startlocationmed
    else
        return uzy0, startlocation, (p,q2), startlocationmed, sol1, (p,q1)
    end
end
 
""" solve system with controled influencers/media to given fixed targets """

function solvefixedtargetsfast(ts, targets, startlocation, (p, q), r, uzy0; savedt=0.05, atol = 1e-6, rtol = 1e-3,dtmin = 0.001,countercontrol="no", stubborntarget=nothing, startlocationmed = nothing,returnsol=false)
    n_targets = size(targets,2)
    followersum=0.
    speedpenalty=0.
    J = q.J
    if countercontrol !="no"
        stubbornspeed = norm(stubborntarget - startlocation)/ts[end]
        q= merge(q, (;controltarget2 = stubborntarget, controlspeed2 = stubbornspeed))
    end

    if countercontrol =="med"
        startlocation = startlocationmed
    end

    if returnsol!=false
        Ps = Any[]
        sols = Any[]
    end
    
    for n in 1:n_targets
        target = targets[:,n]'
        Dt = ts[n+1]-ts[n]

        speed = norm(target - startlocation)/Dt
        #restrict speedbound?
        q1 = merge(q, (;controltarget = target, controlspeed = speed))

        # solve ODE with added influencer
        prob = ODEProblem(f,uzy0,(0.0,Dt),(p,q1))
        @time sol = DifferentialEquations.solve(prob, dtmin = dtmin, force_dtmin = true, saveat = 0:savedt:Dt ,save_start=true, abstol = atol, reltol = rtol)
        #check stiff solver solution, how often is solution saved? 
        followersum += sum([sum(sol(t).x[1][:,:,:,J])/sum(sol(t).x[1]) for t in sol.t])*savedt
        speedpenalty +=speed^2*Dt
        #prepare next simulation
        startlocation = target
        uzy0 = sol(Dt)
        if returnsol!=false
            push!(Ps, (p,q1))
            push!(sols,sol)
        end
    end

    speedpenalty = sqrt(speedpenalty) #L2 norm of speed over [0,T]
    if returnsol==false
        return followersum, speedpenalty
    else
        return followersum, speedpenalty, sols, Ps
    end
end

vecof2vecs(in::Vector) = collect(eachrow(reshape(in, :, 2)))

""" prepare and solve system with controlled influencers/media to given fixed targets """
function solvefixedtargets(targets;  p=PDEconstructcoarse(),q=parametersstronginf(), r=parameterscontrol(ntarg = 1 ,start="zero"), scenario="control",countercontrol = "no",stubborntarget=[1.5 1.5])
    #options of countercontrol "med" "no" "inf"
    (; ntarg, Tmax,tequil,dtmin) = r
    ts = [0. Tmax/ntarg*collect(1:ntarg)...]
    targets = reshape(targets,(2,ntarg))  #reshape into correct input format

    uzy0, startlocation, (p,q2), startlocationmed, sol1, (p,q1) = prep(tequil;dtmin=dtmin, p=p, q=q, r=r, scenario=scenario,countercontrol = countercontrol,returnsol=true)
    
    followersum, speedpenalty, sols, Ps= solvefixedtargetsfast(ts, targets, startlocation, (p, q2),r, uzy0, dtmin = dtmin,countercontrol=countercontrol, stubborntarget=stubborntarget, startlocationmed = startlocationmed, returnsol=true)
    prepend!(sols,[sol1])
    prepend!(Ps, [(p,q1)])
    return followersum, speedpenalty, sols, Ps
end




function solveopt(; p = PDEconstructcoarse(), q= parametersstronginf(), r=parameterscontrol(), alg=:LN_COBYLA, mtime = 1000, meval=-1, ftol_rel=1e-4, xtol_rel = 1e-2, x0=zeros(2*r.ntarg),countercontrol="no",stubborntarget=nothing, multistart =false)
        
    (; ntarg, Tmax,tequil,dtmin,start,maximize,alpha) = r
    uzy0, startlocation, (p,q),startlocationmed =  prep(tequil; p=p, q=q, r=r, dtmin = dtmin,countercontrol=countercontrol)

    x_list = Vector{Float64}[]
    followersum_list = Float64[]
    penalty_list = Float64[]
    iter = 0

    ts = [0. Tmax/ntarg*collect(1:ntarg)...]
    

    function Objfun(x::Vector, grad::Vector)
        iter += 1
        println("iteration $iter")

        targets = reshape(x,(2,ntarg))  #reshape into correct input format

        
        followersum, speedpenalty = solvefixedtargetsfast(ts, targets, startlocation, (p, q),r, uzy0, dtmin = dtmin,countercontrol=countercontrol, stubborntarget=stubborntarget, startlocationmed = startlocationmed)
        speedpenalty = alpha*speedpenalty
        push!(x_list, x)
        push!(followersum_list, followersum)
        push!(penalty_list, speedpenalty)
        println("followersum $followersum at x $x")
        println("speedpenalty $speedpenalty")
        GC.gc(true)
        #check_memory()
        if maximize=="follower"
            return followersum - alpha*speedpenalty
        elseif maximize =="counter_follower"
            return -followersum-speedpenalty
        else
            println("This obj is unknown")
        end
    end

    #number of parameters
    ndim = 2*ntarg 
    if multistart==false
        opt = Opt(alg, ndim)#
    else
        opt = Opt(:G_MLSL_LDS, ndim)#
        lopt = Opt(alg, ndim) #local optimizer
        lopt.ftol_rel = ftol_rel
        lopt.xtol_rel = xtol_rel
        lopt.maxtime = 20000
        opt.local_optimizer = lopt
    end
    #lb ub are arrays of length ndim that bound the parameters from below and above
    lb = ones(ntarg*2)*p.domain[1,1]#, zeros(ntarg)...]
    ub = ones(ntarg*2)*p.domain[1,2]#, Tmax*ones(ntarg)...]
    opt.lower_bounds = lb::Union{AbstractVector,Real} 
    opt.upper_bounds = ub::Union{AbstractVector,Real}

    
    opt.max_objective = Objfun

    #stopping criteria
    opt.maxtime = mtime
    opt.maxeval =  meval
    opt.ftol_rel = ftol_rel
    opt.xtol_rel = xtol_rel

    (maxf,maxx,ret) = optimize(opt, x0)
    numevals = opt.numevals # the number of function evaluations
    println("got $maxf at $maxx after $numevals iterations (returned $ret)")
    @save string("data/opt_strategy",string(ntarg),string(alg),maximize,start,string(multistart),countercontrol,".jld2") maxf maxx ret numevals x_list followersum_list penalty_list 
    return maxf,maxx, ret, numevals, x_list,followersum_list, penalty_list
end 


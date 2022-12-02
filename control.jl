using NLopt

function parameters_optcont(;ntarg=3, speedbound = 10, Tmax = 5,tequil = 5.)
    r = (; ntarg, speedbound, Tmax,tequil)
    return r
end


#prepare equilibration of density
function prep(tequil = 5.; p = PDEconstruct(), q= parameters_control(),scenario = "optimalcontrol",  savedt=0.05, atol = 1e-6, rtol = 1e-3)
    uzy0, _, controlled,q = constructinitial(scenario,(p,q))
    q1 = (; q..., controlled)

    # Solve the ODE
    prob1 = ODEProblem(f,uzy0,(0.0,tequil),(p,q1))
    @time sol1 = DifferentialEquations.solve(prob1, nothing, saveat = 0:savedt:tequil,save_start=true, abstol = atol, reltol = rtol)


    #add new influencer
    u,z,y = sol2uyz(sol1,tequil)
    q2 = merge(q1, (;J=q1.J+1,controlled = [q1.controlled..., 1]))
    (;N_x, N_y, grid_points,N, dV) = p
    J= q2.J


    u2 = zeros(N_x, N_y, 2, J)
    u2[:,:,:,1:J-1] = u
    startlocation = y[:,1]' #position of influencer in right corner
    # 1/sum(u2,dims=(1,2,4))[1,1,2,1] * reshape(sum(u2, dims=4)[:,:,2,:],1,N)*grid_points #[0 0]
    y2 = zeros(2, J)
    y2[:,1:J-1] = y
    #TODO maybe make constant speed such that starts and ends within 5. time steps
    y2[:,J] = startlocation
    uzy0 = ArrayPartition(u2,z,y2)

    return uzy0, startlocation, (p,q2)
end


function solvefixedtargetsfast(ts, targets, startlocation, (p, q), uzy0; savedt=0.05, atol = 1e-6, rtol = 1e-3)
    bound = Int(round((p.N_x)*0.7))
    n_targets = size(targets,2)
    followersum= 0.
    masscorner = 0.
    dV=p.dV
    J = q.J

    for n in 1:n_targets
        target = targets[:,n]'
        Dt = ts[n+1]-ts[n]
        speed = norm(target - startlocation)/Dt
        q1 = merge(q, (;controltarget = target, controlspeed = speed))

        # solve ODE with added influencer
        prob = ODEProblem(f,uzy0,(0.0,Dt),(p,q1))
        @time sol = DifferentialEquations.solve(prob, nothing,  saveat = 0:savedt:Dt ,save_start=true, abstol = atol, reltol = rtol)


        followersum += sum([sum(sol(t).x[1][:,:,:,J]) for t in sol.t])*dV*savedt
        masscorner += sum([sum(sol(t).x[1][bound:end,bound:end,:,:]) for t in sol.t]) *dV*savedt
        #prepare next simulation
        startlocation = target
        uzy0 = sol(Dt)
    end

    return followersum, masscorner
end

vecof2vecs(in::Vector) = collect(eachrow(reshape(in, :, 2)))


function solveopt(; p = PDEconstructcoarse(), q= parameters_control(), r=parameters_optcont(), alg=:LN_COBYLA, mtime = 1000, meval=-1)

    iter = 0
    (; ntarg, speedbound, Tmax,tequil) = r
    uzy0, startlocation, (p,q) =  prep(tequil; p=p, q=q)

    x_list = Vector{Float64}[]
    followersum_list = Float64[]
    cornersum_list = Float64[]

    function constrainttime(result::Vector, x::Vector, grad::Matrix)
        for j in 1:ntarg-2
            result[j] =  x[2*ntarg+j]-x[2*ntarg+j+1] #has to be <=0
        end
    end

    function constraintspeed(result::Vector, x::Vector, grad::Matrix)
        xy = reshape(x[1:2*ntarg], (2,ntarg))
        xs = [startlocation[1], xy[1,:]...]
        ys = [startlocation[2], xy[2,:]...]
        ts = [0., x[2*ntarg+1:end]...,  Tmax]
        for j in 1:ntarg
            result[j] = (sqrt((xs[j+1]-xs[j])^2 + (ys[j+1]-ys[j])^2))/(ts[j+1]-ts[j]) - speedbound
        end
    end

    function followers(x::Vector, grad::Vector)
        iter += 1
        println("iteration $iter")

        targets = reshape(x[1:2*ntarg],(2,ntarg))  #reshape into correct input format
        ts = [0., x[2*ntarg+1:end]...,Tmax]

        followersum, masscorner = solvefixedtargetsfast(ts, targets, startlocation, (p, q), uzy0)
        push!(x_list, x)
        push!(followersum_list, followersum)
        push!(cornersum_list, masscorner)

        GC.gc()
        GC.gc(true)
        check_memory()
        #GC.gc(true)
        #GC.gc(true)

        return followersum
    end

    #number of parameters
    ndim = 3*ntarg-1 #2ntarg for targets and ntarg-1 for time points
    opt = Opt(alg, ndim)
    #lb ub are arrays of length ndim that bound the parameters from below and above
    lb = [ones(ntarg*2)*p.domain[1,1]..., zeros(ntarg-1)...]
    ub = [ones(ntarg*2)*p.domain[1,2]..., Tmax*ones(ntarg-1)...]
    opt.lower_bounds = lb::Union{AbstractVector,Real}
    opt.upper_bounds = ub::Union{AbstractVector,Real}


    opt.max_objective = followers
    inequality_constraint!(opt, constrainttime, 1e-8*ones(ntarg-2)) #constraint time points to be ordered
    inequality_constraint!(opt, constraintspeed,  1e-8*ones(ntarg)) #constrain speed

    #stopping criteria
    opt.maxtime = mtime
    opt.maxeval =  meval
    x0 = [zeros(2*ntarg)..., collect(1:ntarg-1)*(Tmax/(ntarg))...]
    (maxf,maxx,ret) = optimize(opt, x0)
    numevals = opt.numevals # the number of function evaluations
    println("got $maxf at $maxx after $numevals iterations (returned $ret)")
    @save string("data/opt_strategy.jld2") maxf maxx ret numevals x_list followersum_list cornersum_list
    return maxf,maxx,ret, numevals, x_list,followersum_list,cornersum_list
end


function check_memory(limit=32, verbose = true)
    gb = open("/proc/$(getpid())/statm") do io
        parse(Int, split(read(io, String))[1]) * 4096 / 1024 / 1024 / 1024  # in GB
    end
    if verbose
        run(`ps -p $(getpid()) -o pid,comm,vsize,rss,size`);
        @show gb
    end
    if gb > limit # GB
        @warn "Exceeded memory bounds"
        throw(ErrorException("Exceeded memory bounds"))
    end
end

using LinearAlgebra
using StatsBase
 
using JLD2
 



function ABMinfinitialconditions((p,q))
    (; n, J) = q

    # agent opinions
    x = rand(n,2).*4 .-2 
    # media opinions
    media =[-1. -1.; 1. 1.]
    # agents following different influencers
    xI = Any[]
    push!(xI,  intersect(findall(x-> x>0,x[:,1]), findall(x-> x>0,x[:,2])))
    push!(xI,  intersect(findall(x-> x<=0,x[:,1]), findall(x-> x>0,x[:,2])))
    push!(xI,  intersect(findall(x-> x>0,x[:,1]), findall(x-> x<=0,x[:,2])) )
    push!(xI,  intersect(findall(x-> x<=0,x[:,1]), findall(x-> x<=0,x[:,2]))) 

    #follower network
    fol=zeros(n, J)
    for i in 1:J
        fol[xI[i],i] .=1
    end

    #initial opinions of influencers
    inf = zeros(J,2)
    for i in 1:J
        inf[i,:] = sum(x[xI[i],:],dims=1)/size(xI[i],1)#+ rand(2)'*0.5
    end

    # initialization of political attitude state of all agents
    state = (rand(n).>0.5)*2 .-1
    xM = Any[]
    push!(xM, findall(x->x==-1, state))
    push!(xM, findall(x->x==1, state))

    #initialization of interaction network between agents 
    Net = ones(n,n)  # everyone connected to everyone including self-connections

    # initial number of influencers of different attituted that follow the different influencers
    counts = zeros(J,2)
    for j in 1:J
        for i in 1:2
            counts[j,i] = size(intersect(xI[j], xM[i]),1)
        end
    end

    return x, media, inf, fol, state, Net, counts 
end

function ABMnoinfinitialconditions((p,q))
    (; n) = q

    # agent opinions
    x = rand(n,2).*4 .-2 
    # media opinions
    media =[-1. -1.; 1. 1.]

    #follower network
    fol=ones(n, 1)

    #initial opinions of influencers
    inf = zeros(1,2)
    inf[1,:] = [0 0]

    # initialization of political attitude state of all agents
    state = (rand(n).>0.5)*2 .-1
    xM = Any[]
    push!(xM, findall(x->x==-1, state))
    push!(xM, findall(x->x==1, state))

    #initialization of interaction network between agents 
    Net = ones(n,n)  # everyone connected to everyone including self-connections

    # initial number of influencers of different attituted that follow the different influencers
    counts = [size(xM[1],1) size(xM[2],1)]
    return x, media, inf, fol, state, Net, counts 
end



function attraction(x, Net)
    n = size(Net, 1)
    force = zeros(n,2)
    orderparameter=0.
    for j in 1:n
        L=findall(x->x==1, Net[j,:])
        if isempty(L)
            force[j,:]=[0 0]
        else
            fi = [0 0]
            w_sum=0
            for i in 1:length(L)
                d = x[L[i],:]-x[j,:]
                w = exp(-int_decay*sqrt(d[1]^2 +d[2]^2))
                fi = fi + w*d'
                w_sum = w_sum+w
            end
            force[j,:] = fi/w_sum
        end
        orderparameter+= w_sum
    
    end
    orderparameter = 1/n^2 * orderparameter
    return force, orderparameter
end

function influence(x,media,inf,fol,state,(p,q))
    (; n, b, c, J) = q
    force1 =zeros(size(x))
    force2 =zeros(size(x))
    ord_med = 0.
    ord_inf = 0. 
    for j in 1:n
        if state[j]==1
            force1[j,:] = media[2,:]-x[j,:]
            ord_med += exp(-order_decay*norm(force1[j,:]))
        else  
            force1[j,:] = media[1,:]-x[j,:]
            ord_med += exp(-order_decay*norm(force1[j,:]))
        end

        for k in 1:J
            if fol[j,k]==1
                force2[j,:]=inf[k,:] -x[j,:]
                ord_inf += exp(-order_decay*norm(force2[j,:]))
            end
        end
    end
    force = c*force1 + b*force2
    return force, ord_med/n, ord_inf/n
end

function changeinfluencer(state,x,fol,inf,(p,q))
    (;   eta, n, J) =q
    dt = p.dt
    
    theta =0.1 #threshold for discrete g-function
    
    # compute happiness = fraction of followers with same state
    fraction = zeros(J)
    for i=1:J
        fraction[i]= sum(fol[:,i].*state)/sum(fol[:,i])
    end
    
    #ompute distance of followers to influencers
    dist = zeros(n, J)
    for i=1:J
        for j=1:n
            d = x[j,:]-inf[i,:]
            dist[j,i]= exp(-sqrt(d[1]^2+d[2]^2))
        end
    end
    
    # compute attractiveness of influencer for followers
    attractive = zeros(n, J)
    for j=1:n
        for i=1:J
            g2 = state[j]*fraction[i]
            if g2<theta 
                g2=theta
            end
            attractive[j,i]= eta * dist[j,i]*g2 
        end

        r=rand()
        lambda = sum(attractive[j,:])
        if r<1-exp(-lambda*dt) 
            p = attractive[j,:]/lambda
            r2=rand()
            k=1
            while sum(p[1:k])<r2
                k=k+1
            end
            fol[j,:]=zeros(J)      
            fol[j,k]=1
        end
    end

    return fol
end



function ABMsolve(NT = 100;  p = ABMconstruct(), q=parameters(), scenario="4inf")
    (; dt, domain) = p
    (;n, n_media, J, sigma, sigmahat, sigmatilde, a, b,c, frictionI, frictionM) =q

    if scenario=="4inf"
        x, media, inf, fol, state,Net,counts  = ABMinfinitialconditions((p,q))
    elseif scenario=="noinf"
        x, media, inf, fol, state,Net,counts  =  ABMnoinfinitialconditions((p,q))
        J=1
        b=0
        eta=0
        q=(;q..., J,b,eta)
    end



    xs = [x] #initial condition
    infs = [inf]
    meds = [media]
    xinfs = [fol * collect(1:J)]
    orderparameters = Any[]
    ord_infs = Any[]
    ord_meds = Any[]
    for k in 2:NT+1
        xold = x
        # opinions change due to opinions of friends, influencers and media
        attforce, orderparameter = attraction(xold,Net)
        leader_force, ord_med, ord_inf = influence(xold,media,inf,fol,state,(p,q))
        push!(orderparameters, orderparameter)
        push!(ord_infs, ord_inf)
        push!(ord_meds, ord_med)
        force = a * attforce + leader_force
        x = xold + dt*force + sqrt(dt*sigma)*randn(n,2); 
        # dont allow agents to escape domain
        ind1 = findall(x->x>domain[1,2],x)
        ind2 = findall(x->x<domain[1,1],x)
        x[ind1] .= 2
        x[ind2] .= -2


        # influencer opinions adapt slowly to opinions of followers with friction
        # depending on number of followers
        masscenter=zeros(J,2)
        for i in 1:J
            if sum(fol[:,i])>0 
                masscenter[i,:] =sum(fol[:,i] .* x, dims = 1) /sum(fol[:,i])
                inf[i,:] =  inf[i,:]  + dt/frictionI * (masscenter[i,:]-inf[i,:]) + 1/frictionI*sqrt(dt*sigmahat)*randn(2,1)
            end
        end
        
        # media opinions change very slowly based on opinions of followers with friction
        # depending on number of followers
        masscenter=zeros(n_media,2)
        states = [-1, 1]
        for i in 1:n_media
            x_M = findall(x->x==states[i], state)
            masscenter[i,:] = sum(xold[x_M,:], dims=1)/size(x_M,1)
            media[i,:] = media[i,:]  + dt/frictionM * (masscenter[i,:] -media[i,:]) + 1/frictionM * sqrt(dt*sigmatilde)*randn(2,1)
        end
        
        # individual may jump from one influencer to another
        # jumps according to rate model
        fol = changeinfluencer(state,xold,fol,inf,(p,q))

        xs = push!(xs,copy(x))
        infs = push!(infs, copy(inf))
        meds = push!(meds, copy(media))
        xinfs = push!(xinfs, copy(fol * collect(1:J)))
    end
    attforce, orderparameter = attraction(x,Net)
    leader_force, ord_med, ord_inf = influence(x,media,inf,fol,state,(p,q))
    push!(orderparameters, orderparameter)
    push!(ord_infs, ord_inf)
    push!(ord_meds, ord_med)
    return xs, xinfs, infs, meds, state, (p,q), counts,orderparameters, ord_infs, ord_meds

end


 
function sumgaussian(x, centers;sigma=0.1)
    output = 0
    for i in 1:size(centers,1)
        output = output +  gaussian(x, centers[i,:],sigma=sigma) 
    end
    return output
end 


function ABMsolveplot(NT = 100;  p = ABMconstruct(), q=parameters(), scenario="4inf")
    @time xs, xinfs, infs, meds, state, (p,q), counts,orderparameters, ord_infs, ord_meds = ABMsolve(NT;  p=p, q=q, scenario=scenario)
    (;dt) = p
    ABMplotarray(xs[end], xinfs[end],state,  infs[end], meds[end],  (p,q), dt*(NT-1), scenario=scenario)
    return xs, xinfs, infs, meds, state, (p,q), counts,orderparameters, ord_infs, ord_meds
end




function ABMsolveensemble(NT=200, N=10; savepoints = 4, scenario="4inf", p = ABMconstruct(), q= parameters(),sigma=0.1)
    (; X, Y, domain, dx, dy) = p
    (; n, J) = q
    x_arr = domain[1,1]:dx:domain[1,2]
    y_arr = domain[2,1]:dy:domain[2,2] 
    N_x = size(x_arr,1)
    N_y = size(y_arr,1)
    zs = zeros(2, 2, savepoints, N)
    ys = zeros(2, J, savepoints, N)
    us = zeros(N_x, N_y, 2, J, savepoints, N)
    savetimes = Int.(round.(LinRange(1, NT, savepoints)))
    av_counts = zeros(J,2)
    states = [-1 1]
    ord_infs_mean = zeros(NT+1)
    ord_meds_mean =  zeros(NT+1)
    ord_infs_sqrmean = zeros(NT+1)
    ord_meds_sqrmean =  zeros(NT+1)
    #Threads.@threads 
    for k=1:N
        xs, xinfs, infs, meds, state, _, counts,orderparameters, ord_infs, ord_meds = ABMsolve(NT;  p=p, q=q, scenario=scenario)
        av_counts = av_counts +  counts*(1/N)

        ord_infs_mean += ord_infs*(1/N)
        ord_meds_mean +=  ord_meds*(1/N)    
        ord_infs_sqrmean += ord_infs.^2*(1/N)
        ord_meds_sqrmean +=  ord_meds.^2*(1/N)   


        for m in 1:savepoints
            t = savetimes[m]
            x = xs[t]
            xinf = xinfs[t]
            inf = infs[t]
            media = meds[t]
            for i in 1:2
                for j in 1:J
                    xi  = findall(x-> x==j, xinf)
                    xm = findall(x-> x==states[i], state)
                    choice = intersect(xi, xm)
                    us[:,:,i,j,m,k] = [(1/n)*sumgaussian([X[i,j], Y[i,j]], x[choice,:],sigma=sigma) for i in 1:size(X,1), j in 1:size(X,2)]
                    ys[:,j,m,k] = inf[j,:]
                end
                zs[:,i,m,k] = media[i,:]
            end
        end

    end

    @save string("data/abm_ensemble_",scenario,".jld2") us zs ys ord_infs_mean ord_meds_mean ord_infs_sqrmean ord_meds_sqrmean 
    return us, zs, ys, (p,q), av_counts, ord_infs_mean, ord_meds_mean, ord_infs_sqrmean, ord_meds_sqrmean 
end



function runensembles(N; NT=200, tmax=2., savepoints = 5, q=parameters(),clmax=0.5,sigma=0.1) #0.025 works well with 1000 siulations
    us, zs, ys, (p,q), av_counts, ord_infs_mean, ord_meds_mean, ord_infs_sqrmean, ord_meds_sqrmean = solveensemble(tmax, N;savepoints=savepoints, q=q)
    #plotensemble(us, zs, ys,(p,q), tmax; clmax=clmax)
    us2, zs2, ys2, (p2,q2), av_counts2, ord_infs_mean2, ord_meds_mean2, ord_infs_sqrmean2, ord_meds_sqrmean2 = ABMsolveensemble(NT,N; savepoints=savepoints, q=q,sigma=sigma)
    #ABMplotensemble(us2, zs2, ys2, (p2,q2), NT; clmax=clmax)
    plotensemblesnapshots(us, zs, ys, (p,q), us2, zs2, ys2, (p2,q2), tmax; scenario="4inf")
    return us, zs, ys, (p,q), av_counts, ord_infs_mean, ord_meds_mean, ord_infs_sqrmean, ord_meds_sqrmean, us2, zs2, ys2, (p2,q2), av_counts2, ord_infs_mean2, ord_meds_mean2, ord_infs_sqrmean2, ord_meds_sqrmean2
end

function runensembles_noinf(N; NT=100, tmax=1.,savepoints = 5, q=parameters(),clmax=0.5,sigma=0.1)
    scenario="noinf"
    us, zs, ys, (p,q), av_counts = solveensemble(tmax, N; scenario=scenario, savepoints=savepoints, q=q)
    plotensemble(us, zs, ys, (p,q), tmax; clmax=clmax,scenario=scenario)
    us2, zs2, ys2, (p2,q2), av_counts2 = ABMsolveensemble(NT,N; scenario=scenario, savepoints=savepoints, q=q,sigma=sigma)
    ABMplotensemble(us2, zs2, ys2, (p2,q2), NT; clmax=clmax,scenario=scenario)


    #plot order parameter mean
    N = size(ord_infs_mean2,1) #number of timesteps
    (; J,n) = q2
    (;dt) = p2    
    plot(dt*collect(0:N-1), ord_infs_mean2,label="wrt to influencers (ABM)",size=(95*5,60*5),legend=:bottomright,title="Orderparameter",xlabel="t")
    plot!(dt*collect(0:N-1), ord_meds_mean2,label="wrt to media (ABM)")
    N = size(ord_infs_mean,1)
    plot!(tmax/N*collect(0:N-1), ord_infs_mean,label="wrt to influencers (PDE)",size=(95*5,60*5),legend=:bottomright,title="Orderparameter",xlabel="t")
    plot!(tmax/N*collect(0:N-1), ord_meds_mean,label="wrt to media (PDE)")
    savefig(string("img/ensemble_order_",scenario,".png"))
    savefig(string("img/ensemble_order_",scenario,".pdf"))

    return us, zs, ys, (p,q), av_counts, us2, zs2, ys2, (p2,q2), av_counts2
end



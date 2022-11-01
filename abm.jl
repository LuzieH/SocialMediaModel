using LinearAlgebra
using StatsBase
using Plots
using JLD2

#need to rename all functions to be unique from PDE???
# TODO: adapt to different numbers of influencers and controlled influencers!
# maybe change function g? 
# plot: point cloud of agents and Kernel density estimation


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
    for j in 1:n
        L=findall(x->x==1, Net[j,:])
        if isempty(L)
            force[j,:]=[0 0]
        else
            fi = [0 0]
            w_sum=0
            for i in 1:length(L)
                d = x[L[i],:]-x[j,:]
                w = exp(-sqrt(d[1]^2 +d[2]^2))
                fi = fi + w*d'
                w_sum = w_sum+w
            end
            force[j,:] = fi/w_sum
        end
    
    end
    return force
end

function influence(x,media,inf,fol,state,(p,q))
    (; n, b, c, J) = q
    force1 =zeros(size(x))
    force2 =zeros(size(x))
    for j in 1:n
        if state[j]==1
            force1[j,:] = media[2,:]-x[j,:]
        else  
            force1[j,:] = media[1,:]-x[j,:]
        end

        for k in 1:J
            if fol[j,k]==1
                force2[j,:]=inf[k,:] -x[j,:]
            end
        end
    end
    force = c*force1 + b*force2
    return force
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

        # r=rand()
        # lambda = sum(attractive[j,:]) #total jump rate
        # alpha=-log(1-r)/lambda #random number distributed due to exp(lambda)
        # if dt>alpha 
        #     p = attractive[j,:]/lambda
        #     r2=rand()
        #     k=1
        #     while sum(p[1:k])<r2
        #         k=k+1
        #     end
        #     fol[j,:]=zeros(J)      
        #     fol[j,k]=1
        # end

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
    (;n, n_media, J, sigma, sigmahat, sigmatilde, a, b,c, frictionI, frictionM ) =q

    if scenario=="4inf"
        x, media, inf, fol, state,Net,counts  = ABMinfinitialconditions((p,q))
    elseif scenario=="noinf"
        x, media, inf, fol, state,Net,counts  =  ABMnoinfinitialconditions((p,q))
        J=1
        q=(;q..., J)
    end



    xs = [x] #initial condition
    infs = [inf]
    meds = [media]
    xinfs = [fol * collect(1:J)]

    for k in 2:NT+1
        xold = x
        # opinions change due to opinions of friends, influencers and media
        force = a * attraction(xold,Net) + influence(xold,media,inf,fol,state,(p,q))
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

    return xs, xinfs, infs, meds, state, (p,q), counts

end
 
function sumgaussian(x, centers)
    output = 0
    for i in 1:size(centers,1)
        output = output +  gaussian(x, centers[i,:]) 
    end
    return output
end 

function kdeplot(centers, inf, media, (p,q); title = "",labely ="", labelz ="", scenario="4inf",clim=(-Inf, Inf))
    (;X, Y, domain, dx, dy) = p
    n= q.n

    x_arr = domain[1,1]:dx:domain[1,2]
    y_arr = domain[2,1]:dy:domain[2,2]
    evalkde = [(1/n)*sumgaussian([X[i,j], Y[i,j]], centers) for i in 1:size(X,1), j in 1:size(X,2)]
    #K = kde(x, boundary = ((-2,2),(-2,2)))
    subp = heatmap(x_arr, y_arr, evalkde', c=:berlin, title = title,alpha=0.5, clims=clim)
    #scatter!(subp, centers[:,1], centers[:,2], markercolor=:white,markersize=3, lab = "agents")
    if scenario=="4inf"
        scatter!(subp, inf[1,:], inf[2,:], markercolor=:red,markersize=5, lab=labely)
    end
    scatter!(subp, media[1,:], media[2,:], markercolor=:yellow,markersize=5, lab=labelz)
    return subp
end

function ABMsolveplot(NT = 100;  p = ABMconstruct(), q=parameters(), scenario="4inf")
    @time xs, xinfs, infs, meds, state, (p,q), counts = ABMsolve(NT;  p=p, q=q, scenario=scenario)
    (;dt) = p
    ABMplotarray(xs[end], xinfs[end],state,  infs[end], meds[end],  (p,q), dt*NT, scenario=scenario)
    return xs, xinfs, infs, meds, state, (p,q), counts
end

function ABMplotarray(x, xinf, state, inf, media,  (p,q), t; save = true, scenario="4inf")
    (; J,n) =q
    plot_array = Any[]  
    z_labels = ["z₋₁","z₁" ]
    y_labels = ["y₁", "y₂", "y₃", "y₄"]
    dens_labels = ["ρ₋₁₁" "ρ₋₁₂" "ρ₋₁₃" "ρ₋₁₄";"ρ₁₁" "ρ₁₂" "ρ₁₃" "ρ₁₄"]
    states = [-1 1]
    for j in 1:J    
        for i in 1:2
            xi  = findall(x-> x==j, xinf)
            xm = findall(x-> x==states[i], state)
            choice = intersect(xi, xm)

            title = string(dens_labels[i,j],"(", string(round(t, digits=2)), "), prop = ", string(length(choice)/n)) 
            # make a plot and add it to the plot_array
            push!(plot_array, kdeplot(x[choice,:], inf[j,:], media[i,:], (p,q); title = title,labely = y_labels[j], labelz = z_labels[i], scenario=scenario))
        end
    end
    plot(plot_array..., layout=(J,2),size=(1000,J*250)) #|> display

    if save==true
        savefig(string("img/abm_array_",scenario,".png"))
    end
end

function ABMplotsnapshots(xs, infs, meds, state, (p,q), ts; save = true, scenario="4inf")
    (; J) = q
    (;dt) = p
    n_snapshots = length(ts)
    plot_array = Any[]  
    for s in 1:n_snapshots
        t=ts[s]-dt
        x=xs[t]
        inf=infs[t]
        media=meds[t]
        title =string("t = ", string(round(t*dt, digits=2)))
        subp= kdeplot(x, inf', media', (p,q), scenario=scenario, title = title, clim=(0,0.5))
        x1 = findall(x-> x==-1, state)
        scatter!(subp, x[x1,1], x[x1,2], markercolor=:blue,markersize=3, lab = "media 1")
        x2 = findall(x-> x==1, state)
        scatter!(subp, x[x2,1], x[x2,2], markercolor=:white,markersize=3, lab = "media 2")
        push!(plot_array, subp)
    end    
    plot(plot_array..., layout=(n_snapshots,1),size=(90*6,n_snapshots*50*6))#, left_margin=10mm )


    if save==true
        savefig(string("img/abm_snapshots_",scenario,".png"))
    end
end

function ABMplotsingle(x, inf, media, state, (p,q), t; save = true, scenario="4inf")
    (; J) = q
    title =string("t = ", string(round(t, digits=2)))
    subp= kdeplot(x, inf', media', (p,q), scenario=scenario, title = title)
    x1 = findall(x-> x==-1, state)
    scatter!(subp, x[x1,1], x[x1,2], markercolor=:blue,markersize=3, lab = "media 1")
    x2 = findall(x-> x==1, state)
    scatter!(subp, x[x2,1], x[x2,2], markercolor=:white,markersize=3, lab = "media 2")

    plot(subp) #|> display

    if save==true
        savefig(string("img/abm_single_",scenario,".png"))
    end
end



function ABMgifarray(xs, xinfs, state, infs, meds, (p,q); dN=5, scenario="4inf")
    NT=size(xs,1)
    (;dt) = p
    #cl = (0, maximum(maximum(sol(t).x[1]) for t in 0:dt:tmax)) #limits colorbar
 
    abmgif = @animate for t = 1:dN:NT
        ABMplotarray(xs[t], xinfs[t], state, infs[t], meds[t],  (p,q), t*dt; save = false, scenario=scenario)
    end
    Plots.gif(abmgif, string("img/ABM_array",scenario,".gif"), fps = 10)
end

function ABMgifsingle(xs, xinfs, state, infs, meds, (p,q); dN=5, scenario="4inf")
    NT=size(xs,1)
    (;dt) = p
    #cl = (0, maximum(maximum(sol(t).x[1]) for t in 0:dt:tmax)) #limits colorbar
 
    abmgif = @animate for t = 1:dN:NT
        ABMplotsingle(xs[t], infs[t], meds[t],state, (p,q), t*dt; save = false, scenario=scenario)
    end
    Plots.gif(abmgif, string("img/abm_single_",scenario,".gif"), fps = 10)
end


function ABMsolveensemble(NT=100, N=10; savepoints = 4, scenario="4inf", p = ABMconstruct(), q= parameters())
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
    Threads.@threads for k=1:N
        xs, xinfs, infs, meds, state, _, counts = ABMsolve(NT;  p=p, q=q, scenario=scenario)
        av_counts = av_counts +  counts*(1/N)

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
                    us[:,:,i,j,m,k] = [(1/n)*sumgaussian([X[i,j], Y[i,j]], x[choice,:]) for i in 1:size(X,1), j in 1:size(X,2)]
                    ys[:,j,m,k] = inf[j,:]
                end
                zs[:,i,m,k] = media[i,:]
            end
        end

    end

    @save string("data/abm_ensemble",scenario,".jld2") us zs ys
    return us, zs, ys, (p,q), av_counts 
end

function ABMplotensemble(us, zs, ys, (p,q), NT; clmax = clmax, scenario="4inf")
    (; dt) = p
    plotensemble(us, zs, ys, (p,q), dt*NT; title1 = "img/abm_ensemble", title2 = "img/abm_ensemble_influencer", clmax = clmax,scenario=scenario)
end

function runensembles(N; NT=150, tmax=1.5)
    us, zs, ys, (p,q), av_counts = solveensemble(tmax, N)
    plotensemble(us, zs, ys,(p,q), tmax; clmax=0.5)
    us2, zs2, ys2, (p2,q2), av_counts2 = ABMsolveensemble(NT,N)
    ABMplotensemble(us2, zs2, ys2, (p2,q2), NT; clmax=0.5)
    return us, zs, ys, (p,q), av_counts, us2, zs2, ys2, (p2,q2), av_counts2
end

function runensembles_noinf(N; NT=150, tmax=1.5)
    scenario="noinf"
    us, zs, ys, (p,q), av_counts = solveensemble(tmax, N; scenario=scenario, q= parameters(J=1, b=0, eta=0))
    plotensemble(us, zs, ys, (p,q), tmax; clmax=0.5,scenario=scenario)
    us2, zs2, ys2, (p2,q2), av_counts2 = ABMsolveensemble(NT,N; scenario=scenario, q= parameters(J=1, b=0, eta=0))
    ABMplotensemble(us2, zs2, ys2, (p2,q2), NT; clmax=0.5,scenario=scenario)
    return us, zs, ys, (p,q), av_counts, us2, zs2, ys2, (p2,q2), av_counts2
end
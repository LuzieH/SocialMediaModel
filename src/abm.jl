using LinearAlgebra
using StatsBase
using Random
using JLD2
 

""" create ABM initial conditions """
function ABMinit((p,q))
    (; n, J) = q

    # individuals' opinions
    x = rand(n,2).*4 .-2 
    # media opinions
    media =[-1. -1.; 1. 1.]
    # individuals following different influencers in the different quadrants
    xI = Any[]
    push!(xI,  intersect(findall(x-> x>0,x[:,1]), findall(x-> x>0,x[:,2])))
    push!(xI,  intersect(findall(x-> x<=0,x[:,1]), findall(x-> x>0,x[:,2])))
    push!(xI,  intersect(findall(x-> x>0,x[:,1]), findall(x-> x<=0,x[:,2])) )
    push!(xI,  intersect(findall(x-> x<=0,x[:,1]), findall(x-> x<=0,x[:,2]))) 

    #network between individuals and influencers
    folinf=zeros(n, J)
    for i in 1:J
        folinf[xI[i],i] .=1
    end

    #initial opinions of influencers
    inf = zeros(J,2)
    for i in 1:J
        inf[i,:] = sum(x[xI[i],:],dims=1)/size(xI[i],1)#+ rand(2)'*0.5
    end

    # initial medium per individual, either medium -1 or 1
    state = (rand(n).>0.5)*2 .-1
    xM = Any[]
    push!(xM, findall(x->x==-1, state))
    push!(xM, findall(x->x==1, state))

    #initialization of interaction network between individuals
    Net = ones(n,n)  # everyone connected to everyone including self-connections

    return x, media, inf, folinf, state, Net
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
            for i in eachindex(L)
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

function influence(x,media,inf,folinf,state,(p,q))
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
            if folinf[j,k]==1
                force2[j,:]=inf[k,:] -x[j,:]
            end
        end
    end
    force = c*force1 + b*force2
    return force
end

function changeinfluencer(state,x,folinf,inf,(p,q))
    (; eta, n, J) =q
    dt = p.dt
    
    theta =0.1 # threshold for r-function
    
    # compute structural similiarity = fraction of followers with same state
    fraction = zeros(J)
    for i=1:J
        fraction[i]= sum(folinf[:,i].*state)/sum(folinf[:,i])
    end
    
    # compute distance of followers to influencers
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
            folinf[j,:]=zeros(J)      
            folinf[j,k]=1
        end
    end

    return folinf
end


"""Simulate the ABM """
function ABMsolve(NT = 100;  p = ABMconstruct(), q=parameters(), init="4inf",chosenseed=0)
    Random.seed!(chosenseed)
    (;dt, domain) = p
    (;n, n_media, J, sigma, sigmahat, sigmatilde, a,  frictionI, frictionM) = q

    x, media, inf, folinf, state, Net  = ABMinit((p,q))

    xs = [copy(x)] # list of opinions of individuals in time
    infs = [copy(inf)] # list of opinions of influencers in time
    meds = [copy(media)] # list of opinions of media in time
    xinfs = [copy(folinf) * collect(1:J)] # list of which influencer an individual follows in time

    for k in 2:NT+1
        xold = x

        # opinions change due to interaction with opinions of friends, influencers and media
        attforce = attraction(xold,Net)
        leaderforce = influence(xold,media,inf,folinf,state,(p,q))
        force = a * attforce + leaderforce
        x = xold + dt*force + sqrt(dt)*sigma*randn(n,2); 

        # dont allow agents to escape domain
        ind1 = findall(x->x>domain[1,2],x)
        ind2 = findall(x->x<domain[1,1],x)
        x[ind1] .= 2
        x[ind2] .= -2

        # influencer opinions adapt slowly to opinions of followers with friction
        masscenter=zeros(J,2)
        for i in 1:J
            if sum(folinf[:,i])>0 
                masscenter[i,:] =sum(folinf[:,i] .* x, dims = 1) /sum(folinf[:,i])
                inf[i,:] =  inf[i,:]  + dt/frictionI * (masscenter[i,:]-inf[i,:]) + 1/frictionI*sqrt(dt)*sigmahat*randn(2,1)
            end
        end
        
        # media opinions change very slowly based on opinions of followers with friction
        masscenter=zeros(n_media,2)
        states = [-1, 1]
        for i in 1:n_media
            x_M = findall(x->x==states[i], state)
            masscenter[i,:] = sum(xold[x_M,:], dims=1)/size(x_M,1)
            media[i,:] = media[i,:]  + dt/frictionM * (masscenter[i,:] -media[i,:]) + 1/frictionM * sqrt(dt)*sigmatilde*randn(2,1)
        end
        
        # dont allow influencers and media to escape domain
        ind1 = findall(x->x>domain[1,2],inf)
        ind2 = findall(x->x<domain[1,1],inf)
        inf[ind1] .= 2
        inf[ind2] .= -2
        ind1 = findall(x->x>domain[1,2],media)
        ind2 = findall(x->x<domain[1,1],media)
        media[ind1] .= 2
        media[ind2] .= -2

        # individual may jump from one influencer to another
        # jumps according to rate model
        folinf = changeinfluencer(state,xold,folinf,inf,(p,q))

        xs = push!(xs,copy(x))
        infs = push!(infs, copy(inf))
        meds = push!(meds, copy(media))
        xinfs = push!(xinfs, copy(folinf * collect(1:J)))
    end

    return xs, xinfs, infs, meds, state, (p,q)

end

function ABMsolveplot(;NT = 200, ts = [1 11 51 101 151 201],  p = ABMconstruct(), q=parameters(), init="4inf", save=true,seed=0)
    @time xs, xinfs, infs, meds, state, (p,q) = ABMsolve(NT;  p=p, q=q, init=init,chosenseed=seed)
    ABMplotsnapshots(xs, xinfs, infs, meds, state, (p,q), ts; name=init,save=save)
    ABMplotfollowernumbers(xinfs,state,(p,q), save=save)
    return xs, xinfs, infs, meds, state, (p,q)
end

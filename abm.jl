using LinearAlgebra
using StatsBase
using KernelDensity
using Plots
using JLD2

#need to rename all functions to be unique from PDE???
# TODO: adapt to different numbers of influencers and controlled influencers!
# maybe change function g? 
# plot: point cloud of agents and Kernel density estimation

# Prepare grid for distribution over many simulations
#hGrid = 0.2
#maxNr = fix(4/hGrid)
#gridx=-2:hGrid:2
#gridy=-2:hGrid:2
#gridcentersx = -2+hGrid/2:hGrid:2-hGrid/2
#gridcentersy = -2+hGrid/2:hGrid:2-hGrid/2
#histograms for Individuals, Influencer and Media
#histInd = zeros(2,4,maxNr, maxNr) 
#histInf = zeros(maxNr, maxNr)
#histMed = zeros(maxNr, maxNr)

function ABMconstruct()
   
    # setting model and simulation parameters
    dt=0.01  # simulation stepsize
    n = 128 # number of agents
    n_media = 2
    J = 4
    sigma=0.5 # noise on individual agents 
    sigmahat=0 # noise on influencers
    sigmatilde=0 # noise on media
    a=1. #a=1 in paper, interaction strength between agents 
    b=2. # interaction strength between agents and influencers
    c=1. # interaction strength between agents and media
    frictionI = 25 # friction for influencers
    friction = 100  #friction for media
    eta = 15  #rate constant for changing influencer 
    dx = 0.05
    dy = dx
    domain = [-2 2; -2 2]
    X = [x for x in domain[1,1]:dx:domain[1,2], y in domain[2,1]:dy:domain[2,2]]
    Y = [y for x in domain[1,1]:dx:domain[1,2], y in domain[2,1]:dy:domain[2,2]]
    #grid_points = [vec(X) vec(Y)] 

    p = (; dt, n, n_media, J, sigma, sigmahat, sigmatilde, a, b,c, frictionI, friction, eta, dx, dy, domain, X, Y)

    return p
end

function ABMrandominitialconditions(p)
    (; n, n_media, J) = p

    # agent opinions
    x = rand(n,2).*4 .-2 
    # media opinions
    media =[-1. 1.; -1. 1.]
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
        inf[i,:] = sum(x[xI[i],:],dims=1)/size(xI[i],1)
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

function attraction(x, Net)
    n = size(Net, 1)
    force = zeros(n,2)
    for j in 1:n
        L=findall(x->x==1, Net[j,:]) # TODO: maybe simplify for complete network
        if isempty(L)
            force[j,:]=[0 0]
        else
            fi = [0 0]
            wsum=0
            for i in 1:size(L,1)
                d = x[L[i],:]-x[j,:]
                w = exp(-(d[1]^2 +d[2]^2))
                fi = fi + w*d'
                wsum = wsum+w
            end
            force[j,:] = fi/wsum
        end
    
    end
    return force
end

function influence(x,media,inf,fol,state,p)
    (; n, b, c, J) = p
    force1 =zeros(size(x))
    force2 =zeros(size(x))
    for j in 1:n
        if state[j]==1
            force1[j,:] = media[:,2]-x[j,:]
        else  
            force1[j,:] = media[:,1]-x[j,:]
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

function changeinfluencer(state,x,fol,inf,p)
    (; dt, eta, n, J) = p
    
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
            dist[j,i]= exp(-(d[1]^2+d[2]^2))
        end
    end
    
    # compute attractiveness of influencer for followers
    attractive = zeros(n, J)
    for j=1:n
        for i=1:J
            g = state[j]*fraction[i]
            if g<0 
                g=theta
            else
                g=g+theta
            end
            attractive[j,i]= eta * dist[j,i]*g 
        end

        r=rand()
        lambda = sum(attractive[j,:]) #total jump rate
        alpha=-log(1-r)/lambda #random number distributed due to exp(lambda)
        if dt>alpha 
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

function ABMsolve(NT = 100;  p = ABMconstruct())
    
    x, media, inf, fol, state,Net,counts  = ABMrandominitialconditions(p)
    dt, n, n_media, J, sigma, sigmahat, sigmatilde, a, b,c, frictionI, friction, eta= p

    xs = [x] #initial condition
    infs = [inf]
    meds = [media]
    xinfs = [fol * collect(1:J)]

    for k in 2:NT
        xold = x
        # opinions change due to opinions of friends, influencers and media
        force = a * attraction(xold,Net) + influence(xold,media,inf,fol,state,p)
        x = xold + dt*force + sqrt(dt*sigma)*randn(n,2); 
        # note that there are no boundary conditions, agents could escape [-2,2]x[-2,2]


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
            media[i,:] = media[i,:]  + dt/friction * (masscenter[i,:] -media[i,:]) + 1/friction * sqrt(dt*sigmatilde)*randn(2,1)
        end
        
        # individual may jump from one influencer to another
        # jumps according to rate model
        fol = changeinfluencer(state,xold,fol,inf,p)

        xs = push!(xs,copy(x))
        infs = push!(infs, copy(inf))
        meds = push!(meds, copy(media))
        xinfs = push!(xinfs, copy(fol * collect(1:J)))
    end

    return xs, xinfs, infs, meds, state, p, counts

end

function plothist(x)
    dx = 0.25 #0.05
    edges = (-2:dx:2, -2:dx:2)
    data = (x[:,1], x[:,2])
    h = fit(Histogram, data, edges)
    subp = histogram2d(data, bins=edges) #heatmap(h.weights)
    return subp
end

function sumgaussian(x, centers)
    output = 0
    for i in 1:size(centers,1)
        output = output +  gaussian(x, centers[i,:]) 
    end
    return output
end 

function kdeplot(centers, p)
    (;X, Y, domain) = p

    x_arr = domain[1,1]:dx:domain[1,2]
    y_arr = domain[2,1]:dy:domain[2,2]
    evalkde = [sumgaussian([X[i,j], Y[i,j]], centers) for i in 1:size(X,1), j in 1:size(X,2)]
    #K = kde(x, boundary = ((-2,2),(-2,2)))
    subp = heatmap(x_arr, y_arr, evalkde', c=:berlin)
    scatter!(subp, centers[:,1], centers[:,2], markercolor=:yellow,markersize=4)
    return subp
end

function ABMsolveplot(NT = 100;  p = ABMconstruct())
    xs, xinfs, infs, meds, state, p, counts = ABMsolve(NT;  p=p)
    ABMplotarray(xs[end], xinfs[end], state, p)
end

function ABMplotarray(x, xinf, state,  p; save = true)
    (; J) = p
    plot_array = Any[]  
    z_labels = ["z₋₁","z₁" ]
    y_labels = ["y₁", "y₂", "y₃", "y₄"]
    dens_labels = [  "ρ₋₁₁" "ρ₋₁₂" "ρ₋₁₃" "ρ₋₁₄";"ρ₁₁" "ρ₁₂" "ρ₁₃" "ρ₁₄"]
    states = [-1 1]
    for j in 1:J    
        for i in 1:2
            xi  = findall(x-> x==j, xinf)
            xm = findall(x-> x==states[i], state)
            choice = intersect(xi, xm)

            # make a plot and add it to the plot_array
            push!(plot_array, kdeplot(x[choice,:], p))
            #push!(plot_array, plothist(x[choice,:]))
        end
    end
    plot(plot_array..., layout=(4,2),size=(1000,1000)) |> display

    if save==true
        savefig("img/abmarray.png")
    end
end

function ABMgifarray(xs, xinfs, state, p; dN=10)
    NT=size(xs,1)
    #cl = (0, maximum(maximum(sol(t).x[1]) for t in 0:dt:tmax)) #limits colorbar
 
    abmgif = @animate for t = 1:dN:NT
        ABMplotarray(xs[t], xinfs[t], state,  p; save = false)
    end
    Plots.gif(abmgif, "img/ABMevolution.gif", fps = 30)
end
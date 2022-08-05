using LinearAlgebra
using StatsBase
using Plots
using JLD2

#need to rename all functions to be unique from PDE???
# TODO: adapt to different numbers of influencers and controlled influencers!

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

function construct()
   
    # setting model and simulation parameters
    dt=0.01  # simulation stepsize
    n = 128 # number of agents
    n_media = 2
    n_inf = 4
    sigma=0.5 # noise on individual agents 
    sigmahat=0 # noise on influencers
    sigmatilde=0 # noise on media
    a=1 # interaction strength between agents 
    b=2; # interaction strength between agents and influencers
    c=2 # interaction strength between agents and media
    frictionI = 25 # friction for influencers
    friction = 100  #friction for media
    eta = 50  #rate constant for changing influencer 
 
    p = (; dt, n, n_media, n_inf, sigma, sigmahat, sigmatilde, a, b,c, frictionI, friction, eta)

    return p
end

function random_initialconditions(p)
    (; n, n_media, n_inf) = p

    # agent opinions
    x = rand(n,2).*4 .-2 
    # media opinions
    media =[-1. 1.; -1. 1.]
    # agents following different influencers
    x_I1 = intersect(findall(x-> x>0,x[:,1]), findall(x-> x>0,x[:,2])) 
    x_I2 = intersect(findall(x-> x<=0,x[:,1]), findall(x-> x>0,x[:,2]))
    x_I3 = intersect(findall(x-> x>0,x[:,1]), findall(x-> x<=0,x[:,2]))  
    x_I4 = intersect(findall(x-> x<=0,x[:,1]), findall(x-> x<=0,x[:,2]))  

    #follower network
    fol=zeros(n, n_inf)
    fol[x_I1,1] .=1
    fol[x_I2,2] .=1
    fol[x_I3,3] .=1
    fol[x_I4,4] .=1

    #initial opinions of influencers
    inf = zeros(n_inf,2)
    inf[1,:] = sum(x[x_I1,:],dims=1)/size(x_I1,1)
    inf[2,:]= sum(x[x_I2,:],dims=1)/size(x_I2,1)
    inf[3,:] = sum(x[x_I3,:],dims=1)/size(x_I3,1)
    inf[4,:]= sum(x[x_I4,:],dims=1)/size(x_I4,1)


    # initialization of political attitude state of all agents
    state = (rand(n).>0.5)*2 .-1
    x_M1 = findall(x->x==-1, state)
    x_M2 = findall(x->x==1, state)

    #initialization of interaction network between agents 
    Net = ones(n,n)  # everyone connected to everyone including self-connections

    return x, media, x_I1, x_I2, x_I3, x_I4, inf, fol, state, x_M1, x_M2, Net 

end

function attraction(x, Net)
    n = size(Net, 1)
    force = zeros(n,2)
    for j in 1:n
        J=findall(x->x==1, Net[j,:]) # TODO: maybe simplify for complete network
        if isempty(J)
            force[j,:]=[0 0]
        else
            fi = [0 0]
            wsum=0
            for i in 1:size(J,1)
                d = x[J[i],:]-x[j,:]
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
    (; n, b, c, n_inf) = p
    force1 =zeros(size(x))
    force2 =zeros(size(x))
    for j in 1:n
        if state[j]==1
            force1[j,:] = media[:,2]-x[j,:]
        else  
            force1[j,:] = media[:,1]-x[j,:]
        end

        for k in 1:n_inf
            if fol[j,k]==1
                force2[j,:]=inf[k,:] -x[j,:]
            end
        end
    end
    force = c*force1 + b*force2
    return force
end

function changeinfluencer(state,x,fol,inf,p)
    (; dt, eta, n, n_inf) = p
    
    theta =0.1 #threshold for discrete g-function
    
    # compute happiness = fraction of followers with same state
    fraction = zeros(n_inf)
    for i=1:n_inf
        fraction[i]= sum(fol[:,i].*state)/sum(fol[:,i])
    end
    
    #ompute distance of followers to influencers
    dist = zeros(n, n_inf)
    for i=1:n_inf
        for j=1:n
            d = x[j,:]-inf[i,:]
            dist[j,i]= exp(-(d[1]^2+d[2]^2))
        end
    end
    
    # compute attractiveness of influencer for followers
    attractive = zeros(n, n_inf)
    for j=1:n
        for i=1:n_inf
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
            fol[j,:]=zeros(n_inf)      
            fol[j,k]=1
        end
    end

    return fol
end

function ABMsolve(NT = 100;  p = construct())
    
    x, media, x_I1, x_I2, x_I3, x_I4, inf, fol, state, x_M1, x_M2, Net  = random_initialconditions(p)
    dt, n, n_media, n_inf, sigma, sigmahat, sigmatilde, a, b,c, frictionI, friction, eta= p

    xs = [x] #initial condition
    infs = [inf]
    meds = [media]
    xinfs = [fol * collect(1:n_inf)]

    for k in 2:NT
        xold = x
        # opinions change due to opinions of friends, influencers and media
        force = a * attraction(xold,Net) + influence(xold,media,inf,fol,state,p)
        x = xold + dt*force + sqrt(dt*sigma)*randn(n,2); 
        # note that there are no boundary conditions, agents could escape [-2,2]x[-2,2]


        # influencer opinions adapt slowly to opinions of followers with friction
        # depending on number of followers
        masscenter=zeros(n_inf,2)
        for i in 1:n_inf
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
        print(media)
        xs = push!(xs,copy(x))
        infs = push!(infs, copy(inf))
        meds = push!(meds, copy(media))
        xinfs = push!(xinfs, copy(fol * collect(1:n_inf)))
    end

    return xs, xinfs, infs, meds, state, p

end

function plothist(x, choice)
    dx = 0.2 #0.05
    edges = (-2:dx:2, -2:dx:2)
    data = (x[choice,1], x[choice,2])
    h = fit(Histogram, data, edges)
    subp = histogram2d(data, bins=edges) #heatmap(h.weights)
    return subp
end

function plotall(x, xinf, state,  n_inf)
    plot_array = Any[]  
    z_labels = ["z₋₁","z₁" ]
    y_labels = ["y₁", "y₂", "y₃", "y₄"]
    dens_labels = [  "ρ₋₁₁" "ρ₋₁₂" "ρ₋₁₃" "ρ₋₁₄";"ρ₁₁" "ρ₁₂" "ρ₁₃" "ρ₁₄"]
    states = [-1 1]
    for j in 1:n_inf    
        for i in 1:2
            xi  = findall(x-> x==j, xinf)
            xm = findall(x-> x==states[i], state)
            choice = intersect(xi, xm)

            # make a plot and add it to the plot_array
            push!(plot_array, plothist(x, choice))
        end
    end
    plot(plot_array..., layout=(4,2),size=(1000,1000)) |> display

    #p2 = plot_solution(sol(tmax).x[1][:,:,2], sol(tmax).x[2][:,2], x_arr, y_arr; title=string("ρ₋₁(",string(tmax),")"), label="z₋₁")
    #plot(p1, p2, layout=[4 4], size=(1000,4*400)) 
    savefig("img/finaltime_abm.png")
end
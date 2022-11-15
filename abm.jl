using LinearAlgebra
using StatsBase
using Plots
using JLD2
pyplot()

#need to rename all functions to be unique from PDE???
# TODO: adapt to different numbers of influencers and controlled influencers!
# maybe change function g? 
# plot: point cloud of agents and Kernel density estimation

# https://nanx.me/ggsci/reference/pal_locuszoom.html
#plotting parameters
colors_followers = ["#5CB85CFF" "#EEA236FF" "#9632B8FF" "#B8B8B8FF"] 
colors_leaders = ["#D43F3AFF"  "#357EBDFF"]
markers_readers = [:circle  :utriangle :star5  :xcross]#  :diamond :cross ]
size_leaders = 5
size_individuals = 4
cmap = :berlin  #:tempo :bilbao :grayC

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
        b=0
        eta=0
        q=(;q..., J,b,eta)
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
 
function sumgaussian(x, centers;sigma=0.1)
    output = 0
    for i in 1:size(centers,1)
        output = output +  gaussian(x, centers[i,:],sigma=sigma) 
    end
    return output
end 

function kdeplot(centers, inf, media, state, xinf, (p,q); title = "",labelx1 = "", labelx2="",labely ="", labelz ="", scenario="4inf",clim=(-Inf, Inf),color_agents=false,sigma=0.1,ylabel="")
    (;X, Y, domain, dx, dy) = p
    (;J) = q
    n= q.n

    x_arr = domain[1,1]:dx:domain[1,2]
    y_arr = domain[2,1]:dy:domain[2,2]
    evalkde = [(1/n)*sumgaussian([X[i,j], Y[i,j]], centers,sigma=sigma) for i in 1:size(X,1), j in 1:size(X,2)]

    subp = heatmap(x_arr, y_arr, evalkde', c=cmap, title = title,alpha=0.5, clims=clim,ylabel=ylabel)
    if color_agents==false
        scatter!(subp, centers[:,1], centers[:,2], markercolor=:white,markersize=size_individuals, markerstrokewidth=0.25)
    else

        states = [-1 1]

        for j in 1:J
            for i in 1:2
                xi  = findall(x-> x==j, xinf)
                xm = findall(x-> x==states[i], state)
                choice = intersect(xi, xm)
                scatter!(subp, centers[choice,1], centers[choice,2], markercolor=colors_followers[j],markershape=markers_readers[i], markersize=size_individuals, markerstrokewidth=0.25, lab = labelx2)
            end
        end


    end
    if scenario=="4inf"
        scatter!(subp, inf[1,:], inf[2,:], markercolor=colors_leaders[1],markersize=size_leaders, lab=labely)
    end
    scatter!(subp, media[1,:], media[2,:], markercolor=colors_leaders[2],markersize=size_leaders, lab=labelz)
    return subp
end

function ABMsolveplot(NT = 100;  p = ABMconstruct(), q=parameters(), scenario="4inf")
    @time xs, xinfs, infs, meds, state, (p,q), counts = ABMsolve(NT;  p=p, q=q, scenario=scenario)
    (;dt) = p
    ABMplotarray(xs[end], xinfs[end],state,  infs[end], meds[end],  (p,q), dt*(NT-1), scenario=scenario)
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
            push!(plot_array, kdeplot(x[choice,:], inf[j,:], media[i,:], state[choice], xinf[choice,:], (p,q); title = title,labely = y_labels[j], labelz = z_labels[i], scenario=scenario))
        end
    end
    plot(plot_array..., layout=(J,2),size=(1000,J*250)) #|> display

    if save==true
        savefig(string("img/abm_array_",scenario,".png"))
    end
end


function ABMplotsnapshots(xs, xinfs, infs, meds, state, (p,q), ts; save = true, scenario="4inf",sigma=0.1)
    (; J) = q
    (;dt) = p
    n_snapshots = length(ts)
    plot_array = Any[]  
    for s in 1:n_snapshots
        t=ts[s]
        x=xs[t]
        inf=infs[t]
        media=meds[t]
        xinf = xinfs[t]
        ylabel =string("t = ", string(round((t-1)*dt, digits=2)))

        subp= kdeplot(x, inf', media', state, xinf, (p,q), scenario=scenario, clim=(0,0.5),color_agents=true,sigma=sigma,ylabel=ylabel)

        push!(plot_array, subp)
    end    
    gridp=plot(plot_array..., layout=(n_snapshots,1),size=(95*5,n_snapshots*50*5),link=:all)#, left_margin=10mm )
 

    for k=1:n_snapshots-1
        plot!(gridp[k],xformatter=_->"")
    end

    if save==true
        savefig(string("img/abm_snapshots_",scenario,".png"))
    end
end


function ABMplotfollowernumbers(xinfs,(p,q),scenario="4inf")
    N = size(xinfs,1) #number of timesteps
    (; J) = q
    (;dt) = p
    markers = [ "O"  "▽"  "△"  "☆"] #https://docs.julialang.org/en/v1/manual/unicode-input/
    numbers = zeros(N,J)
    for j in 1:J
        for n in 1:N
            numbers[n,j] = size(findall( x -> x == j, xinfs[n]),1)
        end
        if j==1
            plot(dt*collect(1:N), numbers[:,j], label = markers[j],legend=:bottomleft)
        else
            plot!(dt*collect(1:N), numbers[:,j], label = markers[j])
        end
    end
    savefig(string("img/abm_follower_",scenario,".png"))

end


"""
scenario="4inf"
#https://docs.juliaplots.org/latest/gallery/pyplot/generated/pyplot-ref58/ 

N = size(xinfs,1) #number of timesteps
(; J) = q
(;dt) = p
states = [-1 1]
markers = [ "O"  "▽"  "△"  "☆"] #https://docs.julialang.org/en/v1/manual/unicode-input/
colors=[:blue :black]
numbers = zeros(N,2,J)

for j in 1:J
    for i=1:2
        xm = findall(x-> x==states[i], state)
        for n in 1:N
            xinf = xinfs[n]
            xi  = findall(x-> x==j, xinf)
            
            choice = intersect(xi, xm)
            numbers[n,i,j] = size(choice,1)

        end
    end
end


areaplot(dt*collect(1:N), reshape(numbers,(N, J*2)), seriescolor = [:red :red  :green :green :blue :blue :grey :grey], fillalpha = [1 0.5 1 0.5 1 0.5 1 0.5])

savefig(string("img/abm_follower_",scenario,".png"))
"""


function ABMplotsingle(x, inf,xinf, media, state, (p,q), t; save = true, scenario="4inf")
    (; J) = q
    title =string("t = ", string(round(t, digits=2)))
    subp= kdeplot(x, inf', media', state, xinf, (p,q), scenario=scenario, title = title,color_agents=true) 

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
        ABMplotarray(xs[t], xinfs[t], state, infs[t], meds[t],  (p,q), (t-1)*dt; save = false, scenario=scenario)
    end
    Plots.gif(abmgif, string("img/ABM_array",scenario,".gif"), fps = 10)
end

function ABMgifsingle(xs, xinfs, state, infs, meds, (p,q); dN=5, scenario="4inf")
    NT=size(xs,1)
    (;dt) = p
    #cl = (0, maximum(maximum(sol(t).x[1]) for t in 0:dt:tmax)) #limits colorbar
 
    abmgif = @animate for t = 1:dN:NT
        ABMplotsingle(xs[t], infs[t], xinfs[t], meds[t],state, (p,q), (t-1)*dt; save = false, scenario=scenario)
    end
    Plots.gif(abmgif, string("img/abm_single_",scenario,".gif"), fps = 10)
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
                    us[:,:,i,j,m,k] = [(1/n)*sumgaussian([X[i,j], Y[i,j]], x[choice,:],sigma=sigma) for i in 1:size(X,1), j in 1:size(X,2)]
                    ys[:,j,m,k] = inf[j,:]
                end
                zs[:,i,m,k] = media[i,:]
            end
        end

    end

    @save string("data/abm_ensemble_",scenario,".jld2") us zs ys
    return us, zs, ys, (p,q), av_counts 
end

function ABMplotensemble(us, zs, ys, (p,q), NT; clmax = 0.5, scenario="4inf")
    (; dt) = p
    plotensemble(us, zs, ys, (p,q), dt*NT; title1 = "img/abm_ensemble", title2 = "img/abm_ensemble_influencer", clmax = clmax,scenario=scenario)
end

 
function plotensemblesnapshots(us, zs, ys, (p,q), us2, zs2, ys2, (p2,q2), tmax; scenario="4inf")
    # check if size(ys,4)==size(ys2,4) and size(us,5) == size(us2,5)
    N = size(ys,4)
    savepoints = size(us,5)
    savetimes = LinRange(0,tmax,savepoints)

    us_list = [us2, us]
    zs_list = [zs2, zs]
    ys_list = [ys2, ys]
    p_list = [p2, p]
    q_list = [q2, q]

    plot_array = Any[]  

    maxus = dropdims(maximum(dropdims(sum(us,dims=(3,4,6)),dims=(3,4,6))*(1/N),dims=(1,2)),dims=(1,2))
    maxus2 = dropdims(maximum(dropdims(sum(us2,dims=(3,4,6)),dims=(3,4,6))*(1/N),dims=(1,2)),dims=(1,2))
    for k in 1:savepoints
        clmax = maximum([maxus[k],maxus2[k]])
        for i in 1:2
            usi = us_list[i]
            zsi = zs_list[i]
            ysi = ys_list[i]
            p_i = p_list[i]
            qi = q_list[i]
    
            av_u = sum(usi,dims=6)*(1/N)
            av_z = sum(zsi,dims=4)*(1/N)
            av_y = sum(ysi,dims=4)*(1/N)
            if i==1
                legend=false
                ylabel = string("t=", string(round(savetimes[k], digits=2)))
            else
                legend=true
                ylabel=""
            end

            if k==1
                if i==1
                    title="ABM"
                else
                    title = "PDE"
                end
            else
                title=""
            end
            
            subp = plotsingle(av_u[:,:,:,:,k], av_z[:,:,k], av_y[:,:,k], (p_i,qi),savetimes[k]; save=false, scenario=scenario, labely="", labelx="", legend=legend,clim=(0,clmax),title=title,ylabel=ylabel)
            push!(plot_array, subp)

        end
    end
    gridp = plot(plot_array..., layout=(savepoints,2),size=(95*6,savepoints*35*6),link=:all)#, left_margin=10mm )

    #share axis
    for k=1:savepoints
        plot!(gridp[2k],yformatter=_->"")
 
    end



    for k=1:2*savepoints-2
        plot!(gridp[k],xformatter=_->"")
    end

    savefig(string("img/ensemblesnapshots_",scenario,".png"))

end



function runensembles(N; NT=200, tmax=2., savepoints = 5, q=parameters(),clmax=0.5,sigma=0.1) #0.025 works well with 1000 siulations
    us, zs, ys, (p,q), av_counts = solveensemble(tmax, N;savepoints=savepoints, q=q)
    plotensemble(us, zs, ys,(p,q), tmax; clmax=clmax)
    us2, zs2, ys2, (p2,q2), av_counts2 = ABMsolveensemble(NT,N; savepoints=savepoints, q=q,sigma=sigma)
    ABMplotensemble(us2, zs2, ys2, (p2,q2), NT; clmax=clmax)
    return us, zs, ys, (p,q), av_counts, us2, zs2, ys2, (p2,q2), av_counts2
end

function runensembles_noinf(N; NT=100, tmax=1.,savepoints = 5, q=parameters(),clmax=0.5,sigma=0.1)
    scenario="noinf"
    us, zs, ys, (p,q), av_counts = solveensemble(tmax, N; scenario=scenario, savepoints=savepoints, q=q)
    plotensemble(us, zs, ys, (p,q), tmax; clmax=clmax,scenario=scenario)
    us2, zs2, ys2, (p2,q2), av_counts2 = ABMsolveensemble(NT,N; scenario=scenario, savepoints=savepoints, q=q,sigma=sigma)
    ABMplotensemble(us2, zs2, ys2, (p2,q2), NT; clmax=clmax,scenario=scenario)
    return us, zs, ys, (p,q), av_counts, us2, zs2, ys2, (p2,q2), av_counts2
end
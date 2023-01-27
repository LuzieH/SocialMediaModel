using Plots
pyplot()


colors_followers = [ "#44AA99"  "#DDCC77"  "#CC6677" "#88CCEE"  ] 
colors_leaders = ["#BBBBBB" :black ]
steered_marker = :dtriangle
controlled_marker =  :xcross
markers_readers = [:circle  :utriangle]
size_leaders = 10
size_leaders_pde = 6
size_individuals = 6
cmap =reverse(cgrad(:gist_earth)) 
color_noinf = :white

### PDE

function plotsingle(u,z,y,(p,q),t; save=true, scenario="4inf", labely="influencers", labelx="media", clim=(0,Inf), legend=true, ylabel="", title = string("t=", string(round(t, digits=2))))
    (; domain, dx, dy) = p
    (;controlled_inf, controlled_med) = q
    J= q.J


    x_arr = domain[1,1]:dx:domain[1,2]
    y_arr = domain[2,1]:dy:domain[2,2]

    dens = dropdims(sum(u, dims=(3,4)), dims=(3,4))

    subp = heatmap(x_arr,y_arr, dens', title = title,ylabel=ylabel, c=cmap,clim=clim,legend=legend)

    for i in 1:2
        if controlled_med[i]==0
            scatter!(subp, [z[1,i]], [z[2,i]], markercolor=colors_leaders[2],markersize=size_leaders_pde, lab = labelx)
        elseif controlled_med[i]==1
            scatter!(subp, [z[1,i]], [z[2,i]], markercolor=colors_leaders[2],markersize=size_leaders_pde, markershape = controlled_marker,lab="")
        end
    end


    for j in 1:J
        if controlled_inf[j]==0
            scatter!(subp, [y[1,j]], [y[2,j]], markercolor=colors_leaders[1],markersize=size_leaders_pde, lab=labely)
        elseif controlled_inf[j]==1
            scatter!(subp, [y[1,j]], [y[2,j]], markercolor=colors_leaders[1],markersize=size_leaders_pde, markershape = controlled_marker,lab="")
        elseif controlled_inf[j]==3
            scatter!(subp, [y[1,j]], [y[2,j]], markercolor=colors_leaders[1],markersize=size_leaders_pde, markershape = steered_marker,lab="")
        end
    end



    if save==true
        savefig(string("img/pde_single_",scenario,".png"))
        savefig(string("img/pde_single_",scenario,".pdf"))
    end
    return subp
end


plotsnapshots(sol, (p,q), args...; kwargs...) = plotsnapshots([sol], [(p,q)], args...; kwargs...)

function plotsnapshots(sols::Vector, Ps::Vector, ts; save = true, scenario="4inf",followercount=false)
    n_snapshots = length(ts)
    n_sols =  size(sols,1)
    plot_array = Any[]  
    sum_t = 0
    counter =1
    for (sol, P) in zip(sols, Ps)
        for t in ts
            if t>= sum_t && t< sum_t + sol.t[end]  
                u,z,y = sol2uyz(sol, t-sum_t)
                ylabel =string("t = ", string(round(t, digits=2)))
                title=""
                if followercount == true && t>=5
                    (p,q) = P
                    sum_followers = sum(sol(t-sum_t).x[1][:,:,:,end])/sum(sol(t-sum_t).x[1])
                    label =string(L" \Sigma_m \, n_{m,5} =  ",string(round(sum_followers,digits = 3)))
                else
                    label=""
                end
                subp = plotsingle(u,z,y,P,t,save=false,scenario=scenario, labely="", labelx="",ylabel=ylabel,title=title)
                annotate!([(-1.3, 1.7,text(label,10))])
                push!(plot_array, subp)
            elseif counter==n_sols && t>= sum_t && t== sum_t + sol.t[end]
                u,z,y = sol2uyz(sol, t-sum_t)
                ylabel =string("t = ", string(round(t, digits=2)))
                title=""
                if followercount == true && t>=5
                    (p,q) = P
                    sum_followers = sum(sol(t-sum_t).x[1][:,:,:,end])/sum(sol(t-sum_t).x[1])
                    label =string(L" \Sigma_m\, n_{m,5} =  ",string(round(sum_followers,digits = 3)))
                else
                    label=""
                end
                subp = plotsingle(u,z,y,P,t,save=false, scenario=scenario, labely="", labelx="",ylabel=ylabel,title=title)
                annotate!([(-1.3, 1.7,text(label,10))])
                push!(plot_array, subp)
            end
        end
        sum_t+= sol.t[end]
        counter+=1
    end    
    gridp=plot(plot_array..., layout=(n_snapshots,1),size=(95*5,n_snapshots*50*5),link=:all)#, left_margin=10mm )

    for k=1:n_snapshots-1
        plot!(gridp[k],xformatter=_->"")
    end

    if save==true
        savefig(string("img/pde_snapshots_",scenario,".png"))
        savefig(string("img/pde_snapshots_",scenario,".pdf"))
    end
end

gifsingle(sol, (p,q), args...; kwargs...) = gifsingle([sol], [(p,q)], args...; kwargs...)

function gifsingle(sols::Vector, Ps::Vector, dt=0.1; scenario = "4inf",legend=false)
    T = 0
    anim = Animation()
    for (sol, P) in zip(sols, Ps)
        for t in 0:dt:sol.t[end]
            u,z,y = sol2uyz(sol, t)
            plt = plotsingle(u,z,y,P,t+T,save=false, scenario=scenario,legend=legend)
            frame(anim, plt)
        end
        T += sol.t[end]
    end
    Plots.gif(anim, string("img/pde_single_",scenario,".gif"), fps = 10)
end

 

### ABM

function ABMplotsingle(centers, inf, media, state, xinf, (p,q); title = "", labelx2="",labely ="", labelz ="", clim=(-Inf, Inf),color_agents=false,sigma=0.1,ylabel="", kde =true)
    (;X, Y, domain, dx, dy) = p
    (;J) = q
    n= q.n

    x_arr = domain[1,1]:dx:domain[1,2]
    y_arr = domain[2,1]:dy:domain[2,2]
    if kde==true
        evalkde = [(1/n)*sumgaussian([X[i,j], Y[i,j]], centers,sigma=sigma) for i in 1:size(X,1), j in 1:size(X,2)]
        subp = heatmap(x_arr, y_arr, evalkde', c=cmap, title = title,alpha=0.9, clims=clim,ylabel=ylabel)
    else
        subp=plot(ylabel=ylabel)
    end
    if color_agents==false
        scatter!(subp, centers[:,1], centers[:,2], markercolor=color_noinf,markersize=size_individuals, markerstrokewidth=0.5)
    else
        states = [-1 1]
        for j in 1:J
            for i in 1:2
                xi  = findall(x-> x==j, xinf)
                xm = findall(x-> x==states[i], state)
                choice = intersect(xi, xm)
                if J>1
                    scatter!(subp, centers[choice,1], centers[choice,2], markercolor=colors_followers[j],markershape=markers_readers[i], markersize=size_individuals, markerstrokewidth=0.5, lab = labelx2)
                else
                    scatter!(subp, centers[choice,1], centers[choice,2], markercolor=color_noinf,markershape=markers_readers[i], markersize=size_individuals, markerstrokewidth=0.5, lab = labelx2)
                end
            end
        end
    end
    
    for j in 1:J
        scatter!(subp, [inf[1,j]], [inf[2,j]], markercolor=colors_followers[j], lab=labely,markersize=size_leaders,markerstrokewidth=1.5)#*0.5, markerstrokecolor = colors_followers[j])
    end

    for i in 1:2
        scatter!(subp, [media[1,i]], [media[2,i]], markercolor=colors_leaders[2],markersize=size_leaders, markershape = markers_readers[i], lab=labelz)
    end
    return subp
end

function ABMplotsnapshots(xs, xinfs, infs, meds, state, (p,q), ts; save = true, scenario="4inf",sigma=0.1,kde=false)
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
        subp= ABMplotsingle(x, inf', media', state, xinf, (p,q),  clim=(0,0.5),color_agents=true,sigma=sigma,ylabel=ylabel,kde=kde)
        push!(plot_array, subp)
    end    
    gridp=plot(plot_array..., layout=(n_snapshots,1),size=(95*5,n_snapshots*50*5),link=:all)
 

    for k=1:n_snapshots-1
        plot!(gridp[k],xformatter=_->"")
    end

    if save==true
        savefig(string("img/abm_snapshots_",scenario,".png"))
        savefig(string("img/abm_snapshots_",scenario,".pdf"))
    end
end


function plotfollowernumbers(xinfs,state,(p,q);scenario="4inf")
    N = size(xinfs,1) #number of timesteps
    (; J,n) = q
    (;dt) = p
    states = [-1 1]
    numbers = zeros(N,2,J)
    scolors = Any[]
    alphas = Any[]
    labels = Any[]
    for j in 1:J
        for i=1:2
            xm = findall(x-> x==states[i], state)
            for n in 1:N
                xinf = xinfs[n]
                xi  = findall(x-> x==j, xinf)
                
                choice = intersect(xi, xm)
                numbers[n,i,j] = size(choice,1)
            end
            push!(alphas, i*0.5)
            push!(scolors, colors_followers[j])
        end
    end
    props = (1/n)*reshape(numbers,(N, J*2))

    subp = areaplot(dt*collect(1:N), props, seriescolor = permutedims(scolors), fillalpha = permutedims(alphas),title="Proportion of followers",labels = permutedims(labels),size=(95*5,60*5),xlabel="t",legend=false)
    for j in 1:J
        for i in 1:2
            scatter!(subp, [0.05], [sum(props[5,1:2*(j-1) + i])-0.5*props[5,2*(j-1) + i]], markercolor=colors_followers[j],markershape=markers_readers[i],markersize=1.5*size_individuals,legend=false)
        end
    end
    savefig(string("img/abm_follower_",scenario,".png"))
    savefig(string("img/abm_follower_",scenario,".pdf"))
end


function ABMgifsingle(xs, xinfs, state, infs, meds, (p,q); dN=5, scenario="4inf")
    NT=size(xs,1)
    (;dt) = p

    abmgif = @animate for t = 1:dN:NT
        subp = ABMplotsingle(xs[t], infs[t]', meds[t]',state,xinfs[t], (p,q), title = string("t = ", string(round((t-1)*dt, digits=2))),color_agents=true)
        plot(subp)
    end
    Plots.gif(abmgif, string("img/abm_single_",scenario,".gif"), fps = 10)
end

### ENSEMBLE

function plotensemblesnapshots(us, zs, ys, (p,q), us2, zs2, ys2, (p2,q2), tmax; scenario="4inf")
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
                    title = "ABM"
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
    gridp = plot(plot_array..., layout=(savepoints,2),size=(95*6,savepoints*27*6),link=:all)#, left_margin=10mm )
    #share axis
    for k=1:savepoints
        plot!(gridp[2k],yformatter=_->"")
 
    end
    for k=1:2*savepoints-2
        plot!(gridp[k],xformatter=_->"")
    end
    savefig(string("img/ensemblesnapshots_",scenario,".png"))
    savefig(string("img/ensemblesnapshots_",scenario,".pdf"))
end

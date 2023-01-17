using Plots
pyplot()


### parameters
# https://nanx.me/ggsci/reference/pal_locuszoom.html
#plotting parameters
#colors_followers = ["#5CB85CFF" "#EEA236FF" "#9632B8FF" "#B8B8B8FF"] 
#colors_leaders = [ "#357EBDFF" "#D43F3AFF" ]
#https://personal.sron.nl/~pault/
#colors_followers = [ "#EE6677"  "#228833" "#CCBB44" "#66CCEE"  ] # "#4477AA" 
#colors_leaders = ["#BBBBBB" "#4477AA" ]
# https://nanx.me/ggsci/reference/pal_locuszoom.html
colors_followers = [ "#44AA99"  "#DDCC77"  "#CC6677" "#88CCEE"  ] # [cgrad(:gist_earth, 20, categorical = true)[i] for i in [6 15 10 18]] #
colors_leaders = ["#BBBBBB" :black ]
steered_marker = :dtriangle
controlled_marker =  :xcross #:star6
markers_readers = [:circle  :utriangle]#  :diamond :cross ]
size_leaders = 10
size_leaders_pde = 6
size_individuals = 6
cmap =reverse(cgrad(:gist_earth)) #reverse(cgrad(:lapaz)) :batlow  
color_noinf = :white

#####PDE

function plot_solution(rho, z, y, x_arr, y_arr; title="", labelz="", labely="", clim=(-Inf, Inf), scenario="4inf")
    subp = heatmap(x_arr,y_arr, rho', title = title, c=cmap, clims=clim)
    scatter!(subp, z[1,:], z[2,:], markercolor=colors_leaders[2],markersize=size_leaders_pde, lab=labelz)
    if scenario!="noinf"
        scatter!(subp, y[1,:], y[2,:], markercolor=colors_leaders[1],markersize=size_leaders_pde, lab=labely)
    end
    return subp
end

plotsnapshots(sol, (p,q), args...; kwargs...) = plotsnapshots([sol], [(p,q)], args...; kwargs...)

function gifarray(sols::Vector, Ps::Vector, dt=0.1; scenario = "4inf")
    T = 0
    anim = Animation()
    for (sol, P) in zip(sols, Ps)
        for t in 0:dt:sol.t[end]
            u,z,y = sol2uyz(sol, t)
            plt = plotarray(u,z,y, P, t+T; save=false)
            frame(anim, plt)
        end
        T += sol.t[end]
    end
    Plots.gif(anim, string("img/pde_array_",scenario,".gif"), fps = 10)
end


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

function plotarray(u,z,y, (p,q), t; save=true, clmax = maximum(u), scenario="4inf")
    (; domain, dx, dy, dV) = p
    J= q.J

    #u,z,y = sol2uyz(sol, t)

    x_arr = domain[1,1]:dx:domain[1,2]
    y_arr = domain[2,1]:dy:domain[2,2]

    array = Any[]
    z_labels = ["z₋₁","z₁" ]
    y_labels = ["y₁", "y₂", "y₃", "y₄", "y₅"]
    dens_labels = [  "ρ₋₁₁" "ρ₋₁₂" "ρ₋₁₃" "ρ₋₁₄" "ρ₋₁₅";"ρ₁₁" "ρ₁₂" "ρ₁₃" "ρ₁₄" "ρ₁₅"]
    for j in 1:J
        for i in 1:2
            # make a plot and add it to the array
            cl = (0, clmax) #limits colorbar
            title = string(dens_labels[i,j],"(", string(round(t, digits=2)), "), prop = ", string(round(sum(u[:,:,i,j]*dV), digits = 3)))
            push!(array, plot_solution(u[:,:,i,j], z[:,i], y[:,j], x_arr, y_arr; title = title,labely = y_labels[j], labelz = z_labels[i], clim=cl, scenario=scenario))
        end
    end
    plt = plot(array..., layout=(J,2),size=(1000,min(J*250,1000)))
    plt #|> display

    if save==true
        savefig(string("img/pde_array_",scenario,".png"))
        savefig(string("img/pde_array_",scenario,".pdf"))
    end

    return plt
end

#=
function gifarray

    tmax=sol.t[end]
    #cl = (0, maximum(maximum(sol(t).x[1]) for t in 0:dt:tmax)) #limits colorbar

    pdegif = @animate for t = 0:dt:tmax
        u,z,y = sol2uyz(sol, t)
        plotarray(u,z,y, P, t; save=false)
    end
    Plots.gif(pdegif, string("img/pde_array_",scenario,".gif"), fps = 10)
end
=#

gifarray(sol, (p,q), args...; kwargs...) = gifarray([sol], [(p,q)], args...; kwargs...)

function gifarray(sols::Vector, Ps::Vector, dt=0.1; scenario = "4inf")
    T = 0
    anim = Animation()
    for (sol, P) in zip(sols, Ps)
        for t in 0:dt:sol.t[end]
            u,z,y = sol2uyz(sol, t)
            plt = plotarray(u,z,y, P, t+T; save=false)
            frame(anim, plt)
        end
        T += sol.t[end]
    end
    Plots.gif(anim, string("img/pde_array_",scenario,".gif"), fps = 10)
end

function plotsingle(u,z,y,(p,q),t; save=true, scenario="4inf", labely="influencers", labelx="media", clim=(0,Inf), legend=true, ylabel="", title = string("t=", string(round(t, digits=2))))
    (; domain, dx, dy) = p
    (;controlled, controlled_med) = q
    J= q.J
    #u,z,y = sol2uyz(sol, t)

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

    if scenario!="noinf"
        for j in 1:J
            if controlled[j]==0
                scatter!(subp, [y[1,j]], [y[2,j]], markercolor=colors_leaders[1],markersize=size_leaders_pde, lab=labely)
            elseif controlled[j]==1
                scatter!(subp, [y[1,j]], [y[2,j]], markercolor=colors_leaders[1],markersize=size_leaders_pde, markershape = controlled_marker,lab="")
            elseif controlled[j]==3
                scatter!(subp, [y[1,j]], [y[2,j]], markercolor=colors_leaders[1],markersize=size_leaders_pde, markershape = steered_marker,lab="")
            end
        end
    end

    #subp # |> display

    if save==true
        savefig(string("img/pde_single_",scenario,".png"))
        savefig(string("img/pde_single_",scenario,".pdf"))
    end
    return subp
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

function test_f()
    p = PDEconstruct()
    q = parameters()
    uzy0 = initialconditions((p,q))
    duzy = copy(uzy0)
    @time f(duzy, uzy0, (p,q), 0)
    return duzy
end

 
function plotensemble(us, zs, ys, (p,q), tmax; title1 = "img/pde_ensemble", title2 = "img/pde_ensemble_influencer_", clmax = 0.5, scenario="4inf")
    N = size(ys,4)
    savepoints = size(us,5)
    av_u = sum(us,dims=6)*(1/N)
    av_z = sum(zs,dims=4 )*(1/N)
    av_y = sum(ys,dims=4 )*(1/N)
    savetimes = LinRange(0,tmax,savepoints)

    for k in 1:savepoints
        plotarray(av_u[:,:,:,:,k], av_z[:,:,k], av_y[:,:,k], (p,q), savetimes[k]; save=false, clmax = clmax, scenario = scenario)
        savefig(string(title1,string(k),scenario,".png"))
        savefig(string(title1,string(k),scenario,".pdf"))
        if scenario=="4inf"
            (;X, Y, domain, dx, dy) = p
            J= q.J
            x_arr = domain[1,1]:dx:domain[1,2]
            y_arr = domain[2,1]:dy:domain[2,2]

            yall = reshape(ys[:,:,k,:], (2,N*J))'
            evalkde = [sumgaussian([X[i,j], Y[i,j]], yall) for i in 1:size(X,1), j in 1:size(X,2)]
            heatmap(x_arr, y_arr, evalkde', c=cmap, title=string("Distribution of influencers at time ", string(round(savetimes[k], digits=2)))) #|> display
            savefig(string(title2,string(k),".png"))
            savefig(string(title2,string(k),".pdf"))
        end
    end


end


### ABM

function kdeplot(centers, inf, media, state, xinf, (p,q); title = "",labelx1 = "", labelx2="",labely ="", labelz ="", scenario="4inf",clim=(-Inf, Inf),color_agents=false,sigma=0.1,ylabel="", kde =true)
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
    if scenario=="4inf"
        for j in 1:J
            #scatter!(subp, [inf[1,j]], [inf[2,j]], markercolor="#828282",markersize=size_leaders,markerstrokewidth=0.25, lab=labely)
            scatter!(subp, [inf[1,j]], [inf[2,j]], markercolor=colors_followers[j], lab=labely,markersize=size_leaders,markerstrokewidth=1.5)#*0.5, markerstrokecolor = colors_followers[j])
        end

    end

    for i in 1:2
        scatter!(subp, [media[1,i]], [media[2,i]], markercolor=colors_leaders[2],markersize=size_leaders, markershape = markers_readers[i], lab=labelz)
    end
    return subp
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
        savefig(string("img/abm_array_",scenario,".pdf"))
    end
end


function ABMplotsnapshots(xs, xinfs, infs, meds, state, (p,q), ts; save = true, scenario="4inf",sigma=0.1,kde=false)
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

        subp= kdeplot(x, inf', media', state, xinf, (p,q), scenario=scenario, clim=(0,0.5),color_agents=true,sigma=sigma,ylabel=ylabel,kde=kde)

        push!(plot_array, subp)
    end    
    gridp=plot(plot_array..., layout=(n_snapshots,1),size=(95*5,n_snapshots*50*5),link=:all)#, left_margin=10mm )
 

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
            # if j==1
            #     if i==1
            #         push!(labels, "O") #https://docs.julialang.org/en/v1/manual/unicode-input/
            #     else
            #         push!(labels,  "△")
            #     end
            # else
            #     push!(labels, "")
            # end
        end
    end
    props = (1/n)*reshape(numbers,(N, J*2))
    #https://docs.juliaplots.org/latest/gallery/pyplot/generated/pyplot-ref58/ 
    subp = areaplot(dt*collect(1:N), props, seriescolor = permutedims(scolors), fillalpha = permutedims(alphas),title="Proportion of followers",labels = permutedims(labels),size=(95*5,60*5),xlabel="t",legend=false)
    for j in 1:J
        for i in 1:2
            scatter!(subp, [0.05], [sum(props[5,1:2*(j-1) + i])-0.5*props[5,2*(j-1) + i]], markercolor=colors_followers[j],markershape=markers_readers[i],markersize=1.5*size_individuals,legend=false)
        end
    end
    savefig(string("img/abm_follower_",scenario,".png"))
    savefig(string("img/abm_follower_",scenario,".pdf"))
end


function plotorder(orderparameters, ord_infs, ord_meds, (p,q))
    N = size(ord_infs,1) #number of timesteps
    (; J,n) = q
    (;dt) = p    
    plot(dt*collect(0:N-1), ord_infs,label="wrt to influencers",size=(95*5,60*5),legend=:bottomright,title="Orderparameter",xlabel="t")
    plot!(dt*collect(0:N-1), ord_meds,label="wrt to media")
    savefig(string("img/abm_order_",scenario,".png"))
    savefig(string("img/abm_order_",scenario,".pdf"))

end

function ABMplotsingle(x, inf,xinf, media, state, (p,q), t; save = true, scenario="4inf")
    (; J) = q
    title =string("t = ", string(round(t, digits=2)))
    subp= kdeplot(x, inf', media', state, xinf, (p,q), scenario=scenario, title = title,color_agents=true) 

    plot(subp) #|> display

    if save==true
        savefig(string("img/abm_single_",scenario,".png"))
        savefig(string("img/abm_single_",scenario,".pdf"))
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


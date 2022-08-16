using DifferentialEquations, LinearAlgebra
using RecursiveArrayTools
using Plots
using JLD2

# TODO: adapt to different numbers of influencers and controlled influencers! (that follow a certain given trajectory


function generate_K(grid_points)
    N, d = size(grid_points)
    @assert d == 2
    k = zeros(N, N, 2)
    w = zeros(N, N)

    @inbounds for j in 1:N, i in 1:N
        x = grid_points[j,1] - grid_points[i,1] #first component of x' -x
        y = grid_points[j,2] - grid_points[i,2] #2nd component of x' -x
        w[i,j] = ww = exp(-sqrt(x^2 + y^2))
        k[i,j, 1] = ww * x
        k[i,j, 2] = ww * y
    end
    return k, w
end

function g(x)
    if x>0
        return 0.1+x
    else
        return 0.1
    end
end

function gamma(grid_points, rho, y, eta,dV)
    Nx, Ny, I, J = size(rho)
    rate = zeros((Nx*Ny, I, J, J))
    m = dV*sum(rho, dims=(1,2))[1,1,:,:]
    @inbounds for i=1:I #current attitude
        for j in 1:J #current influencer
            for j2 in 1:J #future influencer
                if j2!=j
                    for x in 1:Nx*Ny
                        rate[x,i,j,j2] = exp(-norm(grid_points[x,:]-y[:,j2]))*g((m[i,j2]-m[mod1(i+1,Int(I)),j2])/(m[i,j2]+m[mod1(i+1,Int(I)),j2]))
                    end
                end
            end
        end
    end
    return eta*reshape(rate,(Nx, Ny, I, J, J))
end


function agent_force(rho, K_matrix, W_matrix,  dV)
    force = zeros(size(rho)..., 2)
    @views for d in 1:2
        f = vec(force[:,:,d])
        f .= dV * K_matrix[:,:,d] * vec(rho)
    end

    norms = dV * W_matrix * vec(rho)

    return force./reshape(norms, size(rho))
end

function follower_force(z, grid_points, N_x, N_y)
    pointwise_int =   z' .- grid_points
    force = reshape(pointwise_int,N_x,N_y, 2)
    return force
end


function second_derivative((N_x,N_y), (dx, dy))
    # here matrix should have different shape for Ny not equal Nx
    M = Tridiagonal(ones(N_x-1), fill(-2., N_x), ones(N_x-1))
    # boundary conditions change in time to ensure flux balance, they are added when solving pde
    M[1,1] = -1.
    M[end,end] = -1.
    M .= (1/dx^2)*M

    return M
end

function centered_difference((N_x,N_y), (dx, dy))
    #centered first difference for force, doesnt work for different x, y grids
    C = 1/(2*dx)*Tridiagonal(-ones(N_x-1), zeros(N_x), ones(N_x-1))
    # at the boundary a one-sided difference scheme is used
    C[1,1:2] = 1/(dx)* [-1,1]
    C[end,end-1:end] = 1/(dx)* [-1,1]

    return C
end

function parameters(;J=4, b=2., eta=15.0)
    a = 2. #a=1 in paper, interaction strength between agents 
    #b = 2. # interaction strength between agents and influencers
    c = 1. # interaction strength between agents and media
    #eta = 15.0 #rate constant for changing influencer
    n = 200 #number of agents, important for random initial conditions
    #J = 4 #number of influencers
    n_media = 2
    sigma = 0.5 # noise on individual agents 
    sigmahat = 0 # noise on influencers
    sigmatilde = 0 # noise on media
    frictionI = 25 # friction for influencers
    frictionM = 100  #friction for media

    q = (; n, J, n_media, frictionM, frictionI, a, b, c, eta, sigma, sigmahat, sigmatilde)
    return q
end

function PDEconstruct()
    # Define the constants for the PDE
    dx = 0.05
    dy = dx
    dV = dx*dy
    domain = [-2 2; -2 2]
    N_x = Int((domain[1,2]-domain[1,1])/dx+1)
    N_y = N_x #so far only works if N_y = N_x
    N = N_x*N_y
    X = [x for x in domain[1,1]:dx:domain[1,2], y in domain[2,1]:dy:domain[2,2]]
    Y = [y for x in domain[1,1]:dx:domain[1,2], y in domain[2,1]:dy:domain[2,2]]
    grid_points = [vec(X) vec(Y)] 
    # matrix of K evaluated at gridpoints
    K_matrix, W_matrix  = generate_K(grid_points)

    M = second_derivative((N_x,N_y), (dx, dy))
    C = centered_difference((N_x,N_y), (dx, dy))
 
    p = (; grid_points, dx, dy, dV, X, Y, N, N_x, N_y, domain, K_matrix, W_matrix, C,  M)

    return p
end

function initialconditions(P)
    (; N_x , N_y,  J, n) = P
    rho_0 = zeros(N_x, N_y, 2, 4) 
    mid_y =Int(round(N_y/2))
    mid_x =Int(round(N_x/2))
    rho_0[1:mid_x, 1:mid_y,:,4] .= 1/(2*J*4)
    rho_0[1:mid_x, mid_y+2:end,:,2] .= 1/(2*J*4)
    rho_0[mid_x+2:end, 1:mid_y,:,3] .= 1/(2*J*4)
    rho_0[mid_x+2:end, mid_y+2:end,:,1] .= 1/(2*J*4)
    #cat(1/32*ones(N_x,N_y), 1/32*ones(N_x,N_y), dims=3)
    #rho_0 = fill(1/(2*J*4*4), N_x, N_y, 2, J)
    u0 = rho_0
    z2_0 = [1.,1.]
    z1_0 = [-1.,-1.]
    z0 = [z1_0  z2_0]
    y4_0 = [-1.,-1.]
    y2_0 = [-1.,1.]
    y3_0 = [1.,-1.]
    y1_0 = [1.,1. ]

    y0 = [y1_0  y2_0 y3_0 y4_0]
    counts = n/(2*J)*ones(J,2) # proportion of agents that follow each influencer
    return ArrayPartition(u0,z0,y0), counts
end

#gaussian that integrates to 1 and centered at center
gaussian(x, center, sigma=0.2) = 1/(2*pi*sigma^2) * exp(-1/(2*sigma^2)*norm(x-center)^2)

function random_initialconditions(P)
    
    (; grid_points, N_x , N_y,  dV, J,n) = P
    random_pos = rand(n,2).*4 .-2 #n uniform samples in domain
    rho_0 = zeros(N_x, N_y, 2, J) 
    counts = zeros(J,2)
    y0 = zeros(2,J)

    function add_agent(agenti, infi, state,  rho_0, counts, y0)
        rho_0[:,:,state,infi]+= reshape([gaussian(grid_points[j,:], random_pos[agenti,:]) for j in 1:N_x*N_y], N_x, N_y)
        counts[infi, state]+=1
        y0[:,infi]+=random_pos[agenti,:]
        return rho_0, counts, y0
    end

    for i in 1:n
        state = rand([1,2])
        if random_pos[i,2]>0
            if random_pos[i,1]>0
                rho_0, counts, y0 = add_agent(i, 1, state,  rho_0, counts, y0)
            else
                rho_0, counts, y0 = add_agent(i, 2, state,  rho_0, counts, y0)
            end
        else
            if random_pos[i,1]>0
                rho_0, counts, y0 = add_agent(i, 3, state,  rho_0, counts, y0)
            else
                rho_0, counts, y0 = add_agent(i, 4, state,  rho_0, counts, y0)
            end
        end
    end

    rho_0 = rho_0/(sum(rho_0)*dV)
    u0 = rho_0
    z2_0 = [1.,1.]
    z1_0 = [-1.,-1.]
    z0 = [z1_0  z2_0]
    y0= y0./dropdims(sum(counts, dims=2), dims=2)'
    return ArrayPartition(u0,z0,y0), counts
end

function noinf_initialconditions(P)
    
    (; grid_points, N_x , N_y,  dV, J ,n) = P
    random_pos = rand(n,2).*4 .-2 #n uniform samples in domain
    rho_0 = zeros(N_x, N_y, 2, 1) 
    counts = zeros(1,2)
    y0 = zeros(2,1)

    function add_agent(agenti, infi, state,  rho_0, counts, y0)
        rho_0[:,:,state,infi]+= reshape([gaussian(grid_points[j,:], random_pos[agenti,:]) for j in 1:N_x*N_y], N_x, N_y)
        counts[infi, state]+=1
        y0[:,infi]+=random_pos[agenti,:]
        return rho_0, counts, y0
    end

    for i in 1:n
        state = rand([1,2])

        rho_0, counts, y0 = add_agent(i, 1, state,  rho_0, counts, y0)

    end
    rho_0 = rho_0/(sum(rho_0)*dV)
    u0 = rho_0
    z2_0 = [1.,1.]
    z1_0 = [-1.,-1.]
    z0 = [z1_0  z2_0]
    return ArrayPartition(u0,z0,y0), counts
end

function constructinitial(init,P)
    if init=="random"
        uzy0, counts = random_initialconditions(P)
    elseif init=="noinf"
        uzy0, counts = noinf_initialconditions(P)
    else
        uzy0, counts = initialconditions(P)
    end
    return uzy0, counts
end


function f(duzy,uzy,P,t)
    yield()
    (; grid_points, N_x, N_y,  a, b, c, sigma, eta, K_matrix, W_matrix, dx,dy, dV, C,  M , N, J, frictionM, frictionI) = P
    D = sigma^2 * 0.5
    u, z, y2 = uzy.x
    du, dz, dy2 = duzy.x

    rhosum = sum(u, dims=(3,4))[:,:,1,1]
    rhosum_j = sum(u, dims=3)
    rhosum_i = sum(u, dims=4)
    m_i = dV*sum(u,dims = (1,2,4))
    m_j = dV*sum(u,dims=(1,2,3))
    Fagent = agent_force(rhosum, K_matrix, W_matrix, dV)
    rate_matrix = gamma(grid_points, u, y2, eta,dV)
    for i in 1:2
        for j in 1:J
            rho = @view  u[:,:,i,j]
            drho = @view du[:,:,i,j]
            zi = @view z[:,i]
            dzi = @view dz[:,i]
            yj = @view y2[:,j]
            dyj = @view dy2[:,j]

            force = c * follower_force(zi, grid_points, N_x, N_y) + a * Fagent + b * follower_force(yj, grid_points, N_x, N_y)     
            
            div =  C * (rho .* force[:,:,1]) + (rho .* force[:,:,2]) * C'
            reac = zeros(N_x, N_y)
            for j2=1:J
                if j2!= j
                    reac += -rate_matrix[:,:,i,j,j2] .* rho + rate_matrix[:,:,i,j2,j] .* u[:,:,i,j2]
                end
            end
            
            dif = D*(M*rho + rho*M')  
            #balance fluxes at boundary (part of boundady conditions)
            dif[1,:]+= -D/dx * (force[1,:,1].*rho[1,:]) 
            dif[end,:]+= D/dx * (force[end,:,1].*rho[end,:])
            dif[:,1]+= -D/dy * (force[:,1,2].*rho[:,1])
            dif[:,end]+= D/dy * (force[:,end,2].*rho[:,end])

            drho .=  dif - div + reac

            mean_rhoi = 1/m_i[i] * dV*reshape(rhosum_i[:,:,i,:],1,N)*grid_points
            dzi .= 1/(frictionM) * (mean_rhoi' - zi)

            mean_rhoj = 1/m_j[j] * dV*reshape(rhosum_j[:,:,:,j],1,N)*grid_points
            dyj .= 1/(frictionI) * (mean_rhoj' - yj)
        end
    end

end


function sol2uyz(sol, t)
    u = sol(t).x[1]
    z = sol(t).x[2]
    y = sol(t).x[3]
    return u,z,y
end


function solve(tmax=0.1; alg=nothing, init="random", p = PDEconstruct(), q= parameters())
    
    P = (; p..., q..., init)
    uzy0, counts = constructinitial(init,P)
    
    # Solve the ODE
    prob = ODEProblem(f,uzy0,(0.0,tmax),P)
    @time sol = DifferentialEquations.solve(prob, alg, progress=true,save_everystep=true,save_start=true)
    
    return sol, P, counts
end

function solveplot(tmax=0.1; alg=nothing, init="random", p = PDEconstruct(), q= parameters())
    sol,P, counts = solve(tmax; alg=alg, init=init, p=p, q=q)

    u,z,y = sol2uyz(sol, tmax)

    plotarray(u,z,y, P, tmax)
 
    return sol, P, counts
end


function solveensemble(tmax=0.1, N=10; alg=nothing, init="random", p = PDEconstruct(), q= parameters())
    P = (; p..., q...)
    (; N_x, N_y,J) = P
 
    zs = zeros(2, 2, N)
    ys = zeros(2, J, N)
    us = zeros(N_x, N_y, 2, J,  N)
    ## TODO add threading!
    av_counts = zeros(J,2)
    Threads.@threads for i=1:N
        sol, _ , counts= solve(tmax; alg=alg, init=init, p=p, q=q)
        u,z,y = sol2uyz(sol, tmax)
        av_counts = av_counts +  counts*(1/N)
        us[:,:,:,:,i] = u
        zs[:,:,i] = z
        ys[:,:,i] = y
    end

    @save "data/ensemble.jld2" us zs ys
    return us, zs, ys, P, av_counts 
end

function plotensemble(us, zs, ys, P, tmax; title1 = "img/ensemble.png", title2 = "img/ensemble_influencer.png", clmax = 0.5)
    N = size(ys,3)
    av_u = sum(us,dims=5)*(1/N)
    av_z = sum(zs,dims=3 )*(1/N)
    av_y = sum(ys,dims=3 )*(1/N)

    plotarray(av_u, av_z, av_y, P, tmax; save=false, clmax = clmax)
    savefig(title1)
    
    (;X, Y, domain, dx, dy, J) = P
    x_arr = domain[1,1]:dx:domain[1,2]
    y_arr = domain[2,1]:dy:domain[2,2]
    

    yall = reshape(ys, (2,N*J))'
    evalkde = [sumgaussian([X[i,j], Y[i,j]], yall) for i in 1:size(X,1), j in 1:size(X,2)]
    heatmap(x_arr, y_arr, evalkde', c=:berlin, title="Final distribution of influencers") |> display
    savefig(title2)
end

function psensemble(tmax=0.1, N=10; alg=nothing, init="random")
    us, zs, ys, P, av_counts  = solveensemble(tmax, N; alg=alg, init=init)
    plotensemble(us, zs, ys, P, tmax)
    return us, zs, ys, P, av_counts 
end

function plot_solution(rho, z, y, x_arr, y_arr; title="", labelz="", labely="", clim=(-Inf, Inf))
    subp = heatmap(x_arr,y_arr, rho', title = title, c=:berlin, clims=clim)
    scatter!(subp, z[1,:], z[2,:], markercolor=:yellow,markersize=6, lab=labelz)
    scatter!(subp, y[1,:], y[2,:], markercolor=:red,markersize=6, lab=labely)
    return subp
end


function plotarray(u,z,y, P, t; save=true, clmax = maximum(u))
    (; domain, dx, dy,J) = P

    #u,z,y = sol2uyz(sol, t)

    x_arr = domain[1,1]:dx:domain[1,2]
    y_arr = domain[2,1]:dy:domain[2,2]
    
    array = Any[]  
    z_labels = ["z₋₁","z₁" ]
    y_labels = ["y₁", "y₂", "y₃", "y₄"]
    dens_labels = [  "ρ₋₁₁" "ρ₋₁₂" "ρ₋₁₃" "ρ₋₁₄";"ρ₁₁" "ρ₁₂" "ρ₁₃" "ρ₁₄"]
    for j in 1:J    
        for i in 1:2
            # make a plot and add it to the array
            cl = (0, clmax) #limits colorbar
            title = string(dens_labels[i,j],"(", string(t), "), prop = ", string(round(sum(u[:,:,i,j]*dx*dy), digits = 3))) 
            push!(array, plot_solution(u[:,:,i,j], z[:,i], y[:,j], x_arr, y_arr; title = title,labely = y_labels[j], labelz = z_labels[i], clim=cl))
        end
    end
    plot(array..., layout=(J,2),size=(1000,1000)) |> display

    if save==true
        savefig("img/pdearray.png")
    end
end
 
function gifarray(sol, P, dt=0.01)

    tmax=sol.t[end]
    #cl = (0, maximum(maximum(sol(t).x[1]) for t in 0:dt:tmax)) #limits colorbar
 
    pdegif = @animate for t = 0:dt:tmax
        u,z,y = sol2uyz(sol, t)
        plotarray(u,z,y, P, t; save=false)
    end
    Plots.gif(pdegif, "img/evolution.gif", fps = 100)
end

function plotsingle(u,z,y,P,t; save=true, inf=true)
    (; domain, dx, dy, J) = P
    #u,z,y = sol2uyz(sol, t)

    x_arr = domain[1,1]:dx:domain[1,2]
    y_arr = domain[2,1]:dy:domain[2,2] 

    dens = dropdims(sum(u, dims=(3,4)), dims=(3,4))
    subp = heatmap(x_arr,y_arr, dens', title = string("t=", string(t)), c=:berlin) 
    scatter!(subp, z[1,:], z[2,:], markercolor=:yellow,markersize=4)
    if inf==true
        scatter!(subp, y[1,:], y[2,:], markercolor=:red,markersize=4)|> display
    end

    if save==true
        savefig("img/pdesingle.png")
    end
end

function gifsingle(sol,P, dt=0.01)
    tmax=sol.t[end]
    #    cl = (0, maximum(maximum(sol(t).x[1]) for t in 0:dt:tmax)) #limits colorbar

    pdegif = @animate for t = 0:dt:tmax
        u,z,y = sol2uyz(sol, t)
        plotsingle(u,z,y,P,t,save=false)
    end
    Plots.gif(pdegif, "img/evolutionsingle.gif", fps = 100)
end


function test_f()
    p = PDEconstruct()
    q = parameters()
    P = (; p..., q...)
    uzy0 = initialconditions(P)
    duzy = copy(uzy0)
    @time f(duzy, uzy0, P, 0)
    return duzy
end

function solvenoinf(tmax; alg=nothing)
    sol, P, _ = solve(tmax; alg=alg,  init="noinf", p = PDEconstruct(), q= parameters(J=1, b=0, eta=0))
    return sol, P
end    

function solveplotnoinf(tmax; alg=nothing)
    sol, P = solvenoinf(tmax; alg=alg)
    u,z,y = sol2uyz(sol, tmax)
    plotsingle(u,z,y,P,tmax; save=true)
end
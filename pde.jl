using DifferentialEquations, LinearAlgebra
using RecursiveArrayTools
using JLD2
using Distances
using LaTeXStrings



function generate_K(grid_points)
    N, d = size(grid_points)
    @assert d == 2
    k = zeros(N, N, 2)
    w = zeros(N, N)

    @inbounds for j in 1:N, i in 1:N
        x = grid_points[j,1] - grid_points[i,1] #first component of x' -x
        y = grid_points[j,2] - grid_points[i,2] #2nd component of x' -x
        w[i,j] = ww = exp(-int_decay*sqrt(x^2 + y^2))
        k[i,j, 1] = ww * x
        k[i,j, 2] = ww * y
    end

    return k, w
end

function g(x) 
    if x>0.1
        return x
    else
        return 0.1
    end
end

function gamma(grid_points, rho, y, eta,dV)
    Nx, Ny, I, J = size(rho)
    rate = zeros((Nx*Ny, I, J, J))
    m = dV*sum(rho, dims=(1,2))[1,1,:,:]
    dists = pairwise(Euclidean(), grid_points', y)

    for i in 1:I, j in 1:J, j2 in 1:J
        j == j2 && continue
        @views @. rate[:,i,j,j2] = eta * exp.(-dists[:, j2]) * g((m[i,j2]-m[mod1(i+1,Int(I)),j2])/(m[i,j2]+m[mod1(i+1,Int(I)),j2]))
    end

    return reshape(rate,(Nx, Ny, I, J, J))
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



function uniform_initialconditions((p,q))
    (; N_x , N_y,dV, domain, dx) = p
    (; J, n) = q
    rho_0 = zeros(N_x, N_y, 2, 4)
    mid_y =Int(round(N_y/2))
    mid_x =Int(round(N_x/2))
    start_x = Int(round((domain[1,2] - 2)/dx + 1))
    end_x = N_x - start_x +1
    rho_0[start_x:mid_x, start_x:mid_y,:,4] .= 1
    rho_0[start_x:mid_x, mid_y+1:end_x,:,2] .= 1
    rho_0[mid_x+1:end_x, start_x:mid_y,:,3] .= 1
    rho_0[mid_x+1:end_x, mid_y+1:end_x,:,1] .= 1

    u0 = rho_0/(sum(rho_0)*dV)
    z2_0 = [1.,1.]
    z1_0 = [-1.,-1.]
    z0 = [z1_0  z2_0]
    y4_0 = [-1.,-1.]
    y2_0 = [-1.,1.]
    y3_0 = [1.,-1.]
    y1_0 = [1.,1. ]

    y0 = [y1_0  y2_0 y3_0 y4_0]
    counts = n/(2*J)*ones(J,2) # proportion of agents that follow each influencer
    controlled = zeros(J)
    controlled_med = zeros(2)
    return ArrayPartition(u0,z0,y0), counts, controlled, controlled_med
end

#2D gaussian that integrates to 1 and centered at center
gaussian(x, center; sigma=0.1) = 1/(2*pi*sigma^2)* exp(-1/(2*sigma^2)*norm(x-center)^2)

function inf_initialconditions((p,q))
    (; N_x , N_y,dV, domain, dx, grid_points) = p
    (; J, n) = q
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
    controlled = zeros(J)
    controlled_med = zeros(2)
    return ArrayPartition(u0,z0,y0), counts, controlled,controlled_med
end

function noinf_initialconditions((p,q))
    (; N_x , N_y,dV, domain, dx) = p
    (; J, n) = q
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
    controlled = zeros(1)
    controlled_med = zeros(2)
    return ArrayPartition(u0,z0,y0), counts, controlled,controlled_med
end

function constructinitial(scenario,(p,q))
    if scenario=="4inf"
        uzy0, counts, controlled,controlled_med = inf_initialconditions((p,q))
    elseif scenario=="noinf"
        J=1
        b=0
        eta=0
        q=(;q..., J,b,eta)
        uzy0, counts, controlled,controlled_med = noinf_initialconditions((p,q))
    elseif scenario =="uniform"
        uzy0, counts, controlled,controlled_med = uniform_initialconditions((p,q))
    elseif scenario=="fixedtarget"
        uzy0, counts, controlled,controlled_med = uniform_initialconditions((p,q))
    elseif scenario=="localmax"
        uzy0, counts, controlled,controlled_med = uniform_initialconditions((p,q))
    elseif scenario=="optimalcontrol"
        uzy0, counts, controlled,controlled_med = uniform_initialconditions((p,q))
    end
    return uzy0, counts, controlled,controlled_med,q
end


function f(duzy,uzy,(p,q),t)
    yield()
    (; a, b, c, sigma, eta,  J, frictionM, frictionI, controlled, controlled_med,controlspeed,controltarget, controlspeed2,controltarget2) = q
    (; grid_points, N_x, N_y,N, K_matrix, W_matrix, dx,dy, dV, C,  M) = p
    
    D = sigma^2 * 0.5
    u, z, y2 = uzy.x
    du, dz, dy2 = duzy.x


    rhosum = sum(u, dims=(3,4))[:,:,1,1]
    rhosum_j = sum(u, dims=3)
    rhosum_i = sum(u, dims=4)
    m_i = dV*sum(u,dims = (1,2,4))
    m_j = dV*sum(u,dims=(1,2,3))
    Fagent = agent_force(rhosum, K_matrix, W_matrix, dV)
    rate_matrix = gamma(grid_points, u, y2, eta, dV)
    reac = zeros(N_x, N_y)
    dif  = similar(reac)
    dive = similar(reac)
    rhoforce = zeros(N_x, N_y, 2)
    for i in 1:2
        for j in 1:J
            rho = @view  u[:,:,i,j]
            drho = @view du[:,:,i,j]
            zi = @view z[:,i]
            dzi = @view dz[:,i]
            yj = @view y2[:,j]
            dyj = @view dy2[:,j]

            rhoforce .= a .* Fagent .+ b .* follower_force(yj, grid_points, N_x, N_y) .+ c .* follower_force(zi, grid_points, N_x, N_y)
            rhoforce .= rho .* rhoforce

            @views mul!(dive, C, rhoforce[:,:,1])
            @views mul!(dive, rhoforce[:,:,2], C', 1, 1)
            #@assert dive == C * rhoforce[:,:,1] + rhoforce[:,:,2] * C'


            reac .= 0
            for j2=1:J
                if j2!= j
                    @. @views reac += -rate_matrix[:,:,i,j,j2] .* rho + rate_matrix[:,:,i,j2,j] .* u[:,:,i,j2]
                end
            end

            a_AB_BAT!(dif, D, M, rho)
            #@assert dif == D*(M*rho + rho*M')

            #balance fluxes at boundary (part of boundady conditions)
            dif[1,:]+= -D/dx * (rhoforce[1,:,1])
            dif[end,:]+= D/dx * (rhoforce[end,:,1])
            dif[:,1]+= -D/dy * (rhoforce[:,1,2])
            dif[:,end]+= D/dy * (rhoforce[:,end,2])

            drho .=  dif .- dive .+ reac


            if controlled[j] == 0
                mean_rhoj = 1/m_j[j] * dV*reshape(rhosum_j[:,:,:,j],1,N)*grid_points
                dyj .= 1/(frictionI) * (mean_rhoj' - yj)
            
            elseif controlled[j]==2
                #local maximization of followers
                speed = 5
                sigma_localmax = 0.5
                grad_rate_x  = (1/dx)*(gamma(grid_points, u, y2 .+ [dx, 0], eta, dV) - gamma(grid_points, u, y2, eta, dV))
                grad_rate_y = (1/dy)*(gamma(grid_points, u, y2 .+ [0, dy], eta, dV) - gamma(grid_points, u, y2, eta, dV))
                
                gradC_x  = sum(grad_rate_x[:,:,:,1,j] .* (dropdims(rhosum_i,dims=4)- u[:,:,:,j]))*dV
                gradC_y  = sum(grad_rate_y[:,:,:,1,j] .* (dropdims(rhosum_i,dims=4)- u[:,:,:,j]))*dV
                dyj .= speed* [gradC_x, gradC_y] + sigma_localmax * randn(2)
                #todo check if it works, maybe add noise to escape local minima 

            elseif controlled[j] == 1 #controll movement
                if norm(controltarget' - yj) >0
                    dyj .= controlspeed* (controltarget' - yj)./ norm(controltarget' - yj)
                else
                    dyj .= 0
                end
            elseif controlled[j] == 3 #second controlled influencer
                if norm(controltarget2' - yj) >0
                    dyj .= controlspeed2* (controltarget2' - yj)./ norm(controltarget2' - yj)
                else
                    dyj .= 0
                end
            end
            if controlled_med[i] ==1 
                if norm(controltarget' - zi) >0
                    dzi .= controlspeed* (controltarget' - zi)./ norm(controltarget' - zi)
                else
                    dzi .= 0
                end
            else
                mean_rhoi = 1/m_i[i] * dV*reshape(rhosum_i[:,:,i,:],1,N)*grid_points
                dzi .= 1/(frictionM) * (mean_rhoi' - zi)
            end
        end
    end


end

# inplace Y = a*(A*B + B+A')
function a_AB_BAT!(Y, a, A, B)
    mul!(Y, A, B)
    mul!(Y, B, A', a, a)
end


function sol2uyz(sol, t)
    u = sol(t).x[1]
    z = sol(t).x[2]
    y = sol(t).x[3]
    return u,z,y
end


function solve(tmax=0.1; alg=nothing, scenario="4inf", p = PDEconstruct(), q= parameters())

    uzy0, counts, controlled,controlled_med,q = constructinitial(scenario,(p,q))
    q = (; q..., controlled,controlled_med)

    # Solve the ODE
    prob = ODEProblem(f,uzy0,(0.0,tmax),(p,q))
    @time sol = DifferentialEquations.solve(prob, alg, alg_hints = [:stiff], save_start=true)

    return sol, (p,q), counts
end

function computeorderparameters(sol, (p,q); ts = sol.t)
    ord_infs = Any[]
    ord_meds = Any[]
    N = p.N
    grid_points = p.grid_points
    for t in ts
        u,z,y = sol2uyz(sol, t)
        ord_inf =0.
        ord_med =0.

        for l in 1:q.J
            rho = vec(dropdims(sum(u[:,:,:,l],dims=3),dims = (3)))
            for j in 1:N
                xc = grid_points[j,1] - y[1,l] #first component of x' -z_l
                yc = grid_points[j,2] - y[2,l]  #2nd component of x' -z_l
                ord_inf += rho[j]*exp(-order_decay*sqrt(xc^2 + yc^2)) *p.dV
            end
        end
        for m in 1:2
            rho = vec(dropdims(sum(u[:,:,m,:],dims=3),dims = (3)))
            @inbounds for j in 1:N
                xc = grid_points[j,1] - z[1,m] #first component of x' -z_l
                yc = grid_points[j,2] - z[2,m]  #2nd component of x' -z_l
                ord_med += rho[j]*exp(-order_decay*sqrt(xc^2 + yc^2)) *p.dV
            end
        end

        push!(ord_meds,ord_med)
        push!(ord_infs, ord_inf)
    end
    return ord_infs, ord_meds
end



function solvelocalmax(tequil = 5., tmax=10.; alg=nothing, scenario="localmax", p = PDEconstruct(), q= parameters_control(), savedt=0.05, atol = 1e-6, rtol = 1e-3)
    uzy0, _, controlled,controlled_med,q = constructinitial(scenario,(p,q))
    q1 = (; q..., controlled,controlled_med)

    # Solve the ODE
    prob1 = ODEProblem(f,uzy0,(0.0,tequil),(p,q1))

    @time sol1 = DifferentialEquations.solve(prob1, alg, alg_hints = [:stiff], saveat = 0:savedt:tequil,save_start=true, abstol = atol, reltol = rtol)

    #add new influencer
    u,z,y = sol2uyz(sol1,tequil)
    q2 = merge(q1, (;J=q1.J+1,controlled = [q1.controlled..., 2]))
    (;N_x, N_y, grid_points,N, dV) = p
    J= q2.J
    u2 = zeros(N_x, N_y, 2, J)
    u2[:,:,:,1:J-1] = u
    y2 = zeros(2, J)
    y2[:,1:J-1] = y
    #TODO maybe make constant speed such that starts and ends within 5. time steps
    startlocation =  [0 0] #1/sum(u2,dims=(1,2,4))[1,1,2,1] * reshape(sum(u2, dims=4)[:,:,2,:],1,N)*grid_points
    y2[:,J] = startlocation 
    uzy0 = ArrayPartition(u2,z,y2)

    # solve ODE with added influencer
    prob2 = ODEProblem(f,uzy0,(0.0,tmax-tequil),(p,q2))
    @time sol2 = DifferentialEquations.solve(prob2, alg, alg_hints = [:stiff],  saveat = 0:savedt:(tmax-tequil),save_start=true, abstol = atol, reltol = rtol)
    followersum = sum([sum(sol2(t).x[1][:,:,:,J])/sum(sol2(t).x[1]) for t in sol2.t])*savedt 
    return [sol1, sol2], [(p,q1), (p,q2)],  followersum
end



function plotfollowersum(followersum, Xi, Yi)
    Ngrid = size(followersum,1)
    #average out second targets
    follower_av1 = dropdims(sum(followersum,dims=(3,4)) *1/(Ngrid^2), dims = (3,4))
    #average out first targets
    follower_av2 = dropdims(sum(followersum, dims=(1,2))*1/(Ngrid^2), dims = (1,2))
    heatmap(Xi[:,1],Xi[:,1],follower_av1,c=cmap, title="Objective function, second target averaged out",xlabel="x-position of target",ylabel="y-position of target")
    savefig("img/obj_secondaveragedout.png")
    savefig("img/obj_secondaveragedout.pdf")
    heatmap(Xi[:,1],Xi[:,1],follower_av2,c=cmap,title="Objective function, first target averaged out",xlabel="x-position of target",ylabel="y-position of target")
    savefig("img/obj_firstaveragedout.png")
    savefig("img/obj_firstaveragedout.pdf")
    #TODO plot paths with obj function color on opinion domain space
end

function solveplot(tmax=0.1; alg=nothing, scenario="4inf", p = PDEconstruct(), q= parameters())
    sol,(p,q), counts = solve(tmax; alg=alg, scenario=scenario, p=p, q=q)

    u,z,y = sol2uyz(sol, tmax)

    plt = plotarray(u,z,y, (p,q), tmax)
    plt |> display

    return sol, (p,q), counts
end


function solveensemble(tmax=0.1, N=10; savepoints = 4, alg=nothing, scenario="4inf", p = PDEconstruct(), q= parameters())
    (; N_x, N_y) = p
    J=q.J

    zs = zeros(2, 2, savepoints, N)
    ys = zeros(2, J, savepoints, N)
    us = zeros(N_x, N_y, 2, J, savepoints,  N)
    savetimes = LinRange(0, tmax, savepoints)
    av_counts = zeros(J,2)
    ts = LinRange(0,tmax,100)
    ord_infs_mean = zeros(size(ts))
    ord_meds_mean =  zeros(size(ts))
    ord_infs_sqrmean = zeros(size(ts))
    ord_meds_sqrmean =  zeros(size(ts))
    #Threads.@threads 
    for i=1:N
        sol, _ , counts= solve(tmax; alg=alg, scenario=scenario, p=p, q=q)
        ord_infs, ord_meds = computeorderparameters(sol, (p,q); ts = ts)

        ord_infs_mean += ord_infs*(1/N)
        ord_meds_mean +=  ord_meds*(1/N)    
        ord_infs_sqrmean += ord_infs.^2*(1/N)
        ord_meds_sqrmean +=  ord_meds.^2*(1/N)      

        av_counts = av_counts +  counts*(1/N)
        for j in 1:savepoints
            u,z,y = sol2uyz(sol, savetimes[j])
            us[:,:,:,:,j,i] = u
            zs[:,:,j,i] = z
            ys[:,:,j,i] = y
        end

    end

    @save string("data/pde_ensemble_",scenario,".jld2") us zs ys ord_infs_mean ord_meds_mean ord_infs_sqrmean ord_meds_sqrmean
    return us, zs, ys, (p,q), av_counts, ord_infs_mean, ord_meds_mean, ord_infs_sqrmean, ord_meds_sqrmean
end


function psensemble(tmax=0.1, N=10; alg=nothing, scenario="4inf")
    us, zs, ys, (p,q), av_counts  = solveensemble(tmax, N; alg=alg, scenario=scenario)
    plotensemble(us, zs, ys, (p,q), tmax)
    return us, zs, ys, (p,q), av_counts
end


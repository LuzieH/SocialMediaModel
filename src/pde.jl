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
        w[i,j] = ww = exp(-sqrt(x^2 + y^2))
        k[i,j, 1] = ww * x
        k[i,j, 2] = ww * y
    end

    return k, w
end

"""recommender system function"""
function r(x) 
    if x>0.1
        return x
    else
        return 0.1
    end
end

"""rate function"""
function gamma(grid_points, rho, y, eta,dV)
    Nx, Ny, M, L = size(rho)
    rate = zeros((Nx*Ny, M, L,L))
    m = dV*sum(rho, dims=(1,2))[1,1,:,:]
    dists = pairwise(Euclidean(), grid_points', y)

    for i in 1:M, j in 1:L, j2 in 1:L
        j == j2 && continue
        @views @. rate[:,i,j,j2] = eta * exp.(-dists[:, j2]) * r((m[i,j2]-m[mod1(i+1,Int(M)),j2])/(m[i,j2]+m[mod1(i+1,Int(M)),j2]))
    end

    return reshape(rate,(Nx, Ny, M, L,L))
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


function uniforminit((p,q))
    (; N_x , N_y,dV, domain, dx) = p
    (; L, M) = q
    @assert M==2
    @assert L==4

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
    controlled_inf = zeros(L)
    controlled_med = zeros(2)
    return ArrayPartition(u0,z0,y0), controlled_inf, controlled_med
end

#2D gaussian that integrates to 1 and centered at center
gaussian(x, center; sigma=0.1) = 1/(2*pi*sigma^2)* exp(-1/(2*sigma^2)*norm(x-center)^2)

function randominit((p,q))
    (; N_x , N_y,dV, grid_points,domain) = p
    (; L, n, M) = q
    @assert M==2
    @assert L==4

    random_pos = rand(n,2).*(domain[1,2]-domain[1,1]) .+domain[1,1] # n uniform samples in domain
    rho_0 = zeros(N_x, N_y, 2, L)
    counts = zeros(L,2)
    y0 = zeros(2,L)

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
    controlled_inf = zeros(L)
    controlled_med = zeros(2)
    return ArrayPartition(u0,z0,y0),  controlled_inf,controlled_med
end

function constructinitial(init,(p,q))
    if init=="4inf"
        uzy0,  controlled_inf,controlled_med = uniforminit((p,q))
    elseif init =="uniform"
        uzy0, controlled_inf,controlled_med = uniforminit((p,q))
    elseif init == "random"
        uzy0, controlled_inf,controlled_med = randominit((p,q))
    end
    q = (; q..., controlled_inf,controlled_med)
    return uzy0, q
end


function f(duzy,uzy,(p,q),t)
    yield()
    (; a, b, c, sigma, eta,  L, frictionM, frictionI, controlled_inf, controlled_med,controlspeed1,controltarget1, controlspeed2,controltarget2) = q
    (; grid_points, N_x, N_y,N, K_matrix, W_matrix, dx,dy, dV, C,  MassM) = p
    
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
        for j in 1:L
            rho = @view  u[:,:,i,j]
            drho = @view du[:,:,i,j]
            zi = @view z[:,i]
            dzi = @view dz[:,i]
            yj = @view y2[:,j]
            dyj = @view dy2[:,j]

            # divergence term
            rhoforce .= a .* Fagent .+ b .* follower_force(yj, grid_points, N_x, N_y) .+ c .* follower_force(zi, grid_points, N_x, N_y)
            rhoforce .= rho .* rhoforce
            @views mul!(dive, C, rhoforce[:,:,1]) # dive = C* rhoforce[:,:,1] (matrix-matrix product)
            @views mul!(dive, rhoforce[:,:,2], C', 1, 1) # dive = rhoforce[:,:,2] * C' + dive

            # reaction term
            reac .= 0
            for j2=1:L
                if j2!= j
                    @. @views reac += -rate_matrix[:,:,i,j,j2] * rho + rate_matrix[:,:,i,j2,j] * u[:,:,i,j2]
                end
            end
            
            # diffusion term
            a_AB_BAT!(dif, D, MassM, rho) # inplace dif = D*(M*rho + rho*M')

            #balance fluxes at boundary (part of boundary conditions)
            dif[1,:]+= -D/dx * (rhoforce[1,:,1])
            dif[end,:]+= D/dx * (rhoforce[end,:,1])
            dif[:,1]+= -D/dy * (rhoforce[:,1,2])
            dif[:,end]+= D/dy * (rhoforce[:,end,2])

            drho .=  dif .- dive .+ reac

            # note: missing boundary conditions for influencer and media movement
            if controlled_inf[j] == 0
                mean_rhoj = 1/m_j[j] * dV*reshape(rhosum_j[:,:,:,j],1,N)*grid_points
                dyj .= 1/(frictionI) * (mean_rhoj' - yj)

            elseif controlled_inf[j] == 1 #first controlled influencer, movement defined by specified target and speed
                if norm(controltarget1' - yj) >0
                    dyj .= controlspeed1* (controltarget1' - yj)./ norm(controltarget1' - yj)
                else
                    dyj .= 0
                end
            elseif controlled_inf[j] == 3 #second controlled influencer
                if norm(controltarget2' - yj) >0
                    dyj .= controlspeed2* (controltarget2' - yj)./ norm(controltarget2' - yj)
                else
                    dyj .= 0
                end
            end
            if controlled_med[i] ==1 #controlled medium
                if norm(controltarget1' - zi) >0
                    dzi .= controlspeed1* (controltarget1' - zi)./ norm(controltarget1' - zi)
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


function PDEsolve(tmax=0.1; alg=nothing, init="4inf", p = PDEconstruct(), q= parameters())

    uzy0, q = constructinitial(init,(p,q))
    
    # Solve the ODE
    prob = ODEProblem(f,uzy0,(0.0,tmax),(p,q))
    @time sol = DifferentialEquations.solve(prob, alg, alg_hints = [:stiff], save_start=true)

    return sol, (p,q)
end


function PDEsolveplot(; tmax=2.0, ts = [0. 0.1 0.4 1.0 1.5 2.0], alg=nothing, init="4inf", p = PDEconstruct(), q= parameters(),save=true)
    sol,(p,q) = PDEsolve(tmax; alg=alg, init=init, p=p, q=q)
    PDEplotsnapshots(sol, (p,q), ts;  name=init,save=save)
 
    return sol, (p,q)
end



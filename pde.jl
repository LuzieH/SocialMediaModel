using DifferentialEquations, LinearAlgebra
using RecursiveArrayTools
using Plots

function generate_K(grid_points)
    N, d = size(grid_points)
    @assert d == 2
    k = zeros(N, N, 2)
    w = zeros(N, N)

    @inbounds for j in 1:N, i in 1:N
        x = grid_points[i,1] - grid_points[j,1]
        y = grid_points[i,2] - grid_points[j,2]
        w[i,j] = ww = exp(-sqrt(x^2 + y^2))
        k[i,j, 1] = ww * x
        k[i,j, 2] = ww * y
    end
    return k, w
end

function agent_force(rho::Array{T}, K_matrix, W_matrix,  dV) where T
    force = zeros(T,size(rho)..., 2)
    @views for d in 1:2
        f = vec(force[:,:,d])
        f .= dV * K_matrix[:,:,d] * vec(rho)
    end

    norms = dV * W_matrix * vec(rho)

    return force./reshape(norms, size(rho))
end

function media_force(z, grid_points, N_x, N_y)
    pointwise_int =   z' .- grid_points
    force = reshape(pointwise_int,N_x,N_y, 2)
    return force
end


function second_derivative((N_x,N_y), (dx, dy))
    M = Tridiagonal(ones(N_x-1), fill(-2., N_x), ones(N_x-1))
    # here matrix should have different shape for Ny not equal Nx
    # neumann boundary conditions with zero flux imply (using centered difference) that value of function outside domain
    # equals value of function inside domain (not on boundary), thus the following replacement
    M[1,2] = 2.0
    M[end,end-1] = 2.0
    M .= (1/dx^2)*M

    return M
end

function centered_difference_force((N_x,N_y), (dx, dy))
    #centered first difference for force, doesnt work for different x, y grids
    C = 1/(2*dx)*Tridiagonal(-ones(N_x-1), zeros(N_x), ones(N_x-1))
    C[1,1:2] = 1/(dx)* [-1,1]
    C[end,end-1:end] = 1/(dx)* [-1,1]

    return C
end

function centered_difference_density((N_x,N_y), (dx, dy))
    #centered first difference for density
    Cr = 1/(2*dx)*Tridiagonal(-ones(N_x-1), zeros(N_x), ones(N_x-1))
    Cr[1,2]=0.; Cr[end,end-1]=0.
    return Cr
end

function construct()
    # Define the constants for the PDE
    sigma = 1.0
    D = sigma^2 * 0.5
    a = 2.0
    c = 2.0
    Gamma_0 = 100
    m = [1/2, 1/2] #proportion of agents of type 1,2
    dx = 0.05
    dy = dx
    dV = dx*dy
    # domain = [[-2,2],[-2,2]]
    N_x = Int(4/dx+1)
    N_y = N_x #so far only works if N_y = N_x
    N = N_x*N_y
    X = [x for x in -2:dx:2, y in -2:dy:2]
    Y = [y for x in -2:dx:2, y in -2:dy:2]
    grid_points = [vec(X) vec(Y)]
    # matrix of K evaluated at gridpoints
    K_matrix, W_matrix  = generate_K(grid_points)

    M = second_derivative((N_x,N_y), (dx, dy))
    C = centered_difference_force((N_x,N_y), (dx, dy))
    Cr = centered_difference_density((N_x,N_y), (dx, dy))

    Δ = DiffEqArrayOperator(kern(D * M))

    p = (; Δ, grid_points, N_x, N_y, a, c, K_matrix, W_matrix, dV, C, Cr,  D, M , N, m, Gamma_0)
    return p

end

function initialconditions(p)
    (; N_x , N_y) = p
    #rho_0 = cat(1/32*ones(N_x,N_y), 1/32*ones(N_x,N_y), dims=3)
    rho_0 = fill(1/32, N_x, N_y, 2)
    u0 = rho_0
    z1_0 = [1.,1.]
    z2_0 = [-1.,-1.]
    z0 = [z1_0  z2_0]
    return ArrayPartition(u0,z0)
end

function f(duz,uz,p,t)
    yield()
    (; grid_points, N_x, N_y, a, c, K_matrix, W_matrix, dV, C, Cr,  D, M , N, m, Gamma_0) = p
    u, z = uz.x
    du, dz = duz.x

    rhosum = u[:,:,1] + u[:,:,2]
    Fagent = agent_force(rhosum, K_matrix, W_matrix, dV)
    for i in 1:2
        rho = @view  u[:,:,i]
        drho = @view du[:,:,i]
        zi = @view z[:,i]
        dzi = @view dz[:,i]

        force = c * media_force(zi, grid_points, N_x, N_y) + a * Fagent
        #ensure force is zero at boundary to conserve density
        force[1,:,:] .= force[end,:,:] .= force[:,1,:] .= force[:,end,:].=0
        div = rho .* (C*force[:,:,1] + force[:,:,2]*C') + (Cr*rho).*force[:,:,1]+ (rho*Cr').*force[:,:,2]
        drho .= D*(M*rho + rho*M') - div

        mean_rho = 1/m[i] * dV*reshape(rho,1,N)*grid_points
        dzi .= 1/(Gamma_0*m[i]) * (mean_rho' - zi)
    end

end

function diffusion(duz, uz, p, t)
    (; D, M) = p
    println("f1")
    u, z = uz.x
    du, dz = duz.x
    for i in 1:2
        rho = @view  u[:,:,i]
        drho = @view du[:,:,i]
        drho .= D*(M*rho + rho*M')
    end
end

function advection(duz, uz, p, t)
    println("f2")
    yield()
    (; grid_points, N_x, N_y, a, c, K_matrix, W_matrix, dV, C, Cr, N, m, Gamma_0) = p
    u, z = uz.x
    du, dz = duz.x

    rhosum = u[:,:,1] + u[:,:,2]
    Fagent = agent_force(rhosum, K_matrix, W_matrix, dV)
    for i in 1:2
        rho = @view  u[:,:,i]
        drho = @view du[:,:,i]
        zi = @view z[:,i]
        dzi = @view dz[:,i]

        force = c * media_force(zi, grid_points, N_x, N_y) + a * Fagent
        #ensure force is zero at boundary to conserve density
        force[1,:,:] .= force[end,:,:] .= force[:,1,:] .= force[:,end,:].=0
        div = rho .* (C*force[:,:,1] + force[:,:,2]*C') + (Cr*rho).*force[:,:,1]+ (rho*Cr').*force[:,:,2]
        drho .=  -div

        mean_rho = 1/m[i] * dV*reshape(rho,1,N)*grid_points
        dzi .= 1/(Gamma_0*m[i]) * (mean_rho' - zi)
    end
end

function diffusionop(duz, uz, p, t)
    (;Δ) = p
    println("f1")
    u, z = uz.x
    du, dz = duz.x
    for i in 1:2
        rho = vec(@view  u[:,:,i])
        drho = vec(@view du[:,:,i])
        drho .= Δ * rho
    end
end


function solve(tmax=0.01; alg=nothing)

    p = construct()

    uz0 = initialconditions(p)
    # Solve the ODE
    prob = ODEProblem(f,uz0,(0.0,tmax),p)

    @time sol = DifferentialEquations.solve(prob, alg, progress=true,save_everystep=true,save_start=false)

    p1 = heatmap(sol[end].x[1][:,:,1],title = "rho_1")
    p2 = heatmap(sol[end].x[1][:,:,2],title = "rho_{-1}")
    plot(p1, p2, layout=grid(2,1)) |> display
    return sol
end

function splitsolve(tmax=.01; alg=nothing)
    p = construct()

    uz0 = initialconditions(p)
    # Solve the ODE
    Δ = DiffEqArrayOperator(fullkern(p.D * p.M))
    prob = SplitODEProblem(Δ,advection,uz0,(0.0,tmax),p)

    D = DiffEqArrayOperator(kern(p.D * p.M))

    @time sol = DifferentialEquations.solve(prob, alg, progress=true,save_everystep=true,save_start=false)

    p1 = heatmap(sol[end].x[1][:,:,1],title = "rho_1")
    p2 = heatmap(sol[end].x[1][:,:,2],title = "rho_{-1}")
    plot(p1, p2, layout=grid(2,1)) |> display
    return sol
end

function test_f()
    p = construct()
    uz0 = initialconditions(p)
    duz = copy(uz0)
    @time f(duz, uz0, p, 0)
    return duz
end
#=
function rho_jacobian!(J, u, p, t)
    (;D, C, Cr, K_matrix) = p
    C  = [C, C']
    Cr = [Cr, Cr']
    F  = [F1, F2] # need to compute forces here?
    K  = @views [K_matrix[:,:,1], K_matrix[:,:,2]]
    J .= D # diffusive part
    for k in 1:2
        J .+= Cr[k] .* F[k] #
        J .+= Cr[k] * u .* K[k]
        J[diagind[J]] .+= C[k] * F[k]
        J .+= u .* C[k] * K[k]
    end
end
=#

# blow up 1d kernel to higher (2) dimensions
function kern(M)
    N = size(M, 1)
    MM = zeros(N,N,N,N)

    for i in 1:N, j in 1:N, k in 1:N
        MM[i,k,j,k] += M[i,j]
        MM[k,i,k,j] += M[i,j]
    end


    A = sparse(reshape(MM, N*N, N*N))
    #Z = spzeros(N*N, N*N)
    #Z2 = spzeros
    #[A Z; Z A]
end

function fullkern(M)
    A = kern(M)
    blockdiag(A, A, spzeros(4,4))
end
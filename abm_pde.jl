using DifferentialEquations, LinearAlgebra
using RecursiveArrayTools
using Plots

function generate_K(N, grid_points)  
    W(s,x) = exp(-norm(s-x))
    K(x,s) = W(s,x) * (s-x)   
    k = zeros(N, N, 2)
    w = zeros(N, N)
    for i in 1:N, j in 1:N
        k[i,j,:] = K(grid_points[i,:], grid_points[j,:])
        w[i,j] = W(grid_points[i,:], grid_points[j,:])
    end
    return k, w
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
    dx = 0.1
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
    K_matrix, W_matrix  = generate_K(N, grid_points)
    
    M = second_derivative((N_x,N_y), (dx, dy))
    C = centered_difference_force((N_x,N_y), (dx, dy))
    Cr = centered_difference_density((N_x,N_y), (dx, dy))
 
    p = (; grid_points, N_x, N_y, a, c, K_matrix, W_matrix, dV, C, Cr,  D, M , N, m, Gamma_0)
    return p, X, Y
end

function initialconditions(N_x = 41, N_y = 41)
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
        div = rho .* (C*force[:,:,1] + force[:,:,2]*C') + (Cr*rho).*force[:,:,1]+ (rho*Cr').*force[:,:,2]
        drho .= D*(M*rho + rho*M') + div

        mean_rho = 1/m[i] * reshape(rho,1,N)*grid_points
        dzi .= 1/(Gamma_0*m[i]) * (mean_rho' - zi)
    end
end


function solve(tmax=0.01; alg=nothing)
    uz0 = initialconditions()
    p, X, Y = construct()
    
    # Solve the ODE
    prob = ODEProblem(f,uz0,(0.0,tmax),p)
    
    @time sol = DifferentialEquations.solve(prob, alg, progress=true,save_everystep=true,save_start=false)

    p1 = heatmap(sol[end].x[1][:,:,1],title = "rho_1")
    p2 = heatmap(sol[end].x[1][:,:,2],title = "rho_{-1}")
    plot(p1,p2,layout=grid(2,1)) |> display
    return sol
end

function test_f()
    uz0 = initialconditions()
    duz = copy(uz0)
    p, = construct()
    @time f(duz, uz0, p, 0)
    return duz
end    

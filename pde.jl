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
    C = centered_difference_force((N_x,N_y), (dx, dy))
    Cr = centered_difference_density((N_x,N_y), (dx, dy))
 
    p = (; grid_points, N_x, N_y, domain,  a, c, K_matrix, W_matrix, dx,dy, dV, C, Cr,  D, M , N, m, Gamma_0)

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
    (; grid_points, N_x, N_y, domain, a, c, K_matrix, W_matrix, dx,dy, dV, C, Cr,  D, M , N, m, Gamma_0) = p
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
        #TODO: smoothen the force to zero at boundary!
        force[1,:,:] .= force[end,:,:] .= force[:,1,:] .= force[:,end,:].=0
        div = rho .* (C*force[:,:,1] + force[:,:,2]*C') + (Cr*rho).*force[:,:,1]+ (rho*Cr').*force[:,:,2]
        drho .= D*(M*rho + rho*M') - div

        mean_rho = 1/m[i] * dV*reshape(rho,1,N)*grid_points
        dzi .= 1/(Gamma_0*m[i]) * (mean_rho' - zi)
    end

end


function solve(tmax=0.1; alg=nothing)
    p = construct() 
    uz0 = initialconditions(p)
    
    # Solve the ODE
    prob = ODEProblem(f,uz0,(0.0,tmax),p)
    @time sol = DifferentialEquations.solve(prob, alg, progress=true,save_everystep=true,save_start=false)
    
    return sol, p
end

function solveplot(tmax=0.1; alg=nothing)
    sol, p = solve(0.1; alg=nothing)

    (; domain, dx, dy) = p
    #PLOTTING
    x = domain[1,1]:dx:domain[1,2]
    y = domain[2,1]:dy:domain[2,2]
    p1 = plot_solution(sol[end].x[1][:,:,1], sol[end].x[2][:,1], x, y; title="ρ₁", label="z₁") 
    p2 = plot_solution(sol[end].x[1][:,:,2], sol[end].x[2][:,2], x, y; title="ρ₋₁", label="z₋₁")

    plot(p1, p2, layout=[1 1], size=(1000,400)) |> display
    savefig("finaltime_pde.png")
    return sol, p
end

function test_f()
    p = construct()
    uz0 = initialconditions(p)
    duz = copy(uz0)
    @time f(duz, uz0, p, 0)
    return duz
end

function creategif(sol,p,  tmax=0.1, dt=0.01)

    rho1(t)=sol(t).x[1][:,:,1]
    rho2(t)=sol(t).x[1][:,:,2]
    z1(t)=sol(t).x[2][:,1]
    z2(t)=sol(t).x[2][:,2]

    (; domain, dx, dy) = p
    x = domain[1,1]:dx:domain[1,2]
    y = domain[2,1]:dy:domain[2,2]
    clims = (0, maximum(maximum(sol(t).x[1]) for t in 0:dt:tmax)) #limits colorbar
    pdegif = @animate for t = 0:dt:tmax
        p1 = plot_solution(rho1(t), z1(t), x, y; title="ρ₁", label="z₁",clims) 
        p2 = plot_solution(rho2(t), z2(t), x, y; title="ρ₋₁", label="z₋₁",clims)
        plot(p1, p2, layout=[1 1], size=(1000,400)) 
    end
     
    Plots.gif(pdegif, "evolution.gif", fps = 30)
end

function plot_solution(rho, z, x, y; title="", label="", clims=`:auto`)
    subp =    heatmap(x,y, rho,title = title, c=:berlin, clims=clims)
    scatter!(subp, [z[1]], [z[2]], markercolor=[:yellow],markersize=6, lab=label)
    return subp
end

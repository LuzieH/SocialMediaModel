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
    force = reshape(pointwise_int,N_y,N_x, 2)
    return force
end


function second_derivative((N_x,N_y), (dx, dy))
    Mx = Array(Tridiagonal([1.0 for i in 1:N_x-1],[-2.0 for i in 1:N_x],[1.0 for i in 1:N_x-1]))
    My = copy(Mx)
    # here matrix should have different shape for Ny not equal Nx
    # neumann boundary conditions with zero flux imply (using centered difference) that value of function outside domain
    # equals value of function inside domain (not on boundary), thus the following replacement
    Mx[2,1] = 2.0
    Mx[end-1,end] = 2.0
    My[1,2] = 2.0
    My[end,end-1] = 2.0
    Mx .= (1/dx^2)*Mx
    My .= (1/dy^2)*My
    return Mx, My
end

function centered_difference_force((N_x,N_y), (dx, dy))
    #centered first difference for force
    Cx = Array(Tridiagonal([1.0 for i in 1:N_x-1],[0. for i in 1:N_x],[-1.0 for i in 1:N_x-1]))
    Cy = -copy(Cx)
    Cx .= 1/(2*dx)*Cx
    Cy .= 1/(2*dy)*Cy
    Cy[1,1:2] = 1/(dy)* [-1,1]
    Cy[end,end-1:end] = 1/(dy)* [-1,1]
    Cx[1:2,1] = 1/(dx)* [-1,1]
    Cx[end-1:end,end] = 1/(dx)* [-1,1]
    return Cx, Cy
end

function centered_difference_density((N_x,N_y), (dx, dy))
    #centered first difference for density
    Cx_rho = Array(Tridiagonal([1.0 for i in 1:N_x-1],[0. for i in 1:N_x],[-1.0 for i in 1:N_x-1]))
    Cy_rho = -copy(Cx_rho)
    Cx_rho .= 1/(2*dx)*Cx_rho
    Cy_rho .= 1/(2*dy)*Cy_rho
    Cx_rho[2,1]=0. ; Cx_rho[end-1,end]=0.; Cy_rho[1,2]=0.; Cy_rho[end,end-1]=0.
    return Cx_rho, Cy_rho
end

function construct()
    # Define the constants for the PDE
    sigma = 1.0
    D = sigma^2 * 0.5
    a = 2.0
    c = 2.0
    Gamma_0 = 100
    m_1 = 1/2 #proportion of agents of type 1,2
    m_2 = 1/2 
    dx = 0.25
    dy = dx
    dV = dx*dy
    # domain = [[-2,2],[-2,2]]
    N_x = Int(4/dx+1)
    N_y = N_x #so far only works if N_y = N_x
    N = N_x*N_y
    X = reshape([x for x in -2:dx:2 for y in -2:dy:2],N_y,N_x)
    Y = reshape([y for x in -2:dy:2 for y in -2:dx:2],N_y,N_x)
    grid_points = [reshape(X,1,:);reshape(Y,1,:)]' 
    # matrix of K evaluated at gridpoints
    K_matrix, W_matrix  = generate_K(N, grid_points)
    
    Mx, My = second_derivative((N_x,N_y), (dx, dy))
    Cx, Cy = centered_difference_force((N_x,N_y), (dx, dy))
    Cx_rho, Cy_rho = centered_difference_density((N_x,N_y), (dx, dy))


    # Define the discretized PDE as an ODE function
    function f(duz,uz,p,t)

        grid_points, N_x, N_y, a, c, K_matrix, W_matrix, dV, Cx, Cy, Cx_rho, Cy_rho, D, Mx, My, N, m_1, m_2, Gamma_0 = p
        u, z = uz.x
        du, dz = duz.x

        rho_1 = @view  u[:,:,1]
        rho_2 = @view  u[:,:,2]
        drho_1 = @view du[:,:,1]
        drho_2 = @view du[:,:,2]

        z_1 = @view z[:,1]
        z_2 = @view z[:,2]
        dz_1 = @view dz[:,1]
        dz_2 = @view dz[:,2]

        #@show size(media_interation(rho_1, z_1))

        #@show size(z_1)

        rho = rho_1 + rho_2
        force_1 = c * media_force(z_1, grid_points, N_x, N_y) + a * agent_force(rho, K_matrix, W_matrix, dV)
        force_2 = c * media_force(z_2, grid_points, N_x, N_y) + a * agent_force(rho, K_matrix, W_matrix, dV)
        div_1 = rho_1 .* (Cy*force_1[:,:,2] + force_1[:,:,1]*Cx) + (Cy_rho*rho_1)*force_1[:,:,2]+ (rho_1*Cx_rho)*force_1[:,:,1]
        div_2 = rho_2 .* (Cy*force_2[:,:,2] + force_2[:,:,1]*Cx) + (Cy_rho*rho_2)*force_2[:,:,2]+ (rho_2*Cx_rho)*force_2[:,:,1]
        drho_1 .= D*(My*rho_1 + rho_1*Mx) - div_1
        drho_2 .= D*(My*rho_2 + rho_2*Mx) - div_2
        mean_rho_1 = 1/m_1 * reshape(rho_1,1,N)*grid_points
        mean_rho_2 = 1/m_2 * reshape(rho_2,1,N)*grid_points
        dz_1 .= 1/(Gamma_0*m_1) * (mean_rho_1' - z_1)
        dz_2 .= 1/(Gamma_0*m_2) * (mean_rho_2' - z_2)
    end
    
    p = (grid_points, N_x, N_y, a, c, K_matrix, W_matrix, dV, Cx, Cy, Cx_rho, Cy_rho, D, Mx, My, N, m_1, m_2, Gamma_0)
    return p, X, Y

end

function initialconditions(N_x = 17, N_y = 17)
    rho_0 = cat(1/32*ones(N_y,N_x), 1/32*ones(N_y,N_x), dims=3)
    u0 = rho_0
    z1_0 = [1.,1.]
    z2_0 = [-1.,-1.]
    z0 = [z1_0  z2_0]
    return ArrayPartition(u0,z0)
end

function f(duz,uz,p,t)
    grid_points, N_x, N_y, a, c, K_matrix, W_matrix, dV, Cx, Cy, Cx_rho, Cy_rho, D, Mx, My, N, m_1, m_2, Gamma_0 = p

    u, z = uz.x
    du, dz = duz.x

    rho_1 = @view  u[:,:,1]
    rho_2 = @view  u[:,:,2]
    drho_1 = @view du[:,:,1]
    drho_2 = @view du[:,:,2]

    z_1 = @view z[:,1]
    z_2 = @view z[:,2]
    dz_1 = @view dz[:,1]
    dz_2 = @view dz[:,2]

    #@show size(media_interation(rho_1, z_1))

    #@show size(z_1)

    rho = rho_1 + rho_2
    force_1 = c * media_force(z_1, grid_points, N_x, N_y) + a * agent_force(rho, K_matrix, W_matrix, dV)
    force_2 = c * media_force(z_2, grid_points, N_x, N_y) + a * agent_force(rho, K_matrix, W_matrix, dV)
    div_1 = rho_1 .* (Cy*force_1[:,:,2] + force_1[:,:,1]*Cx) + (Cy_rho*rho_1)*force_1[:,:,2]+ (rho_1*Cx_rho)*force_1[:,:,1]
    div_2 = rho_2 .* (Cy*force_2[:,:,2] + force_2[:,:,1]*Cx) + (Cy_rho*rho_2)*force_2[:,:,2]+ (rho_2*Cx_rho)*force_2[:,:,1]
    drho_1 .= D*(My*rho_1 + rho_1*Mx) - div_1
    drho_2 .= D*(My*rho_2 + rho_2*Mx) - div_2
    mean_rho_1 = 1/m_1 * reshape(rho_1,1,N)*grid_points
    mean_rho_2 = 1/m_2 * reshape(rho_2,1,N)*grid_points
    dz_1 .= 1/(Gamma_0*m_1) * (mean_rho_1' - z_1)
    dz_2 .= 1/(Gamma_0*m_2) * (mean_rho_2' - z_2)
end


function solve(tmax=0.01)
    uz0 = initialconditions()
    p, X, Y = construct()
    
    # Solve the ODE
    prob = ODEProblem(f,uz0,(0.0,tmax),p)
    
    @time sol = DifferentialEquations.solve(prob,progress=true,save_everystep=true,save_start=false)

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

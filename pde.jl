using DifferentialEquations, LinearAlgebra
using RecursiveArrayTools
using Plots

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
        return x
    else
        return 0.1
    end
end

function gamma(grid_points, rho, y, beta,dV)
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
    return beta*reshape(rate,(Nx, Ny, I, J, J))
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

function construct()
    # Define the constants for the PDE
    sigma = 1.0
    D = sigma^2 * 0.5
    a = 2.0
    b = 2.0
    c = 2.0
    beta = 100.0
    Gamma_0 = 100
    gamma_0 = 25
    J=4 # number of influencers
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
 
    p = (; grid_points, N_x, N_y, domain,  a, b, c, beta, K_matrix, W_matrix, dx,dy, dV, C,  D, M , N, J, Gamma_0, gamma_0)

    return p
end

function initialconditions(p)
    (; N_x , N_y,  J) = p
    rho_0 = zeros(N_x, N_y, 2, 4) 
    mid_y =Int(round(N_y/2))
    mid_x =Int(round(N_x/2))
    rho_0[1:mid_x, 1:mid_y,:,1] .= 1/(2*J*4)
    rho_0[1:mid_x, mid_y+2:end,:,2] .= 1/(2*J*4)
    rho_0[mid_x+2:end, 1:mid_y,:,3] .= 1/(2*J*4)
    rho_0[mid_x+2:end, mid_y+2:end,:,4] .= 1/(2*J*4)
    #cat(1/32*ones(N_x,N_y), 1/32*ones(N_x,N_y), dims=3)
    #rho_0 = fill(1/(2*J*4*4), N_x, N_y, 2, J)
    u0 = rho_0
    z1_0 = [1.,1.]
    z2_0 = [-1.,-1.]
    z0 = [z1_0  z2_0]
    y1_0 = [-1.,-1.]
    y2_0 = [-1.,1.]
    y3_0 = [1.,-1.]
    y4_0 = [1.,1. ]

    y0 = [y1_0  y2_0 y3_0 y4_0]
    return ArrayPartition(u0,z0,y0)
end

function f(duzy,uzy,p,t)
    yield()
    (; grid_points, N_x, N_y, domain,  a, b, c, beta, K_matrix, W_matrix, dx,dy, dV, C,  D, M , N, J, Gamma_0, gamma_0) = p
    u, z, y2 = uzy.x
    du, dz, dy2 = duzy.x

    rhosum = sum(u, dims=(3,4))[:,:,1,1]
    rhosum_j = sum(u, dims=3)
    rhosum_i = sum(u, dims=4)
    m_i = dV*sum(u,dims = (1,2,4))
    m_j = dV*sum(u,dims=(1,2,3))
    Fagent = agent_force(rhosum, K_matrix, W_matrix, dV)
    rate_matrix = gamma(grid_points, u, y2, beta,dV)
    for i in 1:2
        for j in 1:J
            rho = @view  u[:,:,i,j]
            drho = @view du[:,:,i,j]
            zi = @view z[:,i]
            dzi = @view dz[:,i]
            yj = @view y2[:,j]
            dyj = @view dy2[:,j]

            force = c * follower_force(zi, grid_points, N_x, N_y) + a * Fagent + b * follower_force(yj, grid_points, N_x, N_y)
            ##ensure force is zero at boundary to conserve density - effectivly together with transparent Neumann BCs, now the fluxes at boundaries cancel
            #force[1,:,:] .= force[end,:,:] .= force[:,1,:] .= force[:,end,:].=0            
            
            div =  C * (rho .* force[:,:,1]) + (rho .* force[:,:,2]) * C'
            reac = zeros(N_x, N_y)
            for j2=1:J
                if j2!= j
                    reac += -rate_matrix[:,:,i,j,j2] .* rho + rate_matrix[:,:,i,j2,j] .* u[:,:,i,j2]
                end
            end
            
            dif = D*(M*rho + rho*M')  
            #balance fluxes at boundary
            dif[1,:]+= -D/dx * (force[1,:,1].*rho[1,:]) 
            dif[end,:]+= D/dx * (force[end,:,1].*rho[end,:])
            dif[:,1]+= -D/dy * (force[:,1,2].*rho[:,1])
            dif[:,end]+= D/dy * (force[:,end,2].*rho[:,end])

            drho .=  dif - div + reac

            mean_rhoi = 1/m_i[i] * dV*reshape(rhosum_i[:,:,i,:],1,N)*grid_points
            dzi .= 1/(Gamma_0*m_i[i]) * (mean_rhoi' - zi)

            mean_rhoj = 1/m_j[j] * dV*reshape(rhosum_j[:,:,:,j],1,N)*grid_points
            dyj .= 1/(gamma_0*m_j[j]) * (mean_rhoj' - yj)
        end
    end

end

function solve(tmax=0.1; alg=nothing)
    p = construct() 
    uzy0 = initialconditions(p)
    
    # Solve the ODE
    prob = ODEProblem(f,uzy0,(0.0,tmax),p)
    @time sol = DifferentialEquations.solve(prob, alg, progress=true,save_everystep=true,save_start=false)
    
    return sol, p
end

function plot_solution(rho, z, y, x_arr, y_arr; title="", labelz="", labely="", clim=(-Inf, Inf))
    subp = heatmap(x_arr,y_arr, rho', title = title, c=:berlin, clims=clim)
    scatter!(subp, [z[1]], [z[2]], markercolor=[:yellow],markersize=6, lab=labelz)
    scatter!(subp, [y[1]], [y[2]], markercolor=[:red],markersize=6, lab=labely)
    return subp
end

function solveplot(tmax=0.1; alg=nothing)
    sol, p = solve(tmax; alg=nothing)

    (; domain, dx, dy,J) = p
    #PLOTTING
    x_arr = domain[1,1]:dx:domain[1,2]
    y_arr = domain[2,1]:dy:domain[2,2]
    
    plot_array = Any[]  
    z_labels = ["z₁", "z₋₁"]
    y_labels = ["y₁", "y₂", "y₃", "y₄"]
    dens_labels = [ "ρ₁₁" "ρ₁₂" "ρ₁₃" "ρ₁₄"; "ρ₋₁₁" "ρ₋₁₂" "ρ₋₁₃" "ρ₋₁₄"]
    for j in 1:J    
        for i in 1:2
            # make a plot and add it to the plot_array
            push!(plot_array, plot_solution(sol(tmax).x[1][:,:,i,j], sol(tmax).x[2][:,i], sol(tmax).x[3][:,j], x_arr, y_arr; title = string(dens_labels[i,j],"(", string(tmax), ")") ,labely = y_labels[j], labelz = z_labels[i]))
        end
    end
    plot(plot_array..., layout=(4,2),size=(1000,1000)) |> display

    #p2 = plot_solution(sol(tmax).x[1][:,:,2], sol(tmax).x[2][:,2], x_arr, y_arr; title=string("ρ₋₁(",string(tmax),")"), label="z₋₁")

    #plot(p1, p2, layout=[4 4], size=(1000,4*400)) 
    savefig("finaltime_pde.png")
    return sol, p
end
 
function creategif(sol,p, dt=0.01)
    #rho1(t)=sol(t).x[1][:,:,1]
    #'rho2(t)=sol(t).x[1][:,:,2]
    #z1(t)=sol(t).x[2][:,1]
    #z2(t)=sol(t).x[2][:,2]
    
    tmax=sol.t[end]
    (; domain, dx, dy, J) = p
    x_arr = domain[1,1]:dx:domain[1,2]
    y_arr = domain[2,1]:dy:domain[2,2]
    cl = (0, maximum(maximum(sol(t).x[1]) for t in 0:dt:tmax)) #limits colorbar
    z_labels = ["z₁", "z₋₁"]
    y_labels = ["y₁", "y₂", "y₃", "y₄"]
    dens_labels = [ "ρ₁₁" "ρ₁₂" "ρ₁₃" "ρ₁₄"; "ρ₋₁₁" "ρ₋₁₂" "ρ₋₁₃" "ρ₋₁₄"]

    pdegif = @animate for t = 0:dt:tmax
        plot_array = Any[]  
        for j in 1:J    
            for i in 1:2
                # make a plot and add it to the plot_array
                push!(plot_array, plot_solution(sol(t).x[1][:,:,i,j], sol(t).x[2][:,i], sol(t).x[3][:,j], x_arr, y_arr; title = string(dens_labels[i,j],"(", string(t), ")") ,labely = y_labels[j], labelz = z_labels[i], clim = cl))
            end
        end
        plot(plot_array..., layout=(4,2),size=(1000,1000)) |> display
    end
    Plots.gif(pdegif, "evolution.gif", fps = 30)
end

function test_f()
    p = construct()
    uzy0 = initialconditions(p)
    duzy = copy(uzy0)
    @time f(duzy, uzy0, p, 0)
    return duzy
end
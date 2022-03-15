using DifferentialEquations, LinearAlgebra
using RecursiveArrayTools
 
# Define the constants for the PDE
const sigma = 1.0
const D = sigma^2 * 0.5
const a = 2.0
const c = 2.0
const Gamma_0 = 100
const dx = 0.1
const dy = 0.1
const dV = dx*dy
# domain = [[-2,2],[-2,2]]
const N_x = Int(4/dx+1)
const N_y = Int(4/dy+1) #so far only works if N_y = N_x
const N = N_x*N_y
const X = reshape([x for x in -2:dx:2 for y in -2:dy:2],N_y,N_x)
const Y = reshape([y for x in -2:dy:2 for y in -2:dx:2],N_y,N_x)
const grid_points = [reshape(X,1,:);reshape(Y,1,:)]' 

w(s,x) = exp(-norm(s-x))
K(s,x) = w(s,x) * (s-x)    

function generate_K()
    k = zeros(N, N, 2)
    for i in 1:N, j in 1:N
        k[i,j,:] = K(grid_points[i,:], grid_points[j,:])
    end
    return k
end

# matrix of K evaluated at gridpoints
const K_matrix  = generate_K()
 
# Define the initial condition as normal arrays
rho_0 = cat(1/32*ones(N_y,N_x), 1/32*ones(N_y,N_x), dims=3)
u0 = rho_0
z1_0 = [1.,1.]
z2_0 = [-1.,-1.]
z0 = [z1_0  z2_0]
const m_1 = 1/2 #proportion of agents of type 1,2
const m_2 = 1/2 

# interaction with other agents
function agent_force(rho)
   # @show size(rho), N
    #pointwise_int = dV*sum(vec(rho).*K_matrix, dims=1)

    pointwise_int = zeros(N, 2)
    for d in 1:2, i in 1:N, j in 1:N
        pointwise_int[i,d] += rho[j] * K_matrix[j, i, d]
    end
    pointwise_int .*= dV

    force = reshape(pointwise_int, N_y, N_x, 2)
    return force
end

# interaction with media
function media_force(z)
    pointwise_int =   z' .- grid_points
    force = reshape(force,N_y,N_x, 2)
    return force
end

const Mx = Array(Tridiagonal([1.0 for i in 1:N_x-1],[-2.0 for i in 1:N_x],[1.0 for i in 1:N_x-1]))
const My = copy(Mx)
# neumann boundary conditions with zero flux imply (using centered difference) that value of function outside domain
# equals value of function inside domain (not on boundary), thus the following replacement
Mx[2,1] = 2.0
Mx[end-1,end] = 2.0
My[1,2] = 2.0
My[end,end-1] = 2.0
Mx .= (1/dx^2)*Mx
My .= (1/dy^2)*My

#centered first difference for force
const Cx = Array(Tridiagonal([1.0 for i in 1:N_x-1],[0. for i in 1:N_x],[-1.0 for i in 1:N_x-1]))
const Cy = -copy(Cx)
Cx .= 1/(2*dx)*Cx
Cy .= 1/(2*dy)*Cy
Cy[1,1:2] = 1/(dy)* [-1,1]
Cy[end,end-1:end] = 1/(dy)* [-1,1]
Cy[1:2,1] = 1/(dx)* [-1,1]
Cy[end-1:end,end] = 1/(dx)* [-1,1]

#centered first difference for density
const Cx_rho = Array(Tridiagonal([1.0 for i in 1:N_x-1],[0. for i in 1:N_x],[-1.0 for i in 1:N_x-1]))
const Cy_rho = -copy(Cx)
Cx_rho .= 1/(2*dx)*Cx
Cy_rho .= 1/(2*dy)*Cy
Cx_rho[2,1]=0. ; Cx_rho[end-1,end]=0.; Cy_rho[1,2]=0.; Cy_rho[end,end-1]=0.

# Define the discretized PDE as an ODE function
function f(duz,uz,p,t)

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

  force_1 = c*media_force(z_1) + a*agent_force(rho_1+ rho_2)
  force_2 = c*media_force(z_2) + a*agent_force(rho_1+ rho_2)
  div_1 = rho_1 .* (Cy*force_1[:,:,2] + force_1[:,:,1]*Cx) 
  div_2 = rho_2 .* (Cy*force_2[:,:,2] + force_2[:,:,1]*Cx)
  drho_1 .= D*(My*rho_1 + rho_1*Mx) - div_1
  drho_2 .= D*(My*rho_2 + rho_2*Mx) - div_2
  mean_rho_1 = 1/m_1 * reshape(rho_1,1,N)*grid_points
  mean_rho_2 = 1/m_2 * reshape(rho_2,1,N)*grid_points
  dz_1 .= 1/(Gamma_0*m_1) * (mean_rho_1' - z_1)
  dz_2 .= 1/(Gamma_0*m_2) * (mean_rho_2' - z_2)
end
 
uz0 = ArrayPartition(u0,z0)


# Solve the ODE
prob = ODEProblem(f,uz0,(0.0,0.01))
@time sol = solve(prob,progress=true,save_everystep=true,save_start=false)

using Plots;
#pyplot()
p1 = surface(X,Y,sol[end].x[1][:,:,1],title = "rho_1")
p2 = surface(X,Y,sol[end].x[1][:,:,2],title = "rho_{-1}")

plot(p1,p2,layout=grid(2,1))

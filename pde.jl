using DifferentialEquations, LinearAlgebra
using RecursiveArrayTools
using Plots
using JLD2
using Distances
pyplot()


function parameters(;
        J=4,  #number of influencers
        n_media = 2, #number of media
        n = 250, #number of agents 
        eta = 15.0, #rate constant for changing influencer  
        controlspeed = 0.3, 
        a = 1., #interaction strength between agents
        b = 4., # interaction strength between agents and influencers
        c = 2., # interaction strength between agents and media 
        sigma = 0.5, # noise on individual agents
        sigmahat = 0., # noise on influencers
        sigmatilde = 0., # noise on media
        frictionI = 10., # friction for influencers
        frictionM = 100.,  #friction for media
        controltarget = [1.5 1.5]
    )

    q = (; n, J, n_media, frictionM, frictionI, a, b, c, eta, sigma, sigmahat, sigmatilde, controlspeed,controltarget)
    return q
end

function parameters_control(;controlspeed = 0.3)
    q = parameters(eta=1.,controlspeed = controlspeed) 
    return q
end

function ABMconstruct(;
        # setting simulation parameters
        dt = 0.01,  # simulation stepsize 
        dx = 0.05,
        dy = dx,
        domain =  2.2*[-1 1; -1 1]
    )
    X = [x for x in domain[1,1]:dx:domain[1,2], y in domain[2,1]:dy:domain[2,2]]
    Y = [y for x in domain[1,1]:dx:domain[1,2], y in domain[2,1]:dy:domain[2,2]]
    dV = dx*dy

    p = (; dt, dx, dy, domain, X, Y, dV)
    return p
end

function PDEconstruct(;
        # Define the constants for the PDE
        dx = 0.05, #0.05
        dy = dx,
        domain = 2.2*[-1 1; -1 1]
    )
    N_x = Int(round((domain[1,2]-domain[1,1])/dx+1))
    N_y = Int(round((domain[2,2]-domain[2,1])/dy+1)) #so far only works if N_y = N_x
    N = N_x*N_y
    dV = dx*dy # to integrate the agent distribution
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


function PDEconstructcoarse()
    return PDEconstruct(;dx = 0.1, dy = 0.1,domain = 2.3*[-1 1; -1 1])
end

function PDEconstructmeso() #for final optimal control computations
    return PDEconstruct(;dx = 44/600, dy = 44/600, domain = 2.2*[-1 1; -1 1])
end

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
    return ArrayPartition(u0,z0,y0), counts, controlled
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
    return ArrayPartition(u0,z0,y0), counts, controlled
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
    return ArrayPartition(u0,z0,y0), counts, controlled
end

function constructinitial(scenario,(p,q))
    if scenario=="4inf"
        uzy0, counts, controlled = inf_initialconditions((p,q))
    elseif scenario=="noinf"
        J=1
        b=0
        eta=0
        q=(;q..., J,b,eta)
        uzy0, counts, controlled = noinf_initialconditions((p,q))
    elseif scenario =="uniform"
        uzy0, counts, controlled = uniform_initialconditions((p,q))
    elseif scenario=="fixedtarget"
        uzy0, counts, controlled = uniform_initialconditions((p,q))
    elseif scenario=="localmax"
        uzy0, counts, controlled = uniform_initialconditions((p,q))
    elseif scenario=="optimalcontrol"
        uzy0, counts, controlled = uniform_initialconditions((p,q))
    end
    return uzy0, counts, controlled,q
end


function f(duzy,uzy,(p,q),t)
    yield()
    (; grid_points, N_x, N_y,N, K_matrix, W_matrix, dx,dy, dV, C,  M) = p
    (; a, b, c, sigma, eta,  J, frictionM, frictionI, controlled, controlspeed,controltarget) = q
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

            mean_rhoi = 1/m_i[i] * dV*reshape(rhosum_i[:,:,i,:],1,N)*grid_points
            dzi .= 1/(frictionM) * (mean_rhoi' - zi)

            if controlled[j] == 0
                mean_rhoj = 1/m_j[j] * dV*reshape(rhosum_j[:,:,:,j],1,N)*grid_points
                dyj .= 1/(frictionI) * (mean_rhoj' - yj)
            
            # elseif controlled[j]==2
            #     #local maximization of followers
            #     speed = 0.1

            #     grad_rate_x  = (1/dx)*(gamma(grid_points, u, y2 .+ [dx 0], eta, dV) - gamma(grid_points, u, y2, eta, dV))
            #     grad_rate_y = (1/dy)*(gamma(grid_points, u, y2 .+ [0 dy], eta, dV) - gamma(grid_points, u, y2, eta, dV))
                
            #     gradC_x  = sum(grad_rate_x[:,:,:,1,j] .* (rhosum_j- u[:,:,:,j]))*dV
            #     gradC_y  = sum(grad_rate_y[:,:,:,1,j] .* (rhosum_j- u[:,:,:,j]))*dV
            #     dyj .= speed* [gradC_x gradC_y]
            #     #todo check if it works, maybe add noise to escape local minima 

            elseif controlled[j] == 1 #controll movement
                if norm(controltarget' - yj) >0
                    dyj .= controlspeed* (controltarget' - yj)./ norm(controltarget' - yj)
                else
                    dyj .= 0
                end
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

    uzy0, counts, controlled,q = constructinitial(scenario,(p,q))
    q = (; q..., controlled)

    # Solve the ODE
    prob = ODEProblem(f,uzy0,(0.0,tmax),(p,q))
    @time sol = DifferentialEquations.solve(prob, alg, save_start=true)

    return sol, (p,q), counts
end

function solvefixedtargets(tequil = 5., tmax=10.1, targets = [[1.5 1.5]]; alg=nothing, scenario="fixedtarget", p = PDEconstruct(), q= parameters_control(), savedt=0.05, atol = 1e-6, rtol = 1e-3)


    #solve until tequil
    uzy0, _, controlled,q = constructinitial(scenario,(p,q))
    q1 = (; q..., controlled)

    # Solve the ODE
    prob1 = ODEProblem(f,uzy0,(0.0,tequil),(p,q1))

    @time sol1 = DifferentialEquations.solve(prob1, alg, saveat = 0:savedt:tequil,save_start=true, abstol = atol, reltol = rtol)
    sols = [sol1]
    Ps =[(p,q1)]

    n_targets = size(targets,1)
    ttotalcontrol = tmax - tequil
    tindcontrol = ttotalcontrol/n_targets

    #add new influencer
    u,z,y = sol2uyz(sol1,tequil)
    q2 = merge(q1, (;J=q1.J+1,controlled = [q1.controlled..., 1]))
    (;N_x, N_y, grid_points,N, dV) = p
    J= q2.J

    
    u2 = zeros(N_x, N_y, 2, J)
    u2[:,:,:,1:J-1] = u
    startlocation =  1/sum(u2,dims=(1,2,4))[1,1,2,1] * reshape(sum(u2, dims=4)[:,:,2,:],1,N)*grid_points #[0 0]
    y2 = zeros(2, J)
    y2[:,1:J-1] = y
    #TODO maybe make constant speed such that starts and ends within 5. time steps
    y2[:,J] = startlocation 
    uzy0 = ArrayPartition(u2,z,y2)
    followersum=0

    for n in 1:n_targets
        target = targets[n]
        speed = norm(target - startlocation)/tindcontrol 
        q3 = merge(q2, (;controltarget = target, controlspeed = speed))

        # solve ODE with added influencer
        prob2 = ODEProblem(f,uzy0,(0.0,tindcontrol ),(p,q3))
        @time sol2 = DifferentialEquations.solve(prob2, alg,  saveat = 0:savedt:tindcontrol ,save_start=true, abstol = atol, reltol = rtol)


        followersum += sum([sum(sol2(t).x[1][:,:,:,J]) for t in sol2.t])*dV*savedt 
        push!(sols,sol2)
        push!(Ps, (p,q3))
        #prepare next simulation
        startlocation = target
        uzy0 = sol2(tindcontrol)
    end


    return sols, Ps,  followersum
end

# function solvelocalmax(tequil = 5., tmax=10.1; alg=nothing, scenario="localmaximization", p = PDEconstruct(), q= parameters_control(), savedt=0.05, atol = 1e-6, rtol = 1e-3)
#     uzy0, _, controlled,q = constructinitial(scenario,(p,q))
#     q1 = (; q..., controlled)

#     # Solve the ODE
#     prob1 = ODEProblem(f,uzy0,(0.0,tequil),(p,q1))

#     @time sol1 = DifferentialEquations.solve(prob1, alg, saveat = 0:savedt:tequil,save_start=true, abstol = atol, reltol = rtol)

#     #add new influencer
#     startlocation =  [0 0] #1/sum(u2,dims=(1,2,4))[1,1,2,1] * reshape(sum(u2, dims=4)[:,:,2,:],1,N)*grid_points
#     u,z,y = sol2uyz(sol1,tequil)
#     q2 = merge(q1, (;J=q1.J+1,controlled = [q1.controlled..., 2]))
#     (;N_x, N_y, grid_points,N, dV) = p
#     J= q2.J
#     u2 = zeros(N_x, N_y, 2, J)
#     u2[:,:,:,1:J-1] = u
#     y2 = zeros(2, J)
#     y2[:,1:J-1] = y
#     #TODO maybe make constant speed such that starts and ends within 5. time steps
#     y2[:,J] = startlocation 
#     uzy0 = ArrayPartition(u2,z,y2)

#     # solve ODE with added influencer
#     prob2 = ODEProblem(f,uzy0,(0.0,tmax-tequil),(p,q2))
#     @time sol2 = DifferentialEquations.solve(prob2, alg,  saveat = 0:savedt:(tmax-tequil),save_start=true, abstol = atol, reltol = rtol)
#     followersum = sum([sum(sol2(t).x[1][:,:,:,J]) for t in sol2.t])*dV*savedt 
#     return [sol1, sol2], [(p,q1), (p,q2)],  followersum
# end

function controlsearch(tequil = 5., tcontrol = 2.5, idx = 0.75, ibound = 1.5; alg=nothing, scenario="optimalcontrol", p = PDEconstructcoarse(), q= parameters_control(), savedt=0.05, atol = 1e-6, rtol = 1e-3)

    #influencer grid points that can be reached at the end of each subinterval
    Xi = [x for x in -ibound:idx:ibound, y in -ibound:idx:ibound]
    Yi = [y for x in -ibound:idx:ibound, y in -ibound:idx:ibound]
    Ngrid = size(Xi,1) #per dimension
    followersum = zeros(Ngrid,Ngrid,Ngrid,Ngrid)

    uzy0, _, controlled,q = constructinitial(scenario,(p,q))
    qc = (; q..., controlled)

    # Solve the ODE
    prob1 = ODEProblem(f,uzy0,(0.0,tequil),(p,qc))

    @time sol1 = DifferentialEquations.solve(prob1, alg, saveat = 0:savedt:tcontrol,save_start=true, abstol = atol, reltol = rtol)
    # todo: maybe coarser savedt and compute follower integral on the fly

    #add new influencer
    u,z,y = sol2uyz(sol1,tequil)
    qadded = merge(qc, (;J=qc.J+1,controlled = [qc.controlled..., 1] ))
    (;N_x, N_y, grid_points,N, dV) = p
    J= qadded.J
    uadded = zeros(N_x, N_y, 2, J)
    uadded[:,:,:,1:J-1] = u
    yadded = zeros(2, J)
    yadded[:,1:J-1] = y
    # start location of new influencer
    # todo: maybe place on symmetry axis and then save computations below due to symmetry
    start =  1/sum(uadded,dims=(1,2,4))[1,1,2,1] * reshape(sum(uadded, dims=4)[:,:,2,:],1,N)*grid_points #this should already be symmetric start point
    # [0 0] #this adds another point symmetry
    yadded[:,J] = start
    uzy0 = ArrayPartition(uadded,z,yadded)

    # solve ODE with added influencer by looping over first endpoint of influencer at time = tcontrol
    for i in 1:Ngrid; for j in 1:Ngrid
        # only compute one side of the symmetry axis
        if i >=j
            #define controlspeed and target
            target1 = [Xi[i,j] Yi[i,j]]
            speed1 = norm(target1 - start)/tcontrol 
            # add to parameter set
            q1 = merge(qadded, (;controltarget = target1, controlspeed = speed1))

            prob2 = ODEProblem(f,uzy0,(0.0,tcontrol),(p,q1))
            @time sol2 = DifferentialEquations.solve(prob2, alg,  saveat = 0:savedt:tcontrol,save_start=true, abstol = atol, reltol = rtol)
            #final states
            uzy1 = sol2(tcontrol)
            followersum[i,j,:,:] .= sum([sum(sol2(t).x[1][:,:,:,J]) for t in sol2.t])*dV*savedt   
            
            #loop over second endpoints starting from previous endpoint
            #Threads.@threads 
            for k in 1:Ngrid; for l in 1:Ngrid
                #define controlspeed and target
                target2 = [Xi[k,l] Yi[k,l]]
                speed2 = norm(target2 - target1)/tcontrol 
                # add to parameter set
                q2 = merge(qadded, (;controltarget = target2, controlspeed = speed2))

                prob3 = ODEProblem(f,uzy1,(0.0,tcontrol),(p,q2))
                @time sol3 = DifferentialEquations.solve(prob3, alg,  saveat = 0:savedt:tcontrol,save_start=true, abstol = atol, reltol = rtol)
                #final states
                followersum[i,j,k,l] += sum([sum(sol3(t).x[1][:,:,:,J]) for t in sol3.t])*dV*savedt              
            end; end
        end
    end;end
    #fill remaining symmetric values
    for i in 1:Ngrid; for j in 1:Ngrid
        # only compute one side of the symmetry axis
        if i <j
            followersum[i,j,:,:] = followersum[j,i,:,:]' #also the second targets have to be mirrored
        end
    end; end
    @save string("data/optimalcontrol.jld2") p qadded followersum Xi Yi
    return Xi, Yi, followersum
end

function plotfollowersum(followersum, Xi, Yi)
    Ngrid = size(followersum,1)
    #average out second targets
    follower_av1 = dropdims(sum(followersum,dims=(3,4)) *1/(Ngrid^2), dims = (3,4))
    #average out first targets
    follower_av2 = dropdims(sum(followersum, dims=(1,2))*1/(Ngrid^2), dims = (1,2))
    heatmap(Xi[:,1],Xi[:,1],follower_av1,c=cmap, title="Objective function, second target averaged out",xlabel="x-position of target",ylabel="y-position of target")
    savefig("img/obj_secondaveragedout.png")
    heatmap(Xi[:,1],Xi[:,1],follower_av2,c=cmap,title="Objective function, first target averaged out",xlabel="x-position of target",ylabel="y-position of target")
    savefig("img/obj_firstaveragedout.png")
    #TODO plot paths with obj function color on opinion domain space
end

function solveplot(tmax=0.1; alg=nothing, scenario="4inf", p = PDEconstruct(), q= parameters())
    sol,(p,q), counts = solve(tmax; alg=alg, scenario=scenario, p=p, q=q)

    u,z,y = sol2uyz(sol, tmax)

    plt = plotarray(u,z,y, (p,q), tmax)
    plt |> display

    return sol, (p,q), counts
end

function vary_parameters_control(tcontrol = 5, tmax=40, savepoints = 4; alg=nothing, scenario="fixedtarget", p = PDEconstruct(), bs=4.0, speeds=0.2, frictions=5:2.5:10., etas =1.:0.5:3. ,  savedt=1., atol = 1e-3, rtol = 1e-2)
    (; N_x, N_y) = p
    J=5
    zs = zeros(2, 2, savepoints,size(etas,1), size(frictions,1), size(speeds,1),size(bs,1))
    ys = zeros(2, J, savepoints,size(etas,1), size(frictions,1), size(speeds,1),size(bs,1))
    us = zeros(N_x, N_y, 2, J, savepoints,size(etas,1), size(frictions,1), size(speeds,1),size(bs,1))
    savetimes = LinRange(0, tmax-tcontrol, savepoints)
    #Threads.@threads 
    for b in 1:size(bs,1)
        for s in 1:size(speeds,1)
            for f in 1:size(frictions,1)
                    for e in 1:size(etas,1)
                        q = parameters(b=bs[b], controlspeed=speeds[s], frictionI = frictions[f], eta = etas[e])
                        sols, Ps, _ = solvecontrolled(tcontrol, tmax; alg=alg, scenario=scenario, p=p, q=q,savedt=savedt, atol = atol, rtol= rtol)
                        sol = sols[2]
                        (p,q) = Ps[2]
                        anim = Animation()
                        for j in 1:savepoints
                            u,z,y = sol2uyz(sol, savetimes[j])
                            us[:,:,:,:,j,e,f,s,b] = u
                            zs[:,:,j,e,f,s,b] = z
                            ys[:,:,j,e,f,s,b] = y
                            plt = plotsingle(u,z,y, (p,q), j; save=false)
                            frame(anim, plt)
                        end
                        Plots.gif(anim, string("img/control_b",string(bs[b]),"speed",string(speeds[s]), "friction",string(frictions[f]), "eta",string(etas[e]),".gif"), fps = 10)
                    end
            end
        end
    end
    @save string("data/parameter_",scenario,".jld2") us zs ys bs speeds frictions etas
    return us,ys,zs
end

#=
for b in 1:size(bs,1)
    for s in 1:size(speeds,1)
        for fi in 1:size(frictions,1)
            anim = Animation()
            for t in 1:10
                plt = plotarray(us[:,:,:,:,t,fi,s,b],zs[:,:,t,fi,s,b],ys[:,:,t,fi,s,b], P, t; save=false)
                frame(anim, plt)
            end
            Plots.gif(anim, string("img/control_b",string(bs[b]),"speed",string(speeds[s]), "friction",string(frictions[fi]),".gif"), fps = 10)          
        end
    end
end
=#

function solveensemble(tmax=0.1, N=10; savepoints = 4, alg=nothing, scenario="4inf", p = PDEconstruct(), q= parameters())
    (; N_x, N_y) = p
    J=q.J

    zs = zeros(2, 2, savepoints, N)
    ys = zeros(2, J, savepoints, N)
    us = zeros(N_x, N_y, 2, J, savepoints,  N)
    savetimes = LinRange(0, tmax, savepoints)
    av_counts = zeros(J,2)
    #Threads.@threads 
    for i=1:N
        sol, _ , counts= solve(tmax; alg=alg, scenario=scenario, p=p, q=q)
        av_counts = av_counts +  counts*(1/N)
        for j in 1:savepoints
            u,z,y = sol2uyz(sol, savetimes[j])
            us[:,:,:,:,j,i] = u
            zs[:,:,j,i] = z
            ys[:,:,j,i] = y
        end

    end

    @save string("data/pde_ensemble_",scenario,".jld2") us zs ys
    return us, zs, ys, (p,q), av_counts
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

        if scenario=="4inf"
            (;X, Y, domain, dx, dy) = p
            J= q.J
            x_arr = domain[1,1]:dx:domain[1,2]
            y_arr = domain[2,1]:dy:domain[2,2]

            yall = reshape(ys[:,:,k,:], (2,N*J))'
            evalkde = [sumgaussian([X[i,j], Y[i,j]], yall) for i in 1:size(X,1), j in 1:size(X,2)]
            heatmap(x_arr, y_arr, evalkde', c=cmap, title=string("Distribution of influencers at time ", string(round(savetimes[k], digits=2)))) #|> display
            savefig(string(title2,string(k),".png"))
        end
    end


end

function psensemble(tmax=0.1, N=10; alg=nothing, scenario="4inf")
    us, zs, ys, (p,q), av_counts  = solveensemble(tmax, N; alg=alg, scenario=scenario)
    plotensemble(us, zs, ys, (p,q), tmax)
    return us, zs, ys, (p,q), av_counts
end

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


function plotsnapshots(sols::Vector, Ps::Vector, ts; save = true, scenario="4inf")
    n_snapshots = length(ts)
    plot_array = Any[]  
    sum_t = 0
    for (sol, P) in zip(sols, Ps)
        for t in ts
            if t>= sum_t && t< sum_t + sol.t[end]  
                u,z,y = sol2uyz(sol, t-sum_t)
                ylabel =string("t = ", string(round(t, digits=2)))
                title=""
                subp = plotsingle(u,z,y,(p,q),t,save=false, scenario=scenario, labely="", labelx="",ylabel=ylabel,title=title)
                push!(plot_array, subp)
            end
        end
        sum_t+= sol.t[end]
    end    
    gridp=plot(plot_array..., layout=(n_snapshots,1),size=(95*5,n_snapshots*50*5),link=:all)#, left_margin=10mm )

    for k=1:n_snapshots-1
        plot!(gridp[k],xformatter=_->"")
    end

    if save==true
        savefig(string("img/pde_snapshots_",scenario,".png"))
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
    J= q.J
    #u,z,y = sol2uyz(sol, t)

    x_arr = domain[1,1]:dx:domain[1,2]
    y_arr = domain[2,1]:dy:domain[2,2]

    dens = dropdims(sum(u, dims=(3,4)), dims=(3,4))

    subp = heatmap(x_arr,y_arr, dens', title = title,ylabel=ylabel, c=cmap,clim=clim,legend=legend )

    scatter!(subp, z[1,:], z[2,:], markercolor=colors_leaders[2],markersize=size_leaders_pde, lab = labelx)

    if scenario!="noinf"
        scatter!(subp, y[1,:], y[2,:], markercolor=colors_leaders[1],markersize=size_leaders_pde, lab=labely)
    end

    #subp # |> display

    if save==true
        savefig(string("img/pde_single_",scenario,".png"))
    end
    return subp
end

gifsingle(sol, (p,q), args...; kwargs...) = gifsingle([sol], [(p,q)], args...; kwargs...)

function gifsingle(sols::Vector, Ps::Vector, dt=0.1; scenario = "4inf")
    T = 0
    anim = Animation()
    for (sol, P) in zip(sols, Ps)
        for t in 0:dt:sol.t[end]
            u,z,y = sol2uyz(sol, t)
            plt = plotsingle(u,z,y,P,t+T,save=false, scenario=scenario)
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

 
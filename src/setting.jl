function parameters(;
        J=4,  #number of influencers
        n_media = 2, #number of media
        n = 250, #number of agents 
        eta = 15.0, #rate constant for changing influencer  
        a = 1., #interaction strength between agents
        b = 4., # interaction strength between agents and influencers
        c = 2., # interaction strength between agents and media 
        sigma = 0.5, # noise on individual agents
        sigmahat = 0., # noise on influencers
        sigmatilde = 0., # noise on media
        frictionI = 10., # friction for influencers
        frictionM = 100.,  #friction for media
        controlspeed1 = sqrt(1.5^2 +1.5^2)/5, #control definition of first influencer or medium
        controltarget1 = [1.5 1.5],
        controlspeed2 = sqrt(1.5^2 +1.5^2)/5, #control definition of another influencer
        controltarget2 = [1.5 1.5]
    )

    q = (; n, J, n_media, frictionM, frictionI, a, b, c, eta, sigma, sigmahat, sigmatilde, controlspeed1,controltarget1, controlspeed2,controltarget2)
    return q
end

function parametersstronginf()
    q = parameters(;eta=1.) 
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

function parameterscontrol(;ntarg=3, speedbound = 2, Tmax = 5,tequil = 5.,dtmin =0.001, start="influencer", maximize="follower",alpha=0.05,mindist=5.)
    r = (; ntarg, speedbound, Tmax,tequil,dtmin,start,maximize,alpha,mindist)
    #options for start: "zero", "mean", "influencer"
    #options for maximize:  "follower", "counter_follower"
    return r
end

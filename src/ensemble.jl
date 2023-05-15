
function sumgaussian(x, centers;sigma=0.1)
    output = 0
    for i in 1:size(centers,1)
        output = output +  gaussian(x, centers[i,:],sigma=sigma) 
    end
    return output
end 


function ABMsolveensemble(NT=200, N=10; savepoints = 4, init="4inf", p = ABMconstruct(), q= parameters(),sigma=0.1,save=true)
    (; X, Y, domain, dx, dy) = p
    (; n, L) = q
    x_arr = domain[1,1]:dx:domain[1,2]
    y_arr = domain[2,1]:dy:domain[2,2] 
    N_x = size(x_arr,1)
    N_y = size(y_arr,1)
    zs = zeros(2, 2, savepoints, N)
    ys = zeros(2, L, savepoints, N)
    us = zeros(N_x, N_y, 2, L, savepoints, N)
    savetimes = Int.(round.(LinRange(1, NT, savepoints)))
    states = [-1 1]

    for k=1:N
        xs, xinfs, infs, meds, state, _ = ABMsolve(NT;  p=p, q=q, init=init,chosenseed=k)

        for m in 1:savepoints
            t = savetimes[m]
            x = xs[t]
            xinf = xinfs[t]
            inf = infs[t]
            media = meds[t]
            for i in 1:2
                for j in 1:L
                    xi  = findall(x-> x==j, xinf)
                    xm = findall(x-> x==states[i], state)
                    choice = intersect(xi, xm)
                    us[:,:,i,j,m,k] = [(1/n)*sumgaussian([X[i,j], Y[i,j]], x[choice,:],sigma=sigma) for i in 1:size(X,1), j in 1:size(X,2)]
                    ys[:,j,m,k] = inf[j,:]
                end
                zs[:,i,m,k] = media[i,:]
            end
        end

    end
    if save==true
        us2 = us;
        zs2 = zs;
        ys2 = ys;
        p2 = p;
        q2 = q;
        @save string("data/abm_ensemble_",init,".jld2") us2 zs2 ys2 p2 q2
    end
    return us, zs, ys, (p,q)
end


function PDEsolveensemble(tmax=0.1, N=10; savepoints = 4, alg=nothing, init="4inf", p = PDEconstruct(), q= parameters(),save=true)
    (; N_x, N_y) = p
    L=q.L

    zs = zeros(2, 2, savepoints, N)
    ys = zeros(2, L, savepoints, N)
    us = zeros(N_x, N_y, 2, L, savepoints,  N)
    savetimes = LinRange(0, tmax, savepoints)

    #Threads.@threads
    local P 
    for i=1:N
        sol, P = PDEsolve(tmax; alg=alg, init=init, p=p, q=q) 

        for j in 1:savepoints
            u,z,y = sol2uyz(sol, savetimes[j])
            us[:,:,:,:,j,i] = u
            zs[:,:,j,i] = z
            ys[:,:,j,i] = y
        end

    end
    if save==true
        (p,q)=P
        @save string("data/pde_ensemble_",init,".jld2") us zs ys p q
    end
    return us, zs, ys, P
end


function runensembles(N; NT=200, tmax=2., savepoints = 5, q=parameters(),sigma=0.1,save=true)
    # pde
    us, zs, ys, (p,q)= PDEsolveensemble(tmax, 1;savepoints=savepoints, q=q,save=save)
    # abm
    us2, zs2, ys2, (p2,q2) = ABMsolveensemble(NT,N; savepoints=savepoints, q=q,sigma=sigma,save=save)
    plotensemblesnapshots(us, zs, ys, (p,q), us2, zs2, ys2, (p2,q), tmax; name="4inf",save=save)
    return us, zs, ys, (p,q), us2, zs2, ys2, (p2,q)
end


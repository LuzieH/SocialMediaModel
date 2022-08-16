
P = (; PDEconstruct()..., parameters()...,ABMconstruct()...)
 
(; grid_points, N_x, N_y,  a, b, c, sigma, eta, K_matrix, W_matrix, dx,dy, dV, C,  M , N, J, frictionM, frictionI) = P
(; dt, n, n_media, J, sigma, sigmahat, sigmatilde, a, b,c, frictionI, frictionM, eta) = P

x = [vec(X) vec(Y)]
n = size(x,1)
Net = ones(n,n)
state = (rand(n).>0.5)*2 .-1

xI = Any[]
push!(xI,  intersect(findall(x-> x>0,x[:,1]), findall(x-> x>0,x[:,2])))
push!(xI,  intersect(findall(x-> x<=0,x[:,1]), findall(x-> x>0,x[:,2])))
push!(xI,  intersect(findall(x-> x>0,x[:,1]), findall(x-> x<=0,x[:,2])) )
push!(xI,  intersect(findall(x-> x<=0,x[:,1]), findall(x-> x<=0,x[:,2]))) 

#follower network
fol=zeros(n, J)
for i in 1:J
    fol[xI[i],i] .=1
end

media =[-1. -1.; 1. 1.]

#initial opinions of influencers
inf = zeros(J,2)
for i in 1:J
    inf[i,:] = sum(x[xI[i],:],dims=1)/size(xI[i],1)
end

xinf = fol * collect(1:J)
###############

uABM = zeros(N_x, N_y, 2, J)
yABM = zeros(2,J)
zABM = zeros(2,2)

states=[-1 1]
for i in 1:2
    for j in 1:J
        xi  = findall(x-> x==j, xinf)
        xm = findall(x-> x==states[i], state)
        choice = intersect(xi, xm)
        uABM[:,:,i,j] = [sumgaussian([X[i,j], Y[i,j]], x[choice,:]) for i in 1:size(X,1), j in 1:size(X,2)]
        yABM[:,j] = inf[j,:]
    end
    zABM[:,i] = media[i,:]
end
uABM= (1/n)*uABM

AABM1 = a * attraction(x,Net)   
force1 =zeros(size(x))
force2 =zeros(size(x))
for j in 1:n
    if state[j]==1
        force1[j,:] = media[2,:]-x[j,:]
    else  
        force1[j,:] = media[1,:]-x[j,:]
    end

    for k in 1:J
        if fol[j,k]==1
            force2[j,:]=inf[k,:] -x[j,:]
        end
    end
end
force = c*force1 + b*force2
BCABM1 = force 

AABM = zeros(81,81,2)
AABM[:,:,1] = reshape(AABM1[:,1],(81,81,1))
AABM[:,:,2] = reshape(AABM1[:,2],(81,81,1))
BCABM= zeros(81,81,2)
BCABM[:,:,1] = reshape(BCABM1[:,1],(81,81,1))
BCABM[:,:,2] = reshape(BCABM1[:,2],(81,81,1))

BPDE = zeros(81,81,2,2,J)
CPDE = zeros(81,81,2,2,J)
for i in 1:2
    for j in 1:J
        BPDE[:,:,:,i,j] = b * follower_force(yABM[:,j], grid_points, N_x, N_y)   
        CPDE[:,:,:,i,j] = c * follower_force(zABM[:,i], grid_points, N_x, N_y)
    end
end

rhosum = sum(uABM, dims=(3,4))[:,:,1,1]
APDE = a* agent_force(rhosum, K_matrix, W_matrix, dV)

BCPDE = BPDE +CPDE

BCPDErep = zeros(81*81,2)
for i in 1:2
    for j in 1:J
        choice = intersect(findall(x->x==states[i], state), findall(x->x==1, fol[:,j]))
        BCPDErep[choice,:] = reshape(BCPDE[:,:,:,i,j], (81*81,2,1,1))[choice,:,1,1]
    end
end

BCPDErep = reshape(BCPDErep, (81,81,2))


heatmap(APDE[:,:,1]-AABM[:,:,1])
heatmap(BCABM[:,:,1]-BCPDErep[:,:,1])

# is the problem in comparin the A force due to not integrating over fine enough grid??
# but both times we are doing a discrete integration...
# something else... test on small example where difference lies!

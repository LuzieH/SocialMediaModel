#todo add to one function, global and then local!



function counterglobal(;mtime = 800000, meval=-1, alg=:LN_COBYLA, ntarg = 3,countercontrol = "med", p = PDEconstructmeso())
    return solveopt(; p = p, q= parameters_control(),r=parameters_optcont(ntarg=ntarg,maximize="counter_follower",start="zero"), alg=alg, mtime = mtime, meval=meval, ftol_rel=1e-4, xtol_rel = 1e-2, multistart =true,countercontrol=countercontrol,stubborntarget=[1.5 1.5],x0 = rand(Uniform(-2,2),2*ntarg))
end

function counterlocal(x0; mtime = 800000, meval=-1, alg=:LN_COBYLA, ntarg = 3, countercontrol = "med", p = PDEconstructmeso())
    return solveopt(;p=p, q= parameters_control(),r=parameters_optcont(ntarg=ntarg,maximize="counter_follower",start="zero"), alg=alg, mtime = mtime, meval=meval, ftol_rel=0., xtol_rel = 0., multistart =false,countercontrol=countercontrol,stubborntarget=[1.5 1.5],x0 = x0)
end


######################### FINAL EXPERIMENTS, started 23.12 at 13:00
# global 10 days = 864000, local 3 days=259200

function GLbob4med(;mtime = 864000, mtime2 = 259200)
    maxf,maxx, ret, numevals, x_list,followersum_list,cornersum_list, penalty_list = Gbob4med(;mtime = mtime)
    maxf2,maxx2, ret2, numevals2, x_list2,followersum_list2,cornersum_list2, penalty_list2 = Lbob4med(maxx; mtime =mtime2)
    return maxf,maxx, ret, numevals, x_list,followersum_list,cornersum_list, penalty_list, maxf2, maxx2, ret2, numevals2, x_list2,followersum_list2,cornersum_list2, penalty_list2
end
#1
#= 8-element Vector{Float64}:
  1.7569358194562925
 -1.0229031117738812
  1.5837631323150598
 -1.7838660406080011
  1.4297214501566768
 -1.831693586844958
  1.3636359489335177
 -1.7497664002376514
 =#


function GLbob4inf(;mtime = 864000, mtime2 = 259200)
    maxf,maxx, ret, numevals, x_list,followersum_list,cornersum_list, penalty_list = Gbob4inf(;mtime = mtime)
    maxf2,maxx2, ret2, numevals2, x_list2,followersum_list2,cornersum_list2, penalty_list2 = Lbob4inf(maxx; mtime =mtime2)
    return maxf,maxx, ret, numevals, x_list,followersum_list,cornersum_list, penalty_list, maxf2, maxx2, ret2, numevals2, x_list2,followersum_list2,cornersum_list2, penalty_list2
end

#1
#= 8-element Vector{Float64}:
 0.5497464162137617
 0.5314500291521418
 0.6703610033055708
 0.6378306837473999
 0.751676050022619
 0.7333139224985369
 0.7993856203742744
 0.7614622845878554 
 =#

function GLcob4med(;mtime = 864000, mtime2 = 259200)
    maxf,maxx, ret, numevals, x_list,followersum_list,cornersum_list, penalty_list = Gcob4med(;mtime = mtime)
    maxf2,maxx2, ret2, numevals2, x_list2,followersum_list2,cornersum_list2, penalty_list2 = Lcob4med(maxx; mtime =mtime2)
    return maxf,maxx, ret, numevals, x_list,followersum_list,cornersum_list, penalty_list, maxf2, maxx2, ret2, numevals2, x_list2,followersum_list2,cornersum_list2, penalty_list2
end
#1
#= 8-element Vector{Float64}:
  1.8106679226546316
 -0.9805088659711532
  1.7642225907551847
 -1.7710551914852335
  1.6338751208236737
 -1.8358368944346843
  1.5683845474974514
 -1.7962804090565507
  =#
function GLcob4inf(;mtime = 864000, mtime2 = 259200)
    maxf,maxx, ret, numevals, x_list,followersum_list,cornersum_list, penalty_list = Gcob4inf(;mtime = mtime)
    maxf2,maxx2, ret2, numevals2, x_list2,followersum_list2,cornersum_list2, penalty_list2 = Lcob4inf(maxx; mtime =mtime2)
    return maxf,maxx, ret, numevals, x_list,followersum_list,cornersum_list, penalty_list, maxf2, maxx2, ret2, numevals2, x_list2,followersum_list2,cornersum_list2, penalty_list2
end

#1
#= 8-element Vector{Float64}:
 0.5604084588036462
 0.5455196695338128
 0.642358069258012
 0.6527252366140318
 0.7454806015269124
 0.734952127347543
 0.8060867535176255
 0.7718781212141405 =#
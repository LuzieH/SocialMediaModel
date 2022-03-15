using Revise
includet("abm_pde.jl")

function transport(p, F, h=1)
  N = size(p, 1)
  D0 = 1/(2*h) * Tridiagonal(-ones(N-1), zeros(N), ones(N-1))

  DP = copy(D0)  # neumann boundary for dP
  DP[1,:] .= 0 
  DP[end,:] .= 0

  DF = copy(D0) # one-sided finite difference for dF
  DF[1,1:2] .= 1/h * [-1, 1]  
  DF[end, end-1:end] .= 1/h * [-1, 1]

  dx =  DP * p    .* F[:,:,1] + p .*  DF * F[:,:,1]
  dy = (DP * p')' .* F[:,:,2] + p .* (DF * F[:,:,2]')'

  div = dx + dy

  return div
end


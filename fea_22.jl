using LinearAlgebra
using Distributions
using NBInclude
using Plots
using NLsolve

@nbinclude("tauchen.ipynb")
#####################
## função recursos ##
#####################
function recursos(k;grid_z, α = 1/3,δ = 0.012)
    k.^(α)*grid_z' .+ (1-δ).*k 
end

Z,P = disc_t(7)
S = exp.(Z)

@time c_policy, k_policy = fea_galerkin(L,K,S,P)

plot(c_policy)
plot(k_policy)

plot(erros)
########################
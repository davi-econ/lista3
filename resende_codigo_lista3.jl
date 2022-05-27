using Distributions
using NBInclude
using Plots
using NLsolve
using LinearAlgebra
# funções a serem usadas
@nbinclude("tauchen.ipynb")

# parametros do modelo
β = beta = 0.987
μ = mu = 2
α = alpha = 1/3
δ = delta = 0.012
ρ = rho = 0.95
σ = sigma = 0.007

# capital de estado estacionario
kee = (α/(1/β -1 + δ))^(1/(1-α))
K = reshape(LinRange(0.75*kee,1.25*kee,500),1,500) |> Matrix
k_norm
# grid de choques
Z, P = disc_t(7)
S = exp.(Z)
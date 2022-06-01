using LinearAlgebra
using Distributions
using NBInclude
using Plots
using NLsolve
using BenchmarkTools


@nbinclude("tauchen.ipynb") # para discretizar os choques
@nbinclude("poli_cheb.ipynb") # Funções para o polinomio de Chebyshev
@nbinclude("phi_consumo.ipynb") # Função de phi e consumo
@nbinclude("colocacao.ipynb") # Elementos finitos + Colocação
@nbinclude("galerkin.ipynb") # Elementos finitos + Galerkin

# parametros do modelo
β = beta = 0.987
μ = mu = 2
α = alpha = 1/3
δ = delta = 0.012
ρ = rho = 0.95 
σ = sigma = 0.007

# capital de estado estacionario
kee = (α/(1/β -1 + δ))^(1/(1-α))
K = LinRange(0.75*kee,1.25*kee,500)
L = LinRange(0.75*kee,1.25*kee,10)
# grid de choques
Z, P = disc_t(7)
S = exp.(Z)


############################
## polinomio de Chebyshev ##
############################
@btime c_pc, k_pc = policy_pc(5,K,S,P)
plot(K,c_pc)
plot(K,k_pc)

######### Questão 2 ###############
a_chute = Vector(LinRange(2,4,length(L)))
for i in 1:(length(S)-1)
   a_chute = hcat(a_chute,Vector(LinRange(2,4,length(L))))
end
######################
#### fea colocação ###
######################
g(A) = res_colocacao(A,L,L,S,P)
a_col = nlsolve(g,a_chute).zero
@benchmark c_col, k_col = fea_colocacao(L,K,S,P)
plot(K,c_col)
plot(K,k_col)
plot(eee_colocacao(a_col,L,K,S,P))

#####################
#### fea galerkin ###
#####################
g(A) = res_galerkin(A,L,K,S,P)
a_gal = nlsolve(g,a_chute).zero
@benchmark c_gal, k_gal =  fea_galerkin(L,K,S,P)
plot(K,c_gal)
plot(K,k_gal)
plot(eee_galerkin(a_gal,L,K,S,P))
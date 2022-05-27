using Distributions
using NBInclude
using LinearAlgebra
using NLsolve
@nbinclude("tauchen.ipynb")
β = beta = 0.987
μ = mu = 2
α = alpha = 1/3
δ = delta = 0.012
ρ = rho = 0.95
σ = sigma = 0.007

kee = (α/(1/β -1 + δ))^(1/(1-α))

Z, P = disc_t(7)
S = exp.(Z)
##############################################
## seguindo o algoritmo como feito em  nota ##
##############################################
function raizes_cheby(ordem)
    m = ordem + 1
    raizes = -cos.((2*LinRange(1,m,m) . -1)*pi/(2*m))
    return raizes
end


## polinomio de chebyshev não recursivo 
function poli_cheb(ordem;capital)
    d = ordem
    t = zeros(length(capital),d+1)
    for j in 0:d
        t[:,j + 1] = cos.(j*acos.(capital))
    end
    return t
end
poli_cheb(2;capital = k_norm)

function consumo(gamas;teis)
    c = teis*gamas
    return c
end


#################################
## tudo isso na função residuo ##
#################################
grid_z = S
d = 1

function res(gamas,d,grid_z,P)
    nz = length(grid_z)
    res = zeros(d+1,nz)

    r_1 = raizes_cheby(d) # raízes para calcular o polinomio
    k_0 = cb_ss(r_1) # capital em nível

    t = poli_cheb(d; capital = r_1) # termos do polinomio

    c_0 = consumo(gamas; teis = t) # consumo com polinomios

    k_1 = k_0.^(1/3)*grid_z' .+ (1 - 0.012).*k_0 .- c_0
    capital_1 = cb_zero(k_1)  # normalizamos p/ [-1,1]

    for estado in 1:nz
        teis_1 = poli_cheb(d; capital = capital_1[:,estado]) # polinomios novos
        c_1 = consumo(gamas;teis = teis_1)
        ulinha = c_1.^(-2)
        dcdk = (1/3)*k_1[:,estado].^(-2/3)*grid_z' .+ (1-0.012) 
        lde = ulinha.*dcdk
        for id in 1:d+1
            res[id,estado] = c_0[id,estado]^(-2) - 0.987*dot(P[estado,:],lde[id,:])
        end
    end

    return res
end

res(ones(2,7),1,S,P)
guess = ones(2,7)
s=1
while s <= 5
    g(gamas) = res(gamas,s,S,P)
    gamaotimo = nlsolve(g,guess).zero
    guess = vcat(gamaotimo, zeros(1,7))
    s = s+1
end
gamaotimo

res(gamaotimo,5,S,P)
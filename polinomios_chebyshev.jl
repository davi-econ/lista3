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
@time while s <= 5
    g(gamas) = res(gamas,s,S,P)
    gamaotimo = nlsolve(g,guess).zero
    guess = vcat(gamaotimo, zeros(1,7))
    s = s+1
end
gamaotimo

res(gamaotimo,5,S,P)
t = poli_cheb(5;capital = k_norm)
c_0 = consumo(gamaotimo; teis = t)
c_0
k_1 = K.^(1/3)*S' .+ (1 - 0.012).*K .- c_0
k_norm_1 = cb_zero(k_1)
gamaotimo
nz = 7
eee = zeros(500,7)
estado = 1

k_norm
for is in 1:ns
    for ik in 1:nk
        c = consumo(gamaotimo;poli_cheb(5;k_norm[ik]))
        consumo[ik,is] = S[is]*K[ik]^(α) + (1-δ)*K[ik] - K[indice[ik,is]]
        u_l = c^(-μ)
        dcdk = S[is]*(α)*(K[ik]^(α-1))+(1-δ)
        u_lk[ik,is] = u_l*dcdk
    end
end
for is in 1:ns
    for ik in 1:nk
        inversa = (β*dot(P[is,:],u_lk[indice[ik,is],:]))^(1/-μ)
        pe = inversa/consumo[ik,is]
        euler = 1 - pe
        eee[ik,is] = log10(abs(euler))
    end
end 


##########################
#### erros de Euler ######
function euler_error_pc(gamas,d,grid_k,grid_z,P)
    nz = length(grid_z)
    nk = length(grid_k)
    r_1 = cb_zero(grid_k)
    t = poli_cheb(d; capital = r_1) # termos do polinomio

    c_0 = consumo(gamas; teis = t) # consumo com polinomios

    k_1 = grid_k.^(1/3)*grid_z' .+ (1 - 0.012).*grid_k .- c_0
    capital_1 = cb_zero(k_1)  # normalizamos p/ [-1,1]

    for estado in 1:nz
        teis_1 = poli_cheb(d; capital = capital_1[:,estado]) # polinomios novos
        c_1 = consumo(gamas;teis = teis_1)
        ulinha = c_1.^(-2)
        dcdk = (1/3)*k_1[:,estado].^(-2/3)*grid_z' .+ (1-0.012) 
        lde = ulinha.*dcdk
        for id in 1:nk
             inversa = (0.987*dot(P[estado,:],lde[id,:]))^(-1/2)
             eee[id,estado] = 1 - inversa/c_0[id,estado]
        end
    end
    eee = log10.(abs.(eee))
    return eee
end

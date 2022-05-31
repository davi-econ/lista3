using LinearAlgebra
using Distributions
using NBInclude
using Plots
using NLsolve

@nbinclude("tauchen.ipynb")
@nbinclude("galerkin.ipynb")
@nbinclude("colocacao.ipynb")

#########################
## Phi, o interpolante ##
#########################
function phi_inf(k;k_i,k_inf)
    (k - k_inf )/(k_i - k_inf)
end

function phi_sup(k;k_i,k_sup)
    (k_sup - k )/(k_sup - k_i)
end

## definir os intervalos de K ##
L = LinRange(0.75kee,1.25kee,10)
K = LinRange(0.75kee,1.25kee,500)

A = ones(10,7)

function phi_matricial(L,K)
    nl = length(L)
    nk = length(K)
    pee = zeros(nk,nl)
    for i in 1:nk
        for j in 2:(nl-1)
            if K[i] >= L[j-1] && K[i] < L[j]
                pee[i,j] = phi_inf(K[i]; k_i = L[j], k_inf = L[j-1])
            elseif K[i] >= L[j] && K[i] < L[j+1]
                pee[i,j] = phi_sup(K[i]; k_i = L[j], k_sup = L[j+1])
            else
                pee[i,j] = 0
            end
        end
        if K[i] >= L[1] && K[i] < L[2]
            pee[i,1] = phi_sup(K[i]; k_i = L[1], k_sup = L[2])
        elseif K[i] < L[1]
            pee[i,1] = 1
        else
            pee[i,1] = 0
        end
        if K[i] <= L[nl] && K[i] > L[nl-1]
            pee[i,nl] = phi_inf(K[i]; k_i = L[nl], k_inf = L[nl-1])
        elseif K[i] > L[nl]
            pee[i,nl] = 1
        else
            pee[i,nl] = 0
        end
    end
    return pee
end
pe = phi_matricial(L,K)

####################
## função consumo ##
####################
function fea_consumo(matriz_a,matriz_phi)
    matriz_phi*matriz_a
end
fea_consumo(A,pe)

####################
## função resíduo ##
#####################
k_1
function res_colocacao(A,L,grid_k,grid_z,P)
    nz = length(grid_z)
    nl = length(L)
    nk = length(grid_k)
    res = zeros(nk,nz)
    k_prov = zeros(nk)

    pe =  phi_matricial(L,grid_k) # termos da colocação
    c_0 = fea_consumo(A,pe) # consumo
    k_1 = grid_k.^(1/3)*grid_z' .+ (1 - 0.012).*grid_k .- c_0

    #k_prov = sum(k_1,dims = 2)/nz
    #phi1 = phi_matricial(L,k_prov)

    for estado in 1:nz
        phi1 = phi_matricial(L,k_1[:,estado])
        c_1 = fea_consumo(A[:,estado],phi1)
        ulinha = c_1.^(-2)
        dcdk = (1/3)*k_1[:,estado].^(-2/3)*grid_z' .+ (1-0.012) 
        lde = ulinha.*dcdk
        for id in 1:nl
            res[id,estado] = c_0[id,estado]^(-2) - 0.987*dot(P[estado,:],lde[id,:])
        end
    end
    return res
end

#####################
## função recursos ##
#####################
function recursos(k;grid_z, α = 1/3,δ = 0.012)
    k.^(α)*grid_z' .+ (1-δ).*k 
end

function fea_colocacao(L,grid_k,grid_z,P)
    A = Vector(LinRange(2,4,length(L)))
    for i in 1:(length(grid_z)-1)
       A= hcat(A,Vector(LinRange(2,4,10)))
    end
    g(A) = res(A,L,L,grid_z,P)
    a_otimo = nlsolve(g,A).zero

    c_0 = fea_consumo(a_otimo, phi_matricial(L,grid_k))
    k_1 = recursos(grid_k;grid_z = grid_z) .- c_0
    
    return c_0, k_1
end
@time c_policy, k_policy = fea_colloc(L,K,S,P)

plot(c_policy)
plot(k_policy)

##########################
#### erros de Euler ######

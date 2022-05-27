



function res(gamas,d)
    res = zeros(d+1)
    
    r_1 = raizes_cheby(d) # raízes para calcular o polinomio
    k_0 = cb_ss(r_1) # capital em nível

    teis = poli_cheb(d; capital = r_1) # termos do polinomio
    
    c_0 = consumo(gamas;teis,ordem = d) # consumo com polinomios


    k_1 = k_0.^(1/3) + (1-0.012).*k_0 - c_0
    capital_1 = cb_zero(k_1)  # normalizamos p/ [-1,1]

    teis_1 = poli_cheb(d; capital = capital_1) # polinomios novos
    c_1 = consumo(gamas;teis = teis_1,ordem = d) # conumo de amanha

    res = 0.987*(c_1).^(-2).*((1/3)*k_1.^(1/3 -1) .+ (1 - 0.012)) .- c_0.^(-2)

    return res
end


gamainicial = 1.0
for d in 1:6
    g(gamas) = res(gamas,d)
    gamainicial = nlsolve(g, hcat(gamainicial, 1.0)).zero
end


## teste
gamainicial
k_norm = cb_zero(K)
t = poli_cheb(5;capital = k_norm)
c_0 = consumo(gamainicial;teis = t,ordem = 5)
k_1 = K.^(1/3) + (1-0.012).*K - reshape(c_0,1,500)
cb_zero(k_1)
nlsolve(g, hcat(gamainicial, 0.0))
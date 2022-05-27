# Funçõe auxiliares, como normalização

## normalização de um vetor interio p/ o intervalo [a,b]
function change_bounds(x;a,b)
    cb = (b-a).*(x .- minimum(x))./(maximum(x) - minimum(x)) .+a
    return cb
end

## normalização de um ponto
function cb_unit(x;a,b,x_min,x_max)
    cb = (b-a)*(x - x_min)/(x_max - x_min) +a
    return cb
end


#### funções especificas para mudar o padrão de um grid para outro intervalo, mas centrado no estado estacionario ou em 0
# em torno do estado estacionario
function cb_ss(x)
    cb = (0.5*kee).*(x .+ 1)./(2) .+ 0.75*kee
    #cb = reshape(cb,1,length(cb)) |> Matrix
    return cb
end
# em torno de zero
function cb_zero(x)
    cb = (2).*(x .- 0.75*kee)./((0.5*kee)) .-1 
    #cb = reshape(cb,1,length(cb)) |> Matrix
    return cb
end


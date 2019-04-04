
# functions for host and parasites dynamic
struct Para
    by::Float64
    ey::Float64
    uy::Float64
    #βy::Float64 #create separate function with βy
end

function parasites(xx::Para, X0::Float64, Y0::Float64 ; timesteps::Int64, iter::Int64)
    densities = zeros(timesteps+1, 2)
    densities[1,:] = [X0, Y0]
    for i in 1:timesteps
        X = densities[i, 1]
        Y = densities[i, 2]

        for j in 1:iter
            dXdt = X *(xx.bx*(1 - X - Y) - xx.ux - xx.c*(xx.βy*Y))
            dYdt = Y *(xx.by * (1 - X - Y) - xx.uy + xx.c*xx.βy*X)
            X = X + dXdt
            Y = Y + dYdt
        end
        if all([X, Y] .> 0)
            densities[i+1,:] = [X, Y]
        else
            break
        end
    end
    return densities
end

using Plots

#paramètres
p = Para(1.0, 0.5, 0.1, 0.2, 0.4, 0.5, 0.5)
#modèle
N = parasites(p, 0.6, 0.1, timesteps = 50, iter = 1)
#représentation graphique
plot(N, xaxis="Time", yaxis="Hosts", title="Figure 3 a)")

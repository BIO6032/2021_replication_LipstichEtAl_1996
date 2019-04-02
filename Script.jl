
# functions for host and parasites dynamic
struct Para
    bx::Float64
    b1::Float64
    b2::Float64
    ux::Float64
    u1::Float64
    u2::Float64
    β1::Float64
    β2::Float64
    c::Float64
end

function parasites(xx::Para, X0::Float64, Y10::Float64, Y20::Float64 ; timesteps::Int64, iter::Int64)
    densities = zeros(timesteps+1, 3)
    densities[1,:] = [X0, Y10, Y20]
    for i in 1:timesteps
        X = densities[i, 1]
        Y1 = densities[i, 2]
        Y2 = densities[i, 3]

        for j in 1:iter
            dXdt = X *(xx.bx*(1 - X - Y1 - Y2) - xx.ux - xx.c*(xx.β1*Y1 + xx.β2*Y2))
            dY1dt = Y1 *(xx.b1 * (1 - X - Y1 - Y2) - xx.u1 + xx.c*xx.β1*X)
            dY2dt = Y2 *(xx.b2 *(1 - X - Y1 - Y2) - xx.u2 + xx.c*xx.β2*X)
            X = X + dXdt
            Y1 = Y1 + dY1dt
            Y2 = Y2 + dY2dt
        end
        if all([X, Y1, Y2] .> 0)
            densities[i+1,:] = [X, Y1, Y2]
        else
            break
        end
    end
    return densities
end

#Vertical case

struct VCparameters
    bi::Float64
    bx::Float64
    ux::Float64
    ui::Float64
end

function Vcase(VC::VCparameters, N::Float64; timesteps::Int64, iter::Int64)
    densities = zeros(timesteps+1, 1)
    densities[1,:] = [N]
    for i in 1:timesteps
        V = densities[i, 1]

        for j in 1:iter
            V = (VC.bi * VC.ux) / (VC.bx * VC.ui)
        end
        if all([V] .>= 0)
            densities[i+1,:] = [V]
        else
            break
        end
    end
    return densities
end


using Plots

#paramètres
p = Para(0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.0, 0.0,4)
#modèle
N = parasites(p, 0.6, 0.1, 0.3, timesteps = 50, iter = 1)
#représentation graphique
plot(N, xaxis="Time", yaxis="Hosts", title="Figure 3 a)")
fon(p, 100.0, 100.0, timesteps = 50, iter = 100)

#Vertical cases
pa = VCparameters(0.1 ,1.0 ,0.2 ,0.2)
N2 = Vcase(pa ,0.1 ,timesteps = 50, iter = 10)
plot(N2, xaxis="Time", yaxis="Vertical cases", title="")

using Plots
##fig3 e/f
timesteps = 1000

#Plotting by vs t
bx = 1.0

function by(bx)
    vertical = zeros(timesteps, 1)
    for t in 1:timesteps
        r1 = rand(Float64)
        r2 = rand(Float64)
        r3 = rand(Float64)

        by = bx*r3*(1-r1*r2)
        vertical[t, 1] = by
    end
    return vertical
end

# call the function by
by_results = by(bx)
show(by_results)
plot(by_results, xaxis = "Time", yaxis = "Vertical cases", title = "Fraction of all new cases acquired vertically")

#Figure 3, a/b
using Plots

struct parameters
    e_i::Float64
    b_i::Float64
    beta_i::Float64
    timesteps::Float64
    iter::Float64
end

r1 = rand(Float64)
r2 = rand(Float64)
r3 = rand(Float64)
ui = rand(Float64)

#constants
c = 0.5 #and test at c = 4
ux = 0.2
bx = 1.0
α = 0.2     #cost of vertical transmission, relative to horizontal. NEED TO FIND
K = 1       #carrying capacity NEED TO FIND

X = 0.1
Y1 = 0.2
Y2 = 0.15

function b_i(bx, r3, r1, r2)
    return bx*r3*(1-r1*r2)
end
bi <- b_i(bx, r3, r1, r2)

function e_i(bx, r3, r1, r2)
    return bx*(1-r3)*(1-r1*r3)
end
ei <- e_i(bx, r3, r1, r2)

function beta_i(r1, α, bi, bx)
    return r1- α*bi/bx
end
βi <- beta_i(r1, α, bi, bx)

#calculating number of hosts using equations 1,2,3
function hosts(bx, bi, ei, X, Y1, Y2, K, ux, ui, c, βi, timesteps = 50, iter = 100)
    densities = zeros(timesteps, 3)
    densities[1,:] = [X, Y1, Y2]
    for t in 1:timesteps
        X = densities [t, 1]
        Y1 = densities [t, 2]
        Y2 = densities [t, 3]
        for i in 1:iter
        dXdt = bx*X + ei*Y1 + ei*Y2)*(1 - (X + Y1 +Y2)/K) - ux*X - c*(βi*Y1 + βi*Y2)*X
        dY1dt = bi*Y1*(1 - (X + Y1 + Y2)/K) - ui*Y1 + c*βi*X*Y1
        dY2dt = bi*Y2*(1-(X + Y1 + Y2)/K) - ui*Y2 + c*βi*X*Y2
        X = X + dXdt
        Y1 = Y1 + dY1dt
        Y2 = Y2 + dY2dt
    end
    if all([X, Y1, Y2] .> 0)
        densities[i+1,:] = [X, Y1, Y2]
    else
        break
    end
    return densities
end

plot(hosts, xaxis = "Time", yaxis = "Hosts", title = "Host abondance over time")

############## Plotting R0 mean vs. time ##################
                ### (FIGURE 4c/d)
# constants
c = 4.0
bx = 1.0
ux = 0.2
α = 0.5              ### TO FIND (unknown)
K = 1.0              ### TO FIND

ui = rand(Float64)
r1 = rand(Float64)
r2 = rand(Float64)
r3 = rand(0:r1)     #returns zeros. find right way to code r3 ε (0, r1)

# bi function
function b_i(bx, r3, r1, r2)
    return bx*r3*(1-r1*r2)
end
bi = b_i(bx, r3, r1, r2)

# βi function
function beta_i(r1, α, r2, r3, bx)
    return r1-α*bi/bx
end
βi = beta_i(r1, α, r2, r3, bx)

# H0 function
function H_0(c, βi, ui, K, ux, bx)
    return c*βi/ui*K*(1-ux/bx)
end
H0 = H_0(c, βi, ui, K, ux, bx)

# V0 function
function V_0(bi, ux, bx, ui)
    return bi*ux/(bx*ui)
end
V0 = V_0(bi, ux, bx, ui)

# calculating R0
timesteps = 10

function R0(H0, V0)
    invasion = zeros(timesteps, 1)
    for t in 1:timesteps

        ui = rand(Float64)
        r1 = rand(Float64)
        r2 = rand(Float64)
        r3 = rand(0:r1)

        R = H0 + V0
        invasion[t,1] = R
    end
    return invasion
end

using Plots
println("Next attempt")
R0_results = R0(H0, V0)
show(R0_results)
# plotting code output
plot(R0_results, xaxis = "time", yaxis = "R0", title = "Average R0 in the Population",
colour=:red)

############## Plotting R0 mean vs. time ##################
                ### (FIGURE 2c/d)
# constants
c = 4.0
bx = 1.0
ux = 0.2
α = 0.5              ### TO FIND (unknown)
K = 1.0              ### TO FIND

# bi function
function bi(bx::Float64, r3::Float64, r1::Float64, r2::Float64)
    return bx*r3*(1-r1*r2)
end

# βi function
function βi(r1::Float64, α::Float64, r2::Float64, r3::Float64, bx::Float64)
    return r1-α*bi/bx
end
#println(βi(r1,α,bi,bx))

# H0 function
function H0(c::Float64, r1::Float64, α::Float64, ui::Float64, K::Float64, ux::Float64, bx::Float64)
    return c*βi/ui*K*(1-ux/bx)
end

# V0 function
function V0(r1::Float64, r2::Float64, r3::Float64, ux::Float64, bx::Float64, ui::Float64)
    return bi*ux/(bx*ui)
end

# calculating R0
function R0(c, ui, K, ux, bx, timesteps, beta_i, b_i)
    densities = zeros(timesteps, 1)
    #densities[1,:] = [R]
    for t in 1:timesteps
        #R = densities[t, 1]
        #plot R0
        ui = rand(Float64)
        r1 = rand(Float64)
        r2 = rand(Float64)
        r3 = rand(Float64)
        R = H0(c, beta_i, ui, K, ux, bx) + V0(b_i, ux, bx, ui)
        densities[t,1] = R
        #if all([R] .> 0)
        #    densities[i+1,:] = [R]
        #else
        #    break
        #end
    end
    return densities
end

using Plots
println("Next attempt")
beta_i = βi(r1,α,b_i,bx)
b_i = bi(bx,r3,r1,r2)
example = R0(c, beta_i, ui, K, ux, bx, b_i, 10)
show(example)
# plotting code output
plot(example, xaxis = "time", yaxis = "R0", title = "Average R0 in the Population", colour=:red)

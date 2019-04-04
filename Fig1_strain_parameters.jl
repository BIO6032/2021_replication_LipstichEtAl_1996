#include("Fig1_XY.jl")
using Plots

struct para_Y
    uy::Float64
    by::Float64
    ey::Float64
    βy::Float64
end

struct host_X
    ux::Float64
    bx::Float64
    ex::Float64
end

# function to create new Y strain (parasite machine)
function new_Y()
    ux=0.2
    uy=rand(Float64)
    by=rand(Float64)
    ey=rand(Float64)
    βy=3*(uy-ux)/(uy-ux+1)
    new_guy = para_Y(uy, by, ey, βy)
    return new_guy
end
#newY = new_Y()             #calling the function once to test

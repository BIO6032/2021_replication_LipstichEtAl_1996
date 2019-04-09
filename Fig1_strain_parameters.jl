struct para_Y
    uy::Float64
    by::Float64
    ey::Float64
    βy::Float64
end

struct beta.y
    r1::Float64
    α::Float64
    by::Float64
    bx::Float64
end

struct host_X
    ux::Float64
    bx::Float64
    ex::Float64
end

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

function parasites(xx::para_Y, X0::Float64, Y0::para_Y ; timesteps::Int64, iter::Int64)
    # calling the new_Y function (for 1000 strains)
    for g in 1:timesteps
        # calling the new_Y function in the container 1
        densities = zeros(timesteps+1, 5)
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
        end
        if all([X, Y] .> 0)
            densities[i+1,:] = [X, Y]
        else
            break
        end
    end

    return densities
end

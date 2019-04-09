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

function new_X()
    ux = 0.2
    bx = 1.0
    ex = 50
    new_host = host_X(ux, bx, ex)
    return new_host
end

#function to generate new parasite
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

function parasites(xx::host_X, yy::para_Y, X0::Float64, Y0::Float64 ; timesteps::Int64, iter::Int64)
    bx = 0.1;
    cc = 0.05;
    densities = zeros(timesteps+1, 2)
    # calling the new_Y function in the container 1
    densities[1,:] = [X0, Y0]
    # calling the new_Y function (for 1000 strains)
    for i in 1:timesteps
        X = densities[i, 1]
        Y = densities[i, 2]
        for j in 1:iter
            dXdt = X * (xx.bx*(1 - X - Y) - xx.ux - cc*(yy.βy*Y)) * (1 / iter)
            dYdt = Y * (yy.by * (1 - X - Y) - yy.uy + cc*yy.βy*X) * (1 / iter)
            X = X + dXdt
            Y = Y + dYdt
            print(dXdt)
        end
        if (mod(i,1000) == 0)
            # call function adding 1000 random strains --> call new_Y() 1000x?
            Y = Y + 100 #something
        end
        if all([X, Y] .> 0)
            densities[i+1,:] = [X, Y]
        else
            break
        end
    end
    return densities
end

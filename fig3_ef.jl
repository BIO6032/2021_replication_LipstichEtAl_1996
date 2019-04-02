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

# Defines the functions used to process the model parameters
"""
Differential equations formula for finding densities of each pop at next step
:param indiv: vector of number of individuals of each group
    (first idx is uninfected hosts; 2:end to each strain)
:param range: window for each time step
:param parameters: model parameters
:returns: differential eqn to solve
"""
function comp_densities(indiv, range, parameters)
    x = indiv[1] # number of uninfected hosts
    y = indiv[2:end] # number of hosts infected with each strain
    regul = (1 .- (sum(indiv)) ./ (range.K))
    #dx = (range.bx * x + sum(range.ei .* y)) * regul - range.ux * x - range.c * sum(range.βy .* y) * x
    dx = (range.bx .* x .+ sum(range.ei .* y)) .* regul .- range.ux .* x .- range.c .* sum(range.βy .* y) .* x
    dy = range.bi .* y .* regul .- range.ui .* y .+ range.c .* range.βy .* x .* y
    return vcat(dx, dy)
end
"Runs the ODE simulation of a population-dynamical model given a set of parameters. This is done for each strain introduction"
function run_simulation()
    @progress "Simulation" for i in 2:length(Y)
        # solve ODE using initial conditions
        prob = ODEProblem(comp_densities, new_U, (windowstart,windowend), parameters)
        solution = solve(prob, saveat=windowstart:1.0:windowend)
        # set values to 0 if negative
        for t in eachindex(solution.t)
            pop = solution.u[t]
            for i in 1:n_parasites
                if (pop.<0)[i]
                    pop[i] = 0
                end
            end
            N[:,Int(solution.t[t]+1)] = pop
        end
        # setup conditions & new parasite for next loop
        global new_U = solution[end]
        new_y = findfirst(x -> x == 0.0, new_U)
        new_U[new_y] = 1.0
        # set limits for next loop
        global windowstart = windowend
        global windowend = windowstart + windowsize
    end
end
"""calculates the evenness of an array (shoutout to Pielou)
    :param n: array
    :returns: evenness values through time
"""
function calculate_evenness(n)
    np = filter(x -> x > eps(), n)
    p = np./sum(np)
    ev = length(p) == 1 ? 0.0 : -sum(p.*log.(p))*(1/log(length(p)))
    return ev
end

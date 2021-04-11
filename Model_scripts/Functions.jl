
function comp_densities(indiv, range, params)
    #=
    differential equations formula for finding densities of each pop at next time step
    indiv: vector of number of individuals of each group (first idx is uninfected hosts; 2:end to each strain)
    range: window for each time step
    params: model parameters
    
    returns: differential eqn to solve
    =#
    x = indiv[1]
    y = indiv[2:end]
    regul = (1-(sum(indiv))/(range.K))
    dx = (range.bx*x + sum(range.ei.*y))*regul - range.ux*x-range.c*sum(range.βy.*y)*x
    dy = range.bi .* y .* regul .- range.ui .* y .+ range.c .* range.βy .* x .* y
    return vcat(dx, dy)
end

function run_simulation()
    # each strain introduction (1000x)
    @progress "Simulation" for i in 2:length(Y)
        # initial conditions
        prob = ODEProblem(comp_densities, new_U, (debut,fin), parameters)
        solution = solve(prob, saveat=debut:1.0:fin)
        for t in eachindex(solution.t)
            pop = solution.u[t]
            for i in 1:100
                if (pop.<0)[i]
                    pop[i] = 0
                end
            end
            N[:,Int(solution.t[t]+1)] = pop
        end

        # set conditions & new parasite for next loop
        global new_U = solution[end]
        new_y = findfirst(x -> x == 0.0, new_U)
        new_U[new_y] = 1.0

        # set limits for next loop
        global debut = fin
        #global fin = Float64((i+1).* 1000.0)
        global fin = debut + duree
    end
end

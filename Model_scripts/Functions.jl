# Defines the functions used to process the model parameters

# Statistical functions ######################################################
"""
Differential eqn formula for finding densities of each pop at next step
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
    dx =
        (range.bx .* x .+ sum(range.ei .* y)) .* regul .- range.ux .* x .-
        range.c .* sum(range.βy .* y) .* x
    dy = range.bi .* y .* regul .- range.ui .* y .+ range.c .* range.βy .* x .* y
    return vcat(dx, dy)
end

"""
Runs the ODE simulation of a population-dynamical model given a set of
parameters. This is done for each strain introduction.
"""
function run_simulation()
    @showprogress for i in 2:length(Y)
        # solve ODE using initial conditions
        prob = ODEProblem(
            comp_densities,
            new_U,
            (windowstart, windowend),
            parameters
        )
        solution = solve(prob; saveat=windowstart:1.0:windowend)
        # set values to 0 if negative
        for t in eachindex(solution.t)
            pop = solution.u[t]
            for i in 1:n_parasites
                if (pop .< 0)[i]
                    pop[i] = 0
                end
            end
            N[:, Int(solution.t[t] + 1)] = pop
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

"""
Calculates the average mortality across the population
(weighted for noise reduction)
    :returns: vector of mortality values through time
"""
function calculate_mortality()
    return sum(Np[:, 2:end] .* ui'; dims=2) ./ sum(Np[:, 2:end]; dims=2);
end

"""
Calculates the average horizontal transmission (βy) across the population
(weighted for noise reduction)
    :returns: vector of βy values through time
"""
function calculate_horizontal_transmission()
    return sum(Np[:, 2:end] .* βy'; dims=2) ./ sum(Np[:, 2:end]; dims=2);
end

"""
Calculates the average the number of new horizontally-acquired cases from a
single infected host before the primary host dies (H0) across the population
(weighted for noise reduction)
    :returns: vector of H0 values through time
"""
function calculate_H0()
    k = 1
    return c * βy_w_avg ./ ui_w_avg .* k .* (1 - ux / bx);
end

"""
Calculates the vertical basic reproductive ratio (V0) across the population
(weighted for noise reduction)
    :returns: vector of H0 values through time
"""
function calculate_V0()
    bi_avg = sum(Np[:, 2:end] .* bi'; dims=2) ./ sum(Np[:, 2:end]; dims=2);
    return bi_avg .* ux ./ (bx * ui_w_avg)
end

"""
Calculates the virulence across the population (weighted for noise reduction)
    :returns: vector of virulence values through time
"""
function calculate_average_virulence()
    return sum(Np[:, 2:end] .* r1'; dims=2) ./ sum(Np[:, 2:end]; dims=2)
end

"""
Calculates the evenness of an Vector (shoutout to Pielou)
    :param n: vector
    :returns: evenness values through time
"""
function calculate_evenness(n)
    np = filter(x -> x > eps(), n)
    p = np ./ sum(np)
    ev = length(p) == 1 ? 0.0 : -sum(p .* log.(p)) * (1 / log(length(p)))
    return ev
end

# Plotting functions #########################################################
"""
Plots the numbers of individuals from each strain, comparing infected
and uninfected hosts as well as the total number of parasites.
    :param labels: vector of labels for each host
    :param plot_title: string of plot title
    :param png_path: path (as str) to write a PNG version of the plot
"""
function plot_population_numbers(
    labels::Array{String},
    plot_title::String,
    png_path::String
)
    # plot the number of infected hosts
    plot(
        Np[:, 2:end];
        c=:blue,
        lw=1.5,
        alpha=0.4,
        title=plot_title,
        xlabel="Time",
        ylabel="Number of individuals",
        label=labels,
        ylims=(0, 100),
    )
    # add the number of uninfected hosts
    plot!(
        Np[:, 1];
        c=:black,
        lw=1.5,
        label="Uninfected"
    )
    # add the total number of parasites
    plot!(
        sum(Np[:, 2:end]; dims=2);
        c=:red,
        label="Total parasites"
    )
    # save figure as a PNG
    png(png_path)
end

"""
Plots mortality data with regards to time.
    :param png_path: path (as str) to write a PNG version of the plot
"""
function plot_mortality(png_path::String)
    plot(
        ui_w_avg;
        c=:black,
        lw=1.5,
        title="Average mortality in the population",
        xlabel="Time",
        ylabel="Mean mortality (ui)",
        leg=false,
        ylims=(0, 1),
    )
    # save figure as a PNG
    png(png_path)
end

"""
Plots the average reproductive rate (R0) with regards to time.
    :param upper_y_limit: upper limit for the y-axis
    :param png_path: path (as str) to write a PNG version of the plot
"""
function plot_reproductive_rate(upper_y_limit::Int, png_path::String)
    plot(
        H0_w + V0_w;
        c=:black,
        lw=1.5,
        title="Average R0 in the population",
        xlabel="Time",
        ylabel="Mean R0",
        leg=false,
        ylims=(0, upper_y_limit),
    )
    # save figure as a PNG
    png(png_path)
end

"""
Plots the vertical basic reproductive ratio (V0) with regards to time.
    :param png_path: path (as str) to write a PNG version of the plot
"""
function plot_vertical_reproductive_ratio(png_path::String)
    plot(
        V0_w;
        c=:black,
        lw=1.5,
        title="Average vertical cases in the population",
        xlabel="Time",
        ylabel="Mean V0",
        leg=false,
        ylims=(0, 1),
    )
    # save figure as a PNG
    png(png_path)
end

"""
Plots average horizontal transmission and virulence with regards to time.
    :param png_path: path (as str) to write a PNG version of the plot
"""
function plot_horizontal_and_virulence(png_path::String)
    # plot the horizontal transmission
    plot(
        βy_data;
        c=:black,
        title="Virulence and beta",
        label="Beta",
        xlabel="Time",
        ylabel="Mean virulence & \n Mean Beta",
        ylims=(0, 1),
    )
    # add virulence to the plot
    plot!(
        vir_w_avg;
        c=:blue,
        lw=1.5,
        label="Virulence"
    )
    # save figure as a PNG
    png(png_path)
end

"""
Calculates population evenness through its matrix.
Plots evenness data with regards to time.
    :param png_path: path (as str) to write a PNG version of the plot
"""
function plot_evenness(png_path::String)
    # calculate the evenness through time
    evenness_data = mapslices(
        calculate_evenness,
        Np[:, 2:end];
        dims=2
    );

    # plots the evenness data
    plot(
        evenness_data,
        c=:black,
        lw=0.5,
        title="Evenness",
        xlabel="Time",
        ylabel="Relative abundance (log)",
        label=false,
        leg=false,
        ylims=(0,1)
    )
    # save figure as a PNG
    png(png_path)
end

########### Defines the parameters & plot for Fig 1 a,c,e,g,i ###########
bi = 0.1; # set the birth rate of infected individuals as a constant 0.1
# import required variables/modules common to all simulations
include("../Model_scripts/fig1params.jl")
include("../Model_scripts/Functions.jl")
# now that we have all the requirements, let's run the simulation
run_simulation();
Np = N';
########### Figure 1a ###########
# plot the number of infected/uninfected hosts
labels = vcat("Infected", ["" for i in 1:1:n_parasites]);
plot_population_numbers(
    hcat(labels...),
    "by=0.1\n \nNumber of infected and uninfected hosts",
    100,
    ncol(Np),
    "$figure_directory/graph_1a.png"
)
########### Figure 1g ###########
# calculate the weighted average mortality
ui_w_avg = calculate_mortality()
# plot the weighted average mortality rate through time
plot_mortality(1, ncol(Np), "$figure_directory/graph_1g.png")
########### Figure 1e ###########
# calculate the weighted average β
βy_w_avg = calculate_horizontal_transmission()
# calculate weighted H0
H0_w = calculate_H0()
# calculate weighted V0
V0_w = calculate_V0()
# plot the average weighted V0
plot_vertical_reproductive_ratio(1, ncol(Np), "$figure_directory/graph_1e.png")
########### Figure 1c ###########
# plot the average weighted R0
plot_reproductive_rate(6, ncol(Np), "$figure_directory/graph_1c.png")
########### Figure 1i ###########
# calculate & plot the population evenness through time
plot_evenness(1, ncol(Np), "$figure_directory/graph_1i.png")

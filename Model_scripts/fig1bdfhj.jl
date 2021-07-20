########### Defines the parameters & plot for Fig 1 b,d,f,h,j ###########
bi = 1.0; # set the birth rate of infected individuals as a constant 1.0
# import required variables/modules common to all simulations
include("../Model_scripts/fig1params.jl")
include("../Model_scripts/Functions.jl")
# now that we have all the requirements, let's run the simulation
run_simulation();
Np = N';
########### Figure 1b ###########
# plot the number of infected/uninfected hosts
labels = vcat("Infected", ["" for i in 1:1:n_parasites]);
plot_population_numbers(
    hcat(labels...),
    "by=1.0\n \nNumber of infected and uninfected hosts",
    "$figure_directory/graph_1b.png"
)
########### Figure 1h ###########
# calculate the weighted average mortality
ui_w_avg = calculate_mortality()
# plot the weighted average mortality rate through time
plot_mortality("$figure_directory/graph_1h.png")
########### Figure 1f ###########
# calculate the weighted average β
βy_w_avg = calculate_horizontal_transmission()
# calculate weighted H0
H0_w = calculate_H0()
# calculate weighted V0
V0_w = calculate_V0()
# plot the average weighted V0
plot_vertical_reproductive_ratio("$figure_directory/graph_1f.png")
########### Figure 1d ###########
# plot the average weighted R0
plot_reproductive_rate(6, "$figure_directory/graph_1d.png")
########### Figure 1j ###########
# calculate & plot the population evenness through time
plot_evenness("$figure_directory/graph_1j.png")

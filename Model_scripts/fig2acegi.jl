########### Defines the parameters & plot for Fig 2 a,c,e,g,i ###########
c = 0.5; # set the host contact rate as a constant 0.5
# import required variables/modules common to all simulations
include("../Model_scripts/fig2params.jl");
include("../Model_scripts/Functions.jl");
# now that we have all the requirements, let's run the simulation
run_simulation();
Np = N';
########### Figure 2a ###########
# plot the number of infected/uninfected hosts
labels = vcat("Infected", ["" for i in 1:1:n_parasites]);
plot_population_numbers(
    hcat(labels...),
    "c=0.5\n \nNumber of infected and uninfected hosts",
    "$figure_directory/graph_2a.png"
)
########### Figure 2e ###########
# calculate the weighted average mortality
ui_w_avg = calculate_mortality()
# calculate weighted V0
V0_w = calculate_V0()
# plot the average V0
plot_vertical_reproductive_ratio("$figure_directory/graph_2e.png")
########### Figure 2c ###########
# calculate the weighted average β
βy_w_avg = calculate_horizontal_transmission()
# calculate weighted H0
H0_w = calculate_H0()
# plot average weighted R0
plot_reproductive_rate(5, "$figure_directory/graph_2c.png")
########### Figure 2g ###########
# calculate the weighted virulence (for noise reduction)
vir_w_avg = calculate_average_virulence()
# plot the average horizontal transmission & weighted virulence
plot_horizontal_and_virulence("$figure_directory/graph_2g.png")
########### Figure 2i ###########
# calculate & plot the population evenness through time
plot_evenness("$figure_directory/graph_2i.png")

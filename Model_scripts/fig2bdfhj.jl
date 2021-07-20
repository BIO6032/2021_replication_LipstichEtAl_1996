########### Defines the parameters & plot for Fig 2 b,d,f,h,j ###########
c = 4.0 # set the host contact rate as a constant 4.0
# import required variables/modules common to all simulations
include("../Model_scripts/fig2params.jl");
include("../Model_scripts/Functions.jl");
# now that we have all the requirements, let's run the simulation
run_simulation();
Np = N'
########### Figure 2b ###########
# plot the number of infected/uninfected hosts
labels = vcat("Infected", ["" for i in 1:1:n_parasites]);
plot_population_numbers(
    hcat(labels...),
    "c=4.0\n \nNumber of infected and uninfected hosts",
    "$figure_directory/graph_2b.png"
)
########### Figure 2d ###########
# calculate the weighted average mortality
ui_w_avg = calculate_mortality()
# calculate the weighted average β
βy_w_avg = calculate_horizontal_transmission()
# calculate weighted H0
H0_w = calculate_H0()
# calculate weighted V0
V0_w = calculate_V0()
# plot the weighted R0
plot_reproductive_rate(5, "$figure_directory/graph_2d.png")
########### Figure 2f ###########
# plot the average V0
plot_vertical_reproductive_ratio("$figure_directory/graph_2f.png")
########### Figure 2h ###########
# calculate the weighted virulence (for noise reduction)
vir_w_avg = calculate_average_virulence()
# plot the average horizontal transmission & weighted virulence
plot_horizontal_and_virulence("$figure_directory/graph_2h.png")
########### Figure 2j ###########
# calculate & plot the population evenness through time
plot_evenness("$figure_directory/graph_2j.png")

########### Defines the parameters & plot for Fig 3 b,d,f,h,j ###########
c = 4.0; # set the host contact rate as a constant 4.0
# import required variables/modules common to all simulations
include("../Model_scripts/fig3params.jl");
include("../Model_scripts/Functions.jl");
# now that we have all the requirements, let's run the simulation
run_simulation();
Np = N';
########### Figure 3b ###########
# plot the number of infected/uninfected hosts
labels = vcat("Infected", ["" for i in 1:1:n_parasites]);
plot_population_numbers(
    hcat(labels...),
    "c=4.0\n \nNumber of infected and uninfected hosts",
    100,
    size(Np)[1],
    "$figure_directory/graph_3b.png"
)
########### Figure 3d ###########
# calculate the weighted average mortality
ui_w_avg = calculate_mortality()
# calculate the weighted average β
βy_w_avg = calculate_horizontal_transmission()
# calculated the weighted H0
H0_w = calculate_H0()
# calculate the weighted V0
V0_w = calculate_V0()
# plot the average weighted R0
plot_reproductive_rate(5, size(Np)[1], "$figure_directory/graph_3d.png")
########### Figure 3f ###########
# plot the average V0
plot_vertical_reproductive_ratio(1, size(Np)[1], "$figure_directory/graph_3f.png")
########### Figure 3h ###########
# calculate the average virulence
vir_w_avg = calculate_average_virulence()
# plot the average horizontal transmission & weighted virulence
plot_horizontal_and_virulence(1, size(Np)[1], "$figure_directory/graph_3h.png")
########### Figure 3j ###########
# calculate & plot the population evenness through time
plot_evenness(1, size(Np)[1], "$figure_directory/graph_3j.png")

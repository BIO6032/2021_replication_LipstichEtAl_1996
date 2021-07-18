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
    population_matrix=Np,
    labels=hcat(labels...),
    plot_title="by=1.0\n \nNumber of infected and uninfected hosts",
    png_path="$figure_directory/graph_1b.png"
)
########### Figure 1h ###########
# calculate the weighted average mortality
ui_w_avg = calculate_mortality(population_matrix=Np, ui=ui)
# plot the weighted average mortality rate through time
plot_mortality(data=ui_w_avg, png_path="$figure_directory/graph_1h.png")
########### Figure 1f ###########
# calculate the weighted average β
βy_w_avg = calculate_horizontal_transmission(population_matrix=Np, βy=βy)
# calculate weighted H0
H0_w = calculate_H0(c=c, βy_w_avg=βy_w_avg, ui_w_avg=ui_w_avg, ux=ux, bx=bx)
# calculate weighted V0
V0_w = calculate_V0(population_matrix=Np, bi=bi, bx=bx, ux=ux, ui_w_avg=ui_w_avg)
# plot the average weighted V0
plot_vertical_reproductive_ratio(
    data=V0_w,
    png_path="$figure_directory/graph_1f.png"
)
########### Figure 1d ###########
# plot the average weighted R0
plot_reproductive_rate(
    data=H0_w + V0_w,
    upper_y_limit=6,
    png_path="$figure_directory/graph_1d.png"
)
########### Figure 1j ###########
# calculate & plot the population evenness through time
plot_evenness(data=Np, png_path="$figure_directory/graph_1j.png")

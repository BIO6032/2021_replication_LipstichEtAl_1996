########### Defines the parameters & plot for Fig 3 a,c,e,g,i ###########
c = 0.5; # set the host contact rate as a constant 0.5
# import required variables/modules common to all simulations
include("../Model_scripts/fig3params.jl");
include("../Model_scripts/Functions.jl");
# now that we have all the requirements, let's run the simulation
run_simulation();
Np = N'
########### Figure 3a ###########
# plot the number of infected/uninfected hosts
labels = vcat("Infected", ["" for i in 1:1:n_parasites]);
plot_population_numbers(
    population_matrix=Np,
    labels=hcat(labels...),
    plot_title="c=0.5\n \nNumber of infected and uninfected hosts",
    png_path="Figure3/graph_3a.png"
)
########### Figure 3c ###########
# calculate the weighted average mortality
ui_w_avg = calculate_mortality(population_matrix=Np, ui=ui)
# calculate the weighted average β
βy_w_avg = calculate_horizontal_transmission(population_matrix=Np, βy=βy)
# calculate weighted H0
H0_w = calculate_H0(c=c, βy_w_avg=βy_w_avg, ui_w_avg=ui_w_avg, ux=ux, bx=bx)
# calculate weighted V0
V0_w = calculate_V0(population_matrix=Np, bi=bi, bx=bx, ux=ux, ui_w_avg=ui_w_avg)
# plot the average weighted R0
plot_reproductive_rate(
    data=H0_w + V0_w,
    upper_y_limit=5,
    png_path="Figure3/graph_3c.png"
)
########### Figure 3e ###########
# plot the average V0
plot_vertical_reproductive_ratio(data=V0_w, png_path="Figure3/graph_3e.png")
########### Figure 3g ###########
# calculate the weighted virulence (for noise reduction)
vir_w_avg = calculate_average_virulence(population_matrix=Np, r1=r1)
# plot the average horizontal transmission & weighted virulence
plot_horizontal_and_virulence(
    βy_data=βy_w_avg,
    vir_data=vir_w_avg,
    png_path="Figure3/graph_3g.png"
)
########### Figure 3i ###########
# calculate & plot the population evenness through time
plot_evenness(data=Np, png_path="Figure3/graph_3i.png")

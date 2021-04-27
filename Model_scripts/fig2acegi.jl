########### Defines the parameters & plot for Fig 2 a,c,e,g,i ###########

c = 0.5; # set the host contact rate as a constant 0.5

# import required variables/modules common to all simulations
include("../Model_scripts/fig2params.jl");
include("../Model_scripts/Functions.jl")

# now that we have all the requirements, let's run the simulation
run_simulation();

Np = N';

########### Figure 2a ###########

labels = ["" for i = 1:1:n_parasites];
labels2 = vcat("Infected", labels);
labels2 = hcat(labels2...);

# plot the number of infected hosts
plot(
    Np[:,2:end],
    c=:blue,
    lw=1.5,
    alpha=0.4,
    title="c=0.5\n \nNumber of infected and uninfected hosts",
    xlabel="Time",
    ylabel="Number of individuals",
    label=labels2,
    ylims=(0,100)
)

# add the number of uninfected hosts
plot!(
    Np[:,1],
    c=:black,
    lw=1.5,
    label="Uninfected"
)

# add the total number of parasites
plot!(
    sum(Np[:,2:end]; dims=2),
    c=:red,
    label="Total parasites")

# save figure as a PNG
png("Figure2/graph_2a.png")


########### Figure 2e ###########

# weighted average mortality (for noise reduction)
ui_w_avg = sum(Np[:,2:end] .* ui'; dims=2) ./ sum(Np[:,2:end]; dims=2);

# calculate weighted V0
bi_avg = sum(Np[:,2:end] .* bi'; dims=2) ./ sum(Np[:,2:end]; dims=2)
V0_w = bi_avg .* ux ./ (bx * ui_w_avg);

# plot the average V0
plot(
    V0_w,
    c=:black,
    lw=1.5,
    title="Average V0 in the population",
    xlabel="Time",
    ylabel="Mean V0",
    leg=false,
    ylims=(0.0,1.0)
)

# save figure as a PNG
png("Figure2/graph_2e.png")


########### Figure 2c ###########

# calculate weighted average β (for noise reduction)
βy_w_avg = sum(Np[:,2:end] .* βy'; dims=2) ./ sum(Np[:,2:end]; dims=2);

# calculate weighted H0
k = 1
H0_w = c * βy_w_avg ./ ui_w_avg .* k .* (1 - ux / bx);

# calculate weighted R0
R0_w = H0_w + V0_w;

# plot average weighted R0
plot(
    R0_w,
    c=:black,
    lw=1.5,
    title="Average R0 in the population",
    xlabel="Time",
    ylabel="Mean R0",
    leg=false,
    ylims=(0,5)
)

# save figure as a PNG
png("Figure2/graph_2c.png")


########### Figure 2g ###########

# weighted virulence (for noise reduction)
vir_w_avg = sum(Np[:,2:end] .* r1'; dims=2) ./ sum(Np[:,2:end]; dims=2);

# plot the average weighted horizontal transmission
plot(
    βy_w_avg,
    c=:black,
    title="Virulence and beta",
    label="Beta",
    xlabel="Time",
    ylabel="Mean virulence & \n Mean Beta",
    ylims=(0,1)
)

# add the average weighted virulence to the plot
plot!(
    vir_w_avg,
    c=:blue,
    lw=1.5,
    label="Virulence"
)

# save figure as a PNG
png("Figure2/graph_2g.png")


########### Figure 2i ###########

# calculate the evenness through time
evenness_data = mapslices(calculate_evenness, Np[:,2:end]; dims=2);

# plot the evenness through time
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
png("Figure2/graph_2i.png")

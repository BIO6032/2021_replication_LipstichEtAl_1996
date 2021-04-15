########### Defines the parameters & plot for Fig 1 b,d,f,h,j ###########

bi = 1.0; # set the birth rate of infected individuals as a constant 1.0

# import required variables/modules common to all simulations
include("../Model_scripts/fig1params.jl")
include("../Model_scripts/Functions.jl")

# now that we have all the requirements, let's run the simulation
run_simulation();

Np = N';

########### Figure 1b ###########

labels = ["" for i = 1:1:n_parasites];
labels2 = vcat("Infected", labels);
labels2 = hcat(labels2...);

# plot the number of infected hosts
plot(
    Np[:,2:end],
    c=:blue,
    lw=1.5,
    alpha=0.4,
    title = "by=1.0\n \nNumber of infected and uninfected hosts",
    xlabel = "Time",
    ylabel = "Number of individuals",
    label = labels2,
    ylims =(0,100)
)

# add the number of uninfected hosts
plot!(
    Np[:,1],
    c=:black,
    lw=1.5,
    label = "Uninfected"
)

# add the total number of parasites
plot!(
    sum(Np[:,2:end]; dims=2),
    c=:red,
    label = "Total parasites"
)

# save figure as a PNG
png("Figure1/graph_1b.png")


########### Figure 1h ###########

# strains that survive at each time step for ui
survival = (Np .> 0.0)[:,2:end];

# weighted average mortality (for noise reduction)
ui_w_avg = sum(Np[:,2:end] .* ui'; dims=2) ./ sum(Np[:,2:end]; dims=2);

# plot the weighted average mortality rate
plot(
    ui_w_avg,
    c=:black,
    lw=1.5,
    title="Average mortality in the population",
    xlabel="Time",
    ylabel="Mean mortality (ui)",
    leg=false,
    ylims=(0,1)
)

# save figure as a PNG
png("Figure1/graph_1h.png")


########### Figure 1f ###########

# strains that survive for βy
survived_βy = survival .* βy';
βy_avg = sum(survived_βy; dims=2) ./ sum(survival; dims=2);

# weighted (for noise reduction)
βy_w_avg = sum(Np[:,2:end] .* βy'; dims=2) ./ sum(Np[:,2:end]; dims=2);

# calculating weighted H0
k = 1;
H0_w = c * βy_w_avg ./ ui_avg .* k .* (1 - ux / bx);

# calculating weighted V0
bi_avg = sum(Np[:,2:end] .* bi'; dims=2) ./ sum(Np[:,2:end]; dims=2);
V0_w = bi_avg .* ux ./ (bx * ui_avg);

# plot the average weighted V0
plot(
    V0_w,
    c=:black,
    lw=1.5,
    title="Average vertical cases in the population",
    xlabel="Time",
    ylabel="Mean V0",
    leg=false,
    ylims=(0,1)
)

# save figure as a PNG
png("Figure1/graph_1f")


########### Figure 1d ###########

# calculating weighted R0
R0_w = H0_w + V0_w;

# plot the average weighted R0
plot(
    R0_w,
    c=:black,
    lw=1.5,
    title="Average R0 in the population",
    xlabel="Time",
    ylabel="Mean R0",
    leg=false,
    ylims=(0,7)
)

# save figure as a PNG
png("Figure1/graph_1d")


########### Figure 1j ###########

# use the function in Functions.jl to calculate the evenness
evenness_data = mapslices(calculate_evenness, Np[:,2:end]; dims=2);

# plot the evenness through time
plot(
    evenness_data,
    c=:black,
    lw=1.5,
    title="Evenness",
    xlabel="Time",
    ylabel="Relative abundance (log)",
    leg=false,
    ylims=(0,1)
)

# save figure as a PNG
png("Figure1/graph_1j.png")

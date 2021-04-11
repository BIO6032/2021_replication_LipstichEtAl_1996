########### Defines the parameters & plot for Fig 2 b,d,f,h,j ###########

c = 4.0 # set the host contact rate as a constant 4.0

# import required variables/modules common to all simulations
include("../Model_scripts/fig2params.jl")
include("../Model_scripts/Functions.jl")

# now that we have all the requirements, let's run the simulation
run_simulation();

Np = N'

########### Figure 2b ###########

# define labels
labels = ["" for i = 1:1:n_parasites];
labels2 = vcat("Infected", labels);
labels2 = hcat(labels2...);

# plot the number of infected hosts
plot(
    Np[:,2:end],
    c=:blue,
    lw=1.5,
    alpha=0.4,
    title="c=4.0\n \nNumber of infected and uninfected hosts",
    xlabel="Time",
    ylabel="Number of individuals",
    label=labels2,
    ylims=(0,70)
)

# add the number of uninfected hosts
plot!(
    Np[:,1],
    c=:black,
    lw=1.5,
    label="Uninfected"
)
#plot!(sum(Np[:,2:end]; dims=2), label = "Total # parasites")

# save figure as a PNG
png("Figure2/graph_2b.png")


########### Figure 2d ###########

# strains that survive at each time step for ui
survival = (Np.>0.0)[:,2:end];

# weighted average mortality (for noise reduction)
ui_avg = sum(Np[:,2:end].*ui'; dims=2)./sum(Np[:,2:end]; dims=2);

# strains that survive for βy
survived_βy = survival .* βy';
βy_avg = sum(survived_βy; dims=2) ./ sum(survival; dims=2);

# weighted β (for noise reduction)
βy_w_avg = sum(Np[:,2:end] .* βy'; dims=2) ./ sum(Np[:,2:end]; dims=2);

# calculating H0
k = 1;
#H0 = c * βy_avg ./ ui_avg .* k .* (1 - ux / bx);
H0_w = c * βy_w_avg ./ ui_avg .* k .* (1 - ux / bx);

# calculating V0
bi_avg = sum(Np[:,2:end].*bi'; dims=2) ./ sum(Np[:,2:end]; dims=2);
V0 = bi_avg .* ux ./ (bx * ui_avg);
V0_w = bi_avg .* ux ./ (bx * ui_avg);

# calculating R0
R0_w = H0_w + V0_w;
#R0 = H0 + V0;
#plot(R0, title = "Average R0 in the population", xlabel = "Time", ylabel = "Mean R0", leg = false)

# plot the weighted R0
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
png("Figure2/graph_2d")


########### Figure 2f ###########

# plot average V0 (TODO: non-weighted?)
plot(
    V0,
    c=:black,
    lw=1.5,
    title="Average V0 in the population",
    xlabel="Time",
    ylabel="Mean V0",
    leg=false,
    ylims=(0.0,1.0)
)

# save figure as a PNG
png("Figure2/graph_2f")


########### Figure 2h ###########

# mean virulence
vir_avg = (1 .- (bi .+ ei) .* ux ./ (bx .* ui));

#survived_vir = survival .* vir_avg';
#vir_avg = sum(survived_vir; dims=2) ./ sum(survival; dims=2);
#plot!(vir_avg, label = "virulence")

# weighted virulence (for noise reduction)
vir_w_avg = sum(Np[:,2:end] .* vir_avg'; dims=2) ./ sum(Np[:,2:end]; dims=2);

# plot the average horizontal transmission (TODO: NON-WEIGHTED???????????????)
plot(
    βy_avg,
    c=:black,
    title="Virulence and beta",
    label="Beta",
    xlabel="Time",
    ylabel="Mean virulence & \n Mean Beta",
    ylims=(0,1)
)
#plot(βy_w_avg, title = "Mean virulence and beta", label = "beta") <- weighted

# add the average weighted virulence to the plot
plot!(
    vir_w_avg,
    c=:blue,
    lw=1.5,
    label = "Virulence"
    )

# save figure as a PNG
png("Figure2/graph_2h.png")


########### Figure 2j ###########

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
    eg=false,
    ylims=(0,1)
)

# save figure as a PNG
png("Figure2/graph_2j.png")

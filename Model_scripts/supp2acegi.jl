# plot first 5000 time points
# plot the number of infected hosts
plot(
    Np[:, 2:end];
    c=:blue,
    lw=1.5,
    alpha=0.4,
    title="c=0.5\n \nNumber of infected and uninfected hosts",
    xlabel="Time",
    ylabel="Number of individuals",
    label=labels2,
    ylims=(0, 100),
    xlims=(0,5e3)
)
# add the number of uninfected hosts
plot!(Np[:, 1]; c=:black, lw=1.5, label="Uninfected")
# add the total number of parasites
plot!(sum(Np[:, 2:end]; dims=2); c=:red, label="Total parasites")
# save figure as a PNG
png("Figure2/supplemental/supp_2a.png")
# plot the average V0
plot(
    V0_w;
    c=:black,
    lw=1.5,
    title="Average V0 in the population",
    xlabel="Time",
    ylabel="Mean V0",
    leg=false,
    ylims=(0.0, 1.0),
    xlims=(0, 5e3)
)
# save figure as a PNG
png("Figure2/supplemental/supp_2e.png")
# plot average weighted R0
plot(
    R0_w;
    c=:black,
    lw=1.5,
    title="Average R0 in the population",
    xlabel="Time",
    ylabel="Mean R0",
    leg=false,
    ylims=(0, 5),
    xlims=(0, 5e3)
)
# save figure as a PNG
png("Figure2/supplemental/supp_2c.png")
# plot the average weighted horizontal transmission
plot(
    βy_w_avg;
    c=:black,
    title="Virulence and beta",
    label="Beta",
    xlabel="Time",
    ylabel="Mean virulence & \n Mean Beta",
    ylims=(0, 1),
    xlims=(0, 5e3)
)
# add the average weighted virulence to the plot
plot!(vir_w_avg; c=:blue, lw=1.5, label="Virulence")
# save figure as a PNG
png("Figure2/supplemental/supp_2g.png")
# plot the evenness through time
plot(
    evenness_data;
    c=:black,
    lw=0.5,
    title="Evenness",
    xlabel="Time",
    ylabel="Relative abundance (log)",
    label=false,
    leg=false,
    ylims=(0, 1),
    xlims=(0, 5e3)
)
# save figure as a PNG
png("Figure2/supplemental/supp_2i.png")

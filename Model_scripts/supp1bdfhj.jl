# plot first 5000 time points
# plot the number of infected hosts
plot(
    Np[:, 2:end];
    c=:blue,
    lw=1.5,
    alpha=0.4,
    title="by=1.0\n \nNumber of infected and uninfected hosts",
    xlabel="Time",
    ylabel="Number of individuals",
    label=labels2,
    ylims=(0, 100),
    xlims=(0, 5e3)
)
# add the number of uninfected hosts
plot!(Np[:, 1]; c=:black, lw=1.5, label="Uninfected")
# add the total number of parasites
plot!(sum(Np[:, 2:end]; dims=2); c=:red, label="Total parasites")
# save figure as a PNG
png("Figure1/supplemental/supp_1b.png")
# plot the weighted average mortality rate through time
plot(
    ui_w_avg;
    c=:black,
    lw=1.5,
    title="Average mortality in the population",
    xlabel="Time",
    ylabel="Mean mortality (ui)",
    leg=false,
    ylims=(0, 1),
    xlims=(0, 5e3)
)
# save figure as a PNG
png("Figure1/supplemental/supp_1h.png")
# plot the average weighted V0
plot(
    V0_w;
    c=:black,
    lw=1.5,
    title="Average vertical cases in the population",
    xlabel="Time",
    ylabel="Mean V0",
    leg=false,
    ylims=(0, 1),
    xlims=(0, 5e3)
)
# save figure as a PNG
png("Figure1/supplemental/supp_1f.png")
# plot the average weighted R0
plot(
    R0_w;
    c=:black,
    lw=1.5,
    title="Average R0 in the population",
    xlabel="Time",
    ylabel="Mean R0",
    leg=false,
    ylims=(0, 6),
    xlims=(0, 5e3)
)
# save figure as a PNG
png("Figure1/supplemental/supp_1d.png")
# plot the evenness through time
plot(
    evenness_data;
    c=:black,
    lw=0.5,
    title="Evenness",
    xlabel="Time",
    ylabel="Relative abundance (log)",
    leg=false,
    ylims=(0, 1),
    xlims=(0,5e3)
)
# save figure as a PNG
png("Figure1/supplemental/supp_1j.png")

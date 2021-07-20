# plot first 5000 time points (before stabilization)
include("../Model_scripts/Functions.jl");
# plot the number of infected & uninfected hosts
plot_population_numbers(
    hcat(labels...),
    "by=0.1\n \nNumber of infected and uninfected hosts",
    100,
    5e3,
    "$figure_directory/supplemental/supp_1a.png"
)
# plot the average weighted R0
plot_reproductive_rate(6, 5e3, "$figure_directory/supplemental/supp_1c.png")
# plot the average weighted V0
plot_vertical_reproductive_ratio(1, 5e3, "$figure_directory/supplemental/supp_1e.png")
# plot the weighted average mortality rate through time
plot_mortality(1, 5e3, "$figure_directory/supplemental/supp_1g.png")
# plot the evenness through time
plot_evenness(1, 5e3, "$figure_directory/supplemental/supp_1i.png")

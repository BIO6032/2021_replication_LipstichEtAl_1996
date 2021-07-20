# plot first 20,000 time steps (before stabilization)
include("../Model_scripts/Functions.jl");
# plot the number of infected & uninfected hosts
plot_population_numbers(
    hcat(labels...),
    "c=0.5\n \nNumber of infected and uninfected hosts",
    100,
    2e4,
    "$figure_directory/supplemental/supp_3a.png"
)
# plot the average weighted R0
plot_reproductive_rate(5, 2e4, "$figure_directory/supplemental/supp_3c.png")
# plot the average V0
plot_vertical_reproductive_ratio(1, 2e4, "$figure_directory/supplemental/supp_3e.png")
# plot the average horizontal transmission
plot_horizontal_and_virulence(1, 2e4, "$figure_directory/supplemental/supp_3g.png")
# plot the evenness through time
plot_evenness(1, 2e4, "$figure_directory/supplemental/supp_3i.png")

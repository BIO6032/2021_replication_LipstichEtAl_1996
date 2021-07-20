# plot first 6,000 time steps (before stabilization)
include("../Model_scripts/Functions.jl");
# plot the number of infected & uninfected hosts
plot_population_numbers(
    hcat(labels...),
    "c=4.0\n \nNumber of infected and uninfected hosts",
    5,
    6e3,
    "$figure_directory/supplemental/supp_3b.png"
)
# plot the average weighted R0
plot_reproductive_rate(5, 6e3, "$figure_directory/supplemental/supp_3d.png")
# plot the average V0
plot_vertical_reproductive_ratio(1, 6e3, "$figure_directory/supplemental/supp_3f.png")
# plot the average horizontal transmission
plot_horizontal_and_virulence(1, 6e3, "$figure_directory/supplemental/supp_3h.png")
# plot the evenness through time
plot_evenness(1, 6.1e3, "$figure_directory/supplemental/supp_3j.png")

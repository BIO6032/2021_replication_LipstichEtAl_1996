# plot first 5000 time points (before stabilization)
include("../Model_scripts/Functions.jl");
# plot the number of infected & uninfected hosts
plot_population_numbers(
    hcat(labels...),
    "c=4.0\n \nNumber of infected and uninfected hosts",
    100,
    5000,
    "$figure_directory/supplemental/supp_2b.png"
)
# plot the average weighted R0
plot_reproductive_rate(5, 5000, "$figure_directory/supplemental/supp_2d.png")
# plot the average V0
plot_vertical_reproductive_ratio(1, 5000, "$figure_directory/supplemental/supp_2f.png")
# plot the average weighted horizontal transmission
plot_horizontal_and_virulence(1, 5000, "$figure_directory/supplemental/supp_2h.png")
# plot the evenness through time
plot_evenness(1, 5000, "$figure_directory/supplemental/supp_2j.png")

# plot first 5000 time points supplemental
include("../Model_scripts/Functions.jl");
# plot the number of infected & uninfected hosts
plot_population_numbers(
    hcat(labels...),
    "by=1.0\n \nNumber of infected and uninfected hosts",
    100,
    5e3,
    "$figure_directory/supplemental/supp_1b.png"
)
# plot the average weighted R0
plot_reproductive_rate(6, 5e3, "$figure_directory/supplemental/supp_1d.png")
# plot the average weighted V0
plot_vertical_reproductive_ratio(1, 5e3, "$figure_directory/supplemental/supp_1f.png")
# plot the weighted average mortality rate through time
plot_mortality(1, 5e3, "$figure_directory/supplemental/supp_1h.png")
# plot the evenness through time
plot_evenness(1, 5e3, "$figure_directory/supplemental/supp_1j.png")

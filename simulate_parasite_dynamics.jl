using Plots
include("Fig1_strain_parameters.jl")


#using fonctions defined in Fig1_strain_parameters
#parameters
p = new_Y()
h = new_X()
#model
N = parasites(h, p, 100.0, 20.0, timesteps = 5000, iter = 1)
#graph
using Plots
plot(N, 1, 5001, labels = ["X", "Y"], title = "Dynamics of populations X and Y")

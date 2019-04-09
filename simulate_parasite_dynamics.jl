using Plots
include("Fig1_strain_parameters.jl")

#using fonctions defined in Fig1_strain_parameters
# making a new parasite & a new host
p = new_Y()
h = new_X()
#model
N = parasites(h, p, 100.0, 20.0, timesteps = 50, iter = 1000)
#graph
plot(N, labels = ["X", "Y"], title = "Dynamics of populations X and Y")

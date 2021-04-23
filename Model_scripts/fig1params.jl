using DifferentialEquations
using Plots
import Random

########### Defines the model parameters for Figure 1 ###########

n_parasites = 100; # number of parasites

c = 4.0; # host contact rate

ux = 0.2; # mortality rate of uninfected hosts

# mortality rate of infected hosts
Random.seed!(1479);
ui = rand(200:1000, n_parasites) / 1000; # always >= ux

βy = 3 * (ui .- ux) ./ (ui .- ux .+ 1) # horizontal transmission rate

bx = 1.0 # birth rate of uninfected individuals

ei = fill(bx .- bi, n_parasites) # number of infected offspring from an infected host

# setup the model parameters
parameters = (bx = bx, βy = βy, ei = ei, c = c, K = 100.0, ux = ux, bi = bi, ui = ui);

# define the numbers of individuals in each category
Y = zeros(Float64, n_parasites); # number of parasites of each strain
Y[1] = 1.0; # start with one initial parasite
X0 = 80.0; # number of initial uninfected hosts
new_U = vcat(X0, Y); # merge both together into one array

# define the limits of the simulation
windowstart = 0.0;
windowsize = 1000.0;
windowend = windowstart + windowsize;
N = zeros(Float64, (n_parasites + 1, (n_parasites - 1) * Int(windowsize) + 1));

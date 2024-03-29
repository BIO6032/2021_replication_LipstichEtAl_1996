using DifferentialEquations
using Plots
using Distributions
using Random: Random
########### Defines the model parameters for Figure 3 ###########
n_parasites = 100; # number of parasites
K = 100; # carrying capacity
ux = 0.2; # mortality rate of uninfected hosts
# mortality rate of infected hosts (always >= ux)
Random.seed!(6376);
ui = rand(200:1000, n_parasites) / 1000;
# define 3 random constraints per strain: r1, r2, r3
# r1 = virulence of the strain
# r2 = fraction of virulence attributable to fecundity loss
# r3 = fraction of offspring of infected hosts whichare infected
r1 = rand(Float64, n_parasites);
r2 = rand(Float64, n_parasites);
r3 = rand() * r1 # in this case, r3 is limited by r1
bx = 1.0; # birth rate of uninfected individuals
bi = bx .* r3 .* (1 .- r1 .* r2); # birth rate of infected individuals
α = rand(Float64, n_parasites); # cost of vertical transmission
βy = r1 .- (α .* bi) / bx; # horizontal transmission rate
# ensure no values in βy are negative by setting them to 0
for (index, value) in enumerate(βy)
    if value < 0
        βy[index] = 0
    end
end
ei = bx .* (1 .- r3) .* (1 .- (r1 .* r2)); # number of infected offspring from an infected host
# set up the model parameters
parameters = (bx=bx, βy=βy, ei=ei, c=c, K=K, ux=ux, bi=bi, ui=ui);
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
figure_directory = "Figure3" # output folder for subplots of Figure 3

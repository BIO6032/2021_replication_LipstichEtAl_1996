using DifferentialEquations
using Plots
import Random

n_parasites = 200;

#parameters
c = 4.0;
ux = 0.2;
Random.seed!(4040);
ui = rand(200:1000, n_parasites)/1000; #always >= ux
βy = 3 * (ui .- ux) ./ (ui .- ux .+ 1) # augmente avec la mortalité u
bx = 1.0
ei = fill(bx .- bi, n_parasites)


debut = 0.0;
duree = 1000.0;
fin = debut + duree;
N = zeros(Float64, (n_parasites+1, (n_parasites-1)*Int(duree)+1));
parameters = (bx = bx, βy = βy, ei = ei, c = c, K = 100, ux = ux, bi = bi, ui = ui);

Y = zeros(Float64, n_parasites);
X0 = 10.0;
new_U = vcat(X0, Y);
Y[1] = 1.0;

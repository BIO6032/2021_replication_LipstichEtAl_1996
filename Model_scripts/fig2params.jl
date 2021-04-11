using DifferentialEquations
using Plots
import Random

n_parasites = 200;

ux = 0.2;
Random.seed!(1234);
ui = rand(200:1000, n_parasites)/1000;
Random.seed!(1235);
r1 = rand(Float64, n_parasites) #virulence
Random.seed!(1236);
r2 = rand(Float64, n_parasites)
Random.seed!(1237);
r3 = rand(Float64, n_parasites)
bx = 1.0; #birth rate of uninfected individuals
bi = bx .* r3 .* (1 .- r1 .* r2); #birth rate infected -- disadvantage of being infected (not as many offspring) aka vertical transmission
# V0 = (bi .* ux) ./ (bx .* ui); #introduction of new strain into uninfected population from existing host (beginning new problem: how fast does it spread) #not used.. previously needed to calculate alpha but now we're assigning random values
α = rand(Float64, n_parasites); #cost of vertical transmission
βy = r1 .- (α.* bi) ./ bx
#ensure no values in βy are negative
for (index, value) in enumerate(βy)
   if value < 0
       βy[index] = 0
   end
end


ei = bx .* (1 .- r3) .* (1 .- (r1 .* r2))

Y = zeros(Float64, n_parasites); #number of parasites
Y[1] = 1.0; #1 parasite at first time step
X0 = 10.0; #number of initial hosts


debut = 0.0;
duree = 1000.0;
fin = debut + duree;
N = zeros(Float64, (n_parasites+1, (n_parasites-1)*Int(duree)+1));
new_U = vcat(X0, Y);
parameters = (bx = bx, βy = βy, ei = ei, c = c, K = 100.0, ux = ux, bi = bi, ui = ui);

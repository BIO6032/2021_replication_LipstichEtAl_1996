# Defines the parameters for the population dynamics model

# define a seed (for reproducibility purposes)
Random.seed!(1234)

########### model parameters ###########
n_parasites = 200
c = 4.0
ux = 0.2;
ux1 = fill(0.2, n_parasites); #ux is ui 1000 because of β //different way to have the same results :[0.2 for x in 1:1000]
ui = rand(200:1000, n_parasites)/1000 # must always be >= à ux
u1 = (ui - ux1);
βy = 3 * u1 ./ (u1.+ 1); # increases with mortality u
bx = 1.0;
ei = fill(0.9, n_parasites);    # constant for fig 2 part 1

start = 0.0
length = 1000.0
fin = start + length
N = zeros(Float64, (n_parasites + 1, (n_parasites - 1) * Int(length) + 1))
new_U = vcat(X0, Y)
parameters = (bx = bx, βy = βy, ei = ei, c = c, K = 80.0, ux = ux, by = by, ui = ui)

Y = zeros(Float64, length(ey));
Y[1] = 1.0;
X0 = 10.0;

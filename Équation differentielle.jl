using DifferentialEquations
using Plots

#Figure 2
c = 4.0;
ux = 0.2;
ux1 = fill(0.2, 1000); #le ux et uy 1000 à cause de la β # une autre façon :[0.2 for x in 1:1000]
uy = rand(200:1000, 1000)/1000 #doit être toujours >= à ux
u1 = (uy - ux1);
βy = 3 * u1 ./ (u1.+ 1); # augmente avec la mortalité u
bx = 1.0;
by = 0.1;
#ey = bx - by
ey = fill(0.9, 1000);    # constant for fig 2 part 1

Y = zeros(Float64, length(ey));
Y[1] = 1.0;
X0 = 10.0;

# differential equations formula for finding densities of each pop at next time step
function fonction(u, p, t)
    x = u[1]
    y = u[2:end]
    regul = (1-(sum(u))/(p.K))
    dx = (p.bx*x + sum(p.ey.*y))*regul - p.ux*x-p.c*sum(p.βy.*y)*x
    dy = p.by .* y .* regul .- p.uy .* y .+ p.c .* p.βy .* x .* y
    return vcat(dx, dy)
end

# each strain introduction (1000x) - ORIGINAL ONE
for i in 1:100
    parameters = (bx = bx, βy = βy, ey = ey, c = c, K = 100.0, ux = ux, by = by, uy = uy)
    #i=1.0
    debut = 0.0
    fin = 100.0
    prob = ODEProblem(fonction, vcat(X0, Y), (debut, fin), parameters)
    solution = solve(prob)
#        for k in 2:1000
            for j in 2:100
            debut = fin
            debut = i .*100.0
            fin = i+1 .*100.0
    		new_u = solution[end]
            new_y = findfirst(x -> x == 0.0, new_u)
            new_u[new_y] = 1.0
            prob = ODEProblem(fonction, new_u, (debut, fin), parameters)
            solution = solve(prob)

            for j in 3:100
            debut = fin
            debut = i .*100.0
            fin = i+1 .*100.0
            new_u = solution[end]
            new_y = findfirst(x -> x == 0.0, new_u)
            new_u[new_y] = 1.0
            prob = ODEProblem(fonction, new_u, (debut, fin), parameters)
            solution = solve(prob)
return solution
        end
    end

end

# SANDRINES TEST
debut = 0.0
fin = 10.0
new_U = vcat(X0, Y)
parameters = (bx = bx, βy = βy, ey = ey, c = c, K = 100.0, ux = ux, by = by, uy = uy)
# each strain introduction (1000x)

for i in 1:10
    print(debut)
    # initial conditions
    #prob = ODEProblem(fonction, vcat(X0, Y), (debut, fin), parameters)
    #solution = solve(prob)
    prob = ODEProblem(fonction, new_U, (debut, fin), parameters)
    solution = solve(prob)

    # set conditions & new parasite for next loop
    global new_U = solution[end]
    new_y = findfirst(x -> x == 0.0, new_U)
    new_U[new_y] = 1.0

    # set limits for next loop
    global debut = fin
    global fin = Float64(i+1) #.* 1000.0

    #end
end

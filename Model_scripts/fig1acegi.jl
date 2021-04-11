bi = 0.1;

include("../Model_scripts/fig1params.jl")

include("../Model_scripts/Functions.jl")
run_simulation()

Np = N';

########## Figure1a ##########

lbls = ["" for i = 1:1:n_parasites];
lbls2 = vcat("Infected", lbls);
lbls2 = hcat(lbls2...);

plot(Np[:,2:end], c=:blue, lw=1.5, alpha=0.4, title = "by=0.1\n \nNumber of infected and uninfected hosts",
    xlabel = "Time", ylabel = "Number of individuals",
    label = lbls2, ylims =(0,100))
plot!(Np[:,1], c=:black, lw=1.5, label = "Uninfected")
#plot!(sum(Np[:,2:end]; dims=2), label = "total parasites")

 png("Figure1/graph_1a.png")

###### Figure1g ######

#strains that survive at each time step for ui
survival = (Np.>0.0)[:,2:end];
survived_ui = survival.*ui';
avg_survived = sum(survived_ui; dims=2)./sum(survival; dims=2);

avg_w_survived = sum(Np[:,2:end].*ui'; dims=2)./sum(Np[:,2:end]; dims=2);
ui_avg = avg_w_survived;

# plot(avg_survived, title = "Average ui in the population",xlabel = "Time", ylabel = "Mean mortality (ui)", leg = false)
plot(avg_w_survived, c=:black, lw=1.5,title = "Average mortality in the population",
    xlabel = "Time", ylabel = "Mean mortality (ui)", leg = false, ylims =(0,1))

 png("Figure1/graph_1g.png")



########### Figure1e ############

#strains that survive for βy
survived_βy = survival.*βy';
avg_survived_βy = sum(survived_βy; dims=2)./sum(survival; dims=2);
βi_avg = avg_survived_βy;

#avec moins de bruit (weighted)
avg_w_survived_βy = sum(Np[:,2:end].*βy'; dims=2)./sum(Np[:,2:end]; dims=2);
βy_avg = avg_w_survived_βy;

k=1;
H0 = c*βi_avg./ui_avg.*k.*(1-ux/bx);

H0_w = c*βy_avg./ui_avg.*k.*(1-ux/bx);

V0 = bi*ux./(bx*ui_avg);

plot(V0, c=:black,lw=1.5, title = "Average vertical cases in the population", xlabel = "Time", ylabel = "Mean V0", leg = false, ylims=(0,1))

 png("Figure1/graph_1e")



########## Figure1c ###########

#calculating R0
R0 = H0 + V0;

R0_w = H0_w + V0;

plot(R0_w, c=:black,lw=1.5, title = "Average R0 in the population", xlabel = "Time", ylabel = "Mean R0", leg = false, ylims=(0,7))

 png("Figure1/graph_1c")



######### Figure1i ###########

#calculating evenness
function pielou(n)
    np = filter(x -> x > eps(), n)
    p = np./sum(np)
    ev = length(p) == 1 ? 0.0 : -sum(p.*log.(p))*(1/log(length(p)))
    return ev
end;
ev = mapslices(pielou, Np[:,2:end]; dims=2);
plot(ev, c=:black, lw=1.5, title = "Evenness", xlabel = "Time", ylabel = "Relative abundance (log)", leg = false, ylims = (0,1))

 png("Figure1/graph_1i.png")

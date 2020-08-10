include("fig2params.jl")
by = 0.1;

include("../Functions.jl")
run_simulation()

Np = N'

########## Figure 2a ##########

lbls = ["" for i = 1:1:n_parasites]
lbls2 = vcat("Infected", lbls)
lbls2 = hcat(lbls2...)

plot(Np[:,2:end], c=:blue, lw=0.4, alpha=0.4, title = "Number of infected and uninfected hosts",
    xlabel = "Time", ylabel = "Number of individuals", label = lbls2, ylims =(0,70))
plot!(Np[:,1], c=:black, lw=0.4, label = "Uninfected")
#plot!(sum(Np[:,2:end]; dims=2), label = "total parasites")

# png("Figure 2/graph_2a.png")


###### FIGURE 2g ######

#strains that survive at each time step for uy
survival = (Np.>0.0)[:,2:end];
survived_ui = survival.*uy';
avg_survived = sum(survived_ui; dims=2)./sum(survival; dims=2);

avg_w_survived = sum(Np[:,2:end].*uy'; dims=2)./sum(Np[:,2:end]; dims=2);
uy_avg = avg_w_survived;

# plot(avg_survived, title = "Average uy in the population",xlabel = "Time", ylabel = "Mean mortality (ui)", leg = false)
plot(avg_w_survived, c=:black,title = "Average ui in the population",
    xlabel = "Time", ylabel = "Mean mortality (ui)", leg = false, ylims =(0,1))

# png("Figure 2/graph_2g.png")



########### Figure 2e ############

#strains that survive for βy
survived_βy = survival.*βy';
avg_survived_βy = sum(survived_βy; dims=2)./sum(survival; dims=2);
βi_avg = avg_survived_βy;

#avec moins de bruit (weighted)
avg_w_survived_βy = sum(Np[:,2:end].*βy'; dims=2)./sum(Np[:,2:end]; dims=2);
βy_avg = avg_w_survived_βy;

k=1;
H0 = c*βi_avg./uy_avg.*k.*(1-ux/bx);

H0_w = c*βy_avg./uy_avg.*k.*(1-ux/bx);

V0 = by*ux./(bx*uy_avg);

plot(V0, c=:black, title = "Average V0 in the population", xlabel = "Time", ylabel = "Mean R0", leg = false, ylims=(0,1))

# png("Figure 2/graph_2e")



########## Figure 2c ###########

#calculating R0
R0 = H0 + V0;

R0_w = H0_w + V0;

plot(R0_w, title = "Average R0 in the population", xlabel = "Time", ylabel = "Mean R0", leg = false, ylims=(0,))




######### Figure 2i ###########

#calculating evenness
function pielou(n)
    np = filter(x -> x > eps(), n)
    p = np./sum(np)
    ev = length(p) == 1 ? 0.0 : -sum(p.*log.(p))*(1/log(length(p)))
    return ev
end;
ev = mapslices(pielou, Np[:,2:end]; dims=2);
plot(ev, c=:black, title = "Evenness", xlabel = "Time", ylabel = "Relative abundance (log)", leg = false, ylims = (0,1))

# png("Figure 2/graph_2i.png")

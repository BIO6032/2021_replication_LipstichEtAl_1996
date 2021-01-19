
include("fig2params.jl")
c = 4.0

include("../Functions.jl")
run_simulation()


########### Figure2b ###########

Np = N'

lbls = ["" for i = 1:1:n_parasites];
lbls2 = vcat("Infected", lbls);
lbls2 = hcat(lbls2...);



plot(Np[:,2:end], c=:blue, lw=1.5, alpha=0.4, title = "c=4.0\n \nNumber of infected and uninfected hosts",
    xlabel = "Time", ylabel = "Number of individuals", label = lbls2, ylims = (0,70))
plot!(Np[:,1], c=:black, lw=1.5, label = "Uninfected")
# plot!(sum(Np[:,2:end]; dims=2), label = "Total # parasites")

 png("Figure2/graph_2b.png")


######## Figure2d #########

#strains that survive at each time step for ui
survival = (Np.>0.0)[:,2:end];
survived_ui = survival.*ui';
avg_survived = sum(survived_ui; dims=2)./sum(survival; dims=2);
avg_w_survived = sum(Np[:,2:end].*ui'; dims=2)./sum(Np[:,2:end]; dims=2);
ui_avg = avg_w_survived;

#strains that survive for βy
survived_βy = survival.*βy';
avg_survived_βy = sum(survived_βy; dims=2)./sum(survival; dims=2);
βi_avg = avg_survived_βy; #mean β

#weighted β
avg_w_survived_βy = sum(Np[:,2:end].*βy'; dims=2)./sum(Np[:,2:end]; dims=2);
βy_avg = avg_w_survived_βy;

survived_bi = survival.*bi';
avg_bi_survived = sum(survived_bi; dims=2)./sum(survival; dims=2);;
avg_w_bi = sum(Np[:,2:end].*bi'; dims=2)./sum(Np[:,2:end]; dims=2)
bi_avg = avg_w_bi;

#calculating R0
k=1;
H0 = c*βi_avg./ui_avg.*k.*(1-ux/bx);

H0_w = c*βy_avg./ui_avg.*k.*(1-ux/bx);

V0_w = bi_avg.*ux./(bx*ui_avg);

R0 = H0 + V0_w;

R0_w = H0_w + V0_w;

# plot(R0, title = "Average R0 in the population", xlabel = "Time", ylabel = "Mean R0", leg = false)
plot(R0_w,c=:black, lw=1.5, title = "Average R0 in the population", xlabel = "Time", ylabel = "Mean R0", leg = false, ylims=(0,7))

 png("Figure2/graph_2d")


####### Figure2f #########

V0= bi_avg.*ux./(bx*ui_avg);

plot(V0,c=:black, lw=1.5, title = "Average V0 in the population", xlabel = "Time", ylabel = "Mean V0", leg = false, ylims = (0.0,1.0))
 png("Figure2/graph_2f")


####### Figure2h ########

#mean virulence
Vir = (1 .- (bi .+ ei).*ux./(bx.*ui));

survived_Vir = survival.*Vir';
avg_survived_Vir = sum(survived_Vir; dims=2)./sum(survival; dims=2);
Vir_avg = avg_survived_Vir;

#weighted virulence
avg_w_survived_Vir = sum(Np[:,2:end].*Vir'; dims=2)./sum(Np[:,2:end]; dims=2);
Vir_avg_w = avg_w_survived_Vir;

#plot(βy_avg, title = "Mean virulence and beta", label = "beta")
#plot!(Vir_avg, label = "virulence")
plot(βi_avg,c=:black, title = "Virulence and beta", label = "Beta", xlabel = "Time", ylabel = "Mean virulence & \n Mean Beta", ylims = (0,1))
plot!(Vir_avg_w,c=:blue, lw=1.5, label = "Virulence")

 png("Figure2/graph_2h.png")


######### Figure2j ############

function pielou(n)
    np = filter(x -> x > eps(), n)
    p = np./sum(np)
    ev = length(p) == 1 ? 0.0 : -sum(p.*log.(p))*(1/log(length(p)))
    return ev
end
ev = mapslices(pielou, Np[:,2:end]; dims=2);
plot(ev,c=:black, lw=1.5, title = "Evenness", xlabel = "Time", ylabel = "Relative abundance (log)", leg = false,ylims = (0,1))

 png("Figure2/graph_2j.png")

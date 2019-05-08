#Figure 2-e pour by = 0.1
using Plots

Vertical_case1 = fill(0.1, 100001)
plot(Vertical_case1,c=:black, ylims=(0,1),leg=false,title = "Fraction of all new cases acquired vertically",
xlabel = "Time", ylabel = "Vertical cases")

png("Figure 2/graph_2e.png")

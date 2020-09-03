### Create plot arrangements with multiple panels for Figures 1,2,3 of Virulence project ###

library(png)
library(grid)
library(gridExtra)
library(ggplot2)

########################## FIGURE 1 ###########################
setwd("/Users/sandrinesoeharjono/Documents/2019_replication_LipstichEtAl_1996/Figure 1")

plot1a <- readPNG('graph_1a.png')
plot1b <- readPNG('graph_1b.png')
plot1c <- readPNG('graph_1c.png')
plot1d <- readPNG('graph_1d.png')
plot1e <- readPNG('graph_1e.png')
plot1f <- readPNG('graph_1f.png')
plot1g <- readPNG('graph_1g.png')
plot1h <- readPNG('graph_1h.png')
plot1i <- readPNG('graph_1i.png')
plot1j <- readPNG('graph_1j.png')

grid.arrange(rasterGrob(plot1a),
             rasterGrob(plot1b),
             rasterGrob(plot1c),
             rasterGrob(plot1d),
             rasterGrob(plot1e),
             rasterGrob(plot1f),
             rasterGrob(plot1g),
             rasterGrob(plot1h),
             rasterGrob(plot1j),
             rasterGrob(plot1j),
             ncol=2)

########################## FIGURE 2 ###########################
setwd("/Users/sandrinesoeharjono/Documents/2019_replication_LipstichEtAl_1996/Figure 2")

plot2a <- readPNG('graph_2a.png')
plot2b <- readPNG('graph_2b.png')
plot2c <- readPNG('graph_2c.png')
plot2d <- readPNG('graph_2d.png')
plot2e <- readPNG('graph_2e.png')
plot2f <- readPNG('graph_2f.png')
plot2g <- readPNG('graph_2g.png')
plot2h <- readPNG('graph_2h.png')
plot2i <- readPNG('graph_2i.png')
plot2j <- readPNG('graph_2j.png')

grid.arrange(rasterGrob(plot2a),
             rasterGrob(plot2b),
             rasterGrob(plot2c),
             rasterGrob(plot2d),
             rasterGrob(plot2e),
             rasterGrob(plot2f),
             rasterGrob(plot2g),
             rasterGrob(plot2h),
             rasterGrob(plot2j),
             rasterGrob(plot2j),
             ncol=2)

########################## FIGURE 3 ###########################
setwd("/Users/sandrinesoeharjono/Documents/2019_replication_LipstichEtAl_1996/Figure 3")

plot3a <- readPNG('graph_3a.png')
plot3b <- readPNG('graph_3b.png')
plot3c <- readPNG('graph_3c.png')
plot3d <- readPNG('graph_3d.png')
plot3e <- readPNG('graph_3e.png')
plot3f <- readPNG('graph_3f.png')
plot3g <- readPNG('graph_3g.png')
plot3h <- readPNG('graph_3h.png')
plot3i <- readPNG('graph_3i.png')
plot3j <- readPNG('graph_3j.png')

grid.arrange(rasterGrob(plot3a),
             rasterGrob(plot3b),
             rasterGrob(plot3c),
             rasterGrob(plot3d),
             rasterGrob(plot3e),
             rasterGrob(plot3f),
             rasterGrob(plot3g),
             rasterGrob(plot3h),
             rasterGrob(plot3j),
             rasterGrob(plot3j),
             ncol=2)

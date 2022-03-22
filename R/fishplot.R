library(tidyverse)
library(fishplot)

# # Principle
# timepoints=c(0,30,75,150)      
# 
# #provide a matrix with the fraction of each population
# #present at each timepoint
# frac.table = matrix(
#   c(100, 45, 00, 00,
#     02, 00, 00, 00,
#     02, 00, 02, 01,
#     98, 00, 95, 40),
#   ncol=length(timepoints))
# 
# parents = c(0,1,1,3)
# 
# #create a fish object
# fish = createFishObject(frac.table,parents,timepoints=timepoints)
# 
# #calculate the layout of the drawing
# fish = layoutClones(fish)
# 
# #draw the plot, using the splining method (recommended)
# #and providing both timepoints to label and a plot title
# fishPlot(fish,shape="spline",title.btm="Sample1",
#          cex.title=1, vlines=c(0,150), 
#          vlab=c("day 0","day 150"))
library(readr)
hct634851 <- read_csv("hct634851.csv")

# hct634851_data <- 
  hct634851 %>% 
  pivot_longer(cols = c(Date_recipient, VAF_recipient))

  hct634851 %>% 
    pivot_longer(cols = c(Date_recipient))
  

a <- "hct634851"

timepoints1=c(as.numeric(hct634851$Time_from_HCT))      
frac.table = matrix( # 28.5, 0.06, 0.0 # 7.5, 84.1, 9.8
  c(0.0001, 0.0001, 
    0.0001, 0.0001,
    4.1, 0.0001,
    4.1, 17),
  ncol=length(timepoints1))
parents = c(0,0) # "RUNX1 P425fs", "ASXL1 G646Wfs*12"
fish = createFishObject(frac.table,parents,timepoints=timepoints1,
                        col = c("darkblue", "red"),
                        clone.labels = c("RUNX1 P425fs", "ASXL1 G646Wfs*12"))
fish = layoutClones(fish)

fishPlot(fish,shape="spline",title.btm= hct634851$HCT_id,
         cex.title=1, vlines=c(timepoints1), col.vline = "grey",
         vlab=c(hct634851$Time_from_HCT), 
         bg.type = "solid",
         bg.col = "white"
         )
drawLegend(fish,cex=1.5,xpos=-40)



# #Case1
# timepoints=c(0,150)      
# frac.table = matrix( # 28.5, 0.06, 0.0 # 7.5, 84.1, 9.8
#   c(28.5, 0.06, 0.0001,
#     7.5, 84.1, 9.8),
#   ncol=length(timepoints))
# parents = c(0,0, 2) # "TET2", "TP53", "KMT2A"
# fish = createFishObject(frac.table,parents,timepoints=timepoints,
#                         col = c("darkblue", "red", "gold"),
#                         clone.labels = c("TET2", "TP53", "KMT2A"))
# fish = layoutClones(fish)
# 
# # jpeg("Case1 fishplot.jpeg", width = 960, height = 480)
# fishPlot(fish,shape="spline",title.btm="Case1",
#          cex.title=1, vlines=c(0,150), col.vline = "grey",
#          vlab=c("ARCH","T"), bg.type = "solid",
#          bg.col = "white")
# drawLegend(fish,cex=1.5,xpos=-40)
# # dev.off()
# 
# #Case3
# timepoints=c(0,150)      
# frac.table = matrix( 
#   c(14.6,
#     53.1),
#   ncol=length(timepoints))
# parents = c(0) # "TP53"
# fish = createFishObject(frac.table,parents,timepoints=timepoints,
#                         col = c("red"),
#                         clone.labels = c("TP53"))
# fish = layoutClones(fish)
# 
# jpeg("Case3 fishplot.jpeg", width = 960, height = 480)
# fishPlot(fish,shape="spline",title.btm="Case3",
#          cex.title=1, vlines=c(0,150), col.vline = "grey",
#          vlab=c("ARCH","T"), bg.type = "solid",
#          bg.col = "white")
# drawLegend(fish,cex=1.5,xpos=-40)
# dev.off()
# 
# #Case6
# timepoints=c(0,150)      
# frac.table = matrix( 
#   c(3, 2, 0.0001,
#     49, 23, 14),
#   ncol=length(timepoints))
# parents = c(0,1,1) # "TP53", "PTPN11", "RUNX1", "KRAS", , "RUNX1"
# fish = createFishObject(frac.table,parents,timepoints=timepoints,
#                         col = c("red", "grey", "#6300A7FF"#, #6300A7FF" "rosybrown1", "lightgreen"
#                         ),
#                         clone.labels = c("TP53", "KRAS", "PTPN11")) # "sandybrown", "lightslateblue" 
# fish = layoutClones(fish)
# 
# jpeg("Case6 fishplot 2.jpeg", width = 960, height = 480)
# fishPlot(fish,shape="spline",title.btm="Case6",
#          cex.title=1, vlines=c(0,150), col.vline = "grey",
#          vlab=c("ARCH","T"), bg.type = "solid",
#          bg.col = "white")
# drawLegend(fish,cex=1.5,xpos=-40)
# dev.off()
# 
# #Case7
# timepoints=c(0,150)      
# frac.table = matrix( 
#   c(3.7, 3.9,
#     38.0, 40.0),
#   ncol=length(timepoints))
# parents = c(0,0) # "TP53", "TP53"
# fish = createFishObject(frac.table,parents,timepoints=timepoints,
#                         col = c("red", "red"),
#                         clone.labels = c("TP53", "TP53"))
# fish = layoutClones(fish)
# 
# jpeg("Case7 fishplot.jpeg", width = 960, height = 480)
# fishPlot(fish,shape="spline",title.btm="Case7",
#          cex.title=1, vlines=c(0,150), col.vline = "grey",
#          vlab=c("ARCH","T"), bg.type = "solid",
#          bg.col = "white")
# drawLegend(fish,cex=1.5,xpos=-40)
# dev.off()
# 
# #Case10
# timepoints=c(0,150)      
# frac.table = matrix( 
#   c(21.4, 4.9, 2.5,
#     37.1, 0.0, 35.7),
#   ncol=length(timepoints))
# parents = c(0,1, 1) # "U2QR1", "BCOR", "TEV6"
# fish = createFishObject(frac.table,parents,timepoints=timepoints,
#                         col = c("grey", "chocolate2", "thistle2"),
#                         clone.labels = c("U2QR1", "BCOR", "TEV6"))
# fish = layoutClones(fish)
# 
# jpeg("Case10 fishplot.jpeg", width = 960, height = 480)
# fishPlot(fish,shape="spline",title.btm="Case10",
#          cex.title=1, vlines=c(0,150), col.vline = "grey",
#          vlab=c("ARCH","T"), bg.type = "solid",
#          bg.col = "white")
# drawLegend(fish,cex=1.5,xpos=-40)
# dev.off()
# 
# #Case13
# timepoints=c(0,150)      
# frac.table = matrix( 
#   c(3.3,
#     27.0),
#   ncol=length(timepoints))
# parents = c(0) # "TET2"
# fish = createFishObject(frac.table,parents,timepoints=timepoints,
#                         col = c("darkblue"),
#                         clone.labels = c("TET2"))
# fish = layoutClones(fish)
# 
# jpeg("Case13 fishplot.jpeg", width = 960, height = 480)
# fishPlot(fish,shape="spline",title.btm="Case13",
#          cex.title=1, vlines=c(0,150), col.vline = "grey",
#          vlab=c("ARCH","T"), bg.type = "solid",
#          bg.col = "white")
# drawLegend(fish,cex=1.5,xpos=-40)
# dev.off()
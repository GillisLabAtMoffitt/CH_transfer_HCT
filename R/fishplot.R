library(tidyverse)
library(lubridate)
library(fishplot)

################################################################################# I ### Load data
path <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", "CH_transfer_HCT")

WES_data <- 
  readxl::read_xlsx(paste0(path, "/data/CICPT2144_cleaned vcf_01.21.22_reviewed w EHR and raw vcf_03.30.22_for Christelle.xlsx"))

# hct634851 <- read_csv("hct634851.csv")
################################################################################# II ### Data cleaning
WES_data1 <- WES_data %>% 
  # Create a HCT id
  mutate(hct_id = paste0("hct_", dense_rank(BMT_date)), .before = 1) %>% 
  # Keep only sequencing from clinical when have 2 platform data for the same date
  arrange(hct_id, Sample_Type, Date_of_Collection, desc(NGS_Type)) %>% 
  distinct(hct_id, Sample_Type, Date_of_Collection, GENE, VARIANT_P, .keep_all = TRUE)



################################################################################# II ### Fish plots

hct_id <- "hct_1"

timepoints <- c(-38, 81, 368, 368, 1011)    

# -38	NA	NA	NA
# 81	ASXL1	G642fs	0.397
# 368	ASXL1	G642fs	0.481
# 368	RUNX1	P425fs	0.041
# 1011	ASXL1	G646Wfs*12	0.17

frac.table = matrix(
  c(0.0001, 0.0001, 
    0.397, 0.0001,
    0.481, 0.0001,
    4.1, 0.041,
    0.17, 0.0001),
  ncol=length(timepoints))
parents = c(0,0)
fish = createFishObject(frac.table,parents,timepoints=timepoints,
                        col = c("darkblue", "red"),
                        clone.labels = c("ASXL1 G642fs", "RUNX1 P425fs"))
fish = layoutClones(fish)

fishPlot(fish,shape="spline",title.btm= hct_id,
         cex.title=1, vlines=c(timepoints), col.vline = "grey",
         vlab=c(timepoints), 
         bg.type = "solid",
         bg.col = "white"
)
drawLegend(fish,cex=1.5,xpos=-40)

##################
hct_id <- "hct_2"

timepoints <- c(-22, 377, 516)    

# -22	NA	NA	NA
# 377	ASXL1	L956fs	0.116
# 516	ASXL1	L956fs	0.01

frac.table = matrix(
  c(0.0001,
    0.116,
    0.01),
  ncol=length(timepoints))
parents = c(0)
fish = createFishObject(frac.table,parents,timepoints=timepoints,
                        col = c("darkblue"),
                        clone.labels = c("AASXL1	L956fs"))
fish = layoutClones(fish)

fishPlot(fish,shape="spline",title.btm= hct_id,
         cex.title=1, vlines=c(timepoints), col.vline = "grey",
         vlab=c(timepoints), 
         bg.type = "solid",
         bg.col = "white"
)
drawLegend(fish,cex=1.5,xpos=-40)

#################

hct_id <- "hct_3"

timepoints <- c(-132, -29, 30, 115, 231, 337)    

# -132	NA	NA	NA
# -29	NF1	R1362X	0.364
# -29	CHEK2	R346C	0.297
# 30	TET2	Q626X	0.125
# 115	TET2	Q626X	0.111
# 231	TET2	Q626X	0.123
# 337	TET2	Q626X	0.124

frac.table = matrix(
  c(0.364, 0.0001, 0.0001,
    0.0001, 0.297, 0.0001,
    0.0001, 0.0001, 0.125,
    0.0001, 0.0001, 0.111,
    0.0001, 0.0001, 0.123,
    0.0001, 0.0001, 0.124),
  ncol=length(timepoints))
parents = c(0,0, 0)
fish = createFishObject(frac.table,parents,timepoints=timepoints,
                        col = c("darkblue", "red", "yellow"),
                        clone.labels = c("NF1	R1362X", "CHEK2	R346C", "TET2	Q626X"))
fish = layoutClones(fish)

fishPlot(fish,shape="spline",title.btm= hct_id,
         cex.title=1, vlines=c(timepoints), col.vline = "grey",
         vlab=c(timepoints), 
         bg.type = "solid",
         bg.col = "white"
)
drawLegend(fish,cex=1.5,xpos=-40)

#################

hct_id <- "hct_4"

timepoints <- c(-233, -23, 86, 363, 727)    

# -233	IDH2	R140Q	0.4
# -233	SRSF2	P95R	0.399
# -233	SETBP1	D868N	0.109
# -233	GNAS	R844H	0.378
# -23	ASXL1	H630fs	0.163
# -23	IDH2	R140Q	0.356
# -23	SRSF2	P95R	0.369
# -23	GNAS	R844H	0.366
# 86	ASXL1	G462fs	0.115
# 363	ASXL1	G462fs	0.268
# 727	ASXL1	G462fs	0.19

frac.table = matrix(
  c(40, 39.9, 10.9, 37.8, 0.0001, 0.0001, 
    35.6, 36.9, 0.0001, 36.6, 16.3, 0.0001, 
    0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 11.5,
    0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 26.8,
    0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 19),
  ncol=length(timepoints))
parents = c(0,0, 0, 2, 0, 0)
fish = createFishObject(frac.table,parents,timepoints=timepoints,
                        col = c("darkblue", "red", "yellow", "pink", "grey", "green"),
                        clone.labels = c("IDH2	R140Q", "SRSF2	P95R", "SETBP1	D868N", "GNAS	R844H", "ASXL1	H630fs", "ASXL1	G462fs"))
fish = layoutClones(fish)

fishPlot(fish,shape="spline",title.btm= hct_id,
         cex.title=1, vlines=c(timepoints), col.vline = "grey",
         vlab=c(timepoints), 
         bg.type = "solid",
         bg.col = "white"
)
drawLegend(fish,cex=1.5,xpos=-40)

#################

hct_id <- "hct_5"

timepoints <- c(-140, -42, 142, 292, 365, 463, 544)    

# -140	DNMT3A	R882C	0.351
# -140	DNMT3A		0.28
# -140	TET2	M1028Nfs*15	0.32
# -140	SRSF2	P95R	0.05
# -140	SRSF2	P95L	0.24
# -42	NA	NA	NA
# 142	SH2B3	480_487del	0.204
# 142	DNMT3A	R882C	0.108
# 142	DNMT3A		0.121
# 142	TET2	M1028Nfs*15	0.1
# 142	SRSF2	P95L	0.1
# 292	DNMT3A	R882C	0.482
# 292	DNMT3A		0.493
# 292	TET2	M1028Nfs*15	0.45
# 292	SRSF2	P95L	0.506
# 292	BCOR	L1646Pfs*6	0.9
# 365	SH2B3	480_487del	0.077
# 365	DNMT3A	R882C	0.333
# 365	DNMT3A		0.374
# 365	TET2	M1028Nfs*15	0.26
# 365	SRSF2	P95L	0.375
# 463	DNMT3A	R882C	0.171
# 463	DNMT3A		0.19
# 463	TET2	M1028Nfs*15	0.18
# 463	SRSF2	P95L	0.18
# 463	NRAS	Q61H	0.1
# 463	BCOR	L1646Pfs*6	0.13
# 544	SH2B3	480_487del	0.041
# 544	DNMT3A	R882C	0.311
# 544	DNMT3A		0.296
# 544	TET2	M1028Nfs*15	0.29
# 544	SRSF2	P95L	0.346
# 544	BCOR	L1646Pfs*6	0.41
# 544	NRAS	Q61H	0.25

frac.table = matrix(
  c(35, 28, 32, 5, 24, 0.0001, 0.0001, 0.0001, #-140
    0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, #-42       ######################Cannot do that
    10.8, 12.1, 10, 0.0001, 10, 20.4, 0.0001, 0.0001, #142
    48.2, 49.3, 45, 0.0001, 50.6, 0.0001, 90, 0.0001, #292
    33.3, 37.4, 26, 0.0001, 37.5, 7.7, 0.0001, 0.0001, #365
    17.1, 19, 18, 0.0001, 18, 0.0001, 13, 10, #463
    31.1, 29.6, 29, 0.0001, 34,6, 4.1, 41, 25
    ),
  ncol=length(timepoints))
parents = c(0,0, 0, 0, 0, 0)
fish = createFishObject(frac.table,parents,timepoints=timepoints,
                        col = c("darkblue", "red", "yellow", "pink", "grey", "green", "white", "black"),
                        clone.labels = c("DNMT3A	R882C", "DNMT3A", "TET2	M1028Nfs*15", "SRSF2	P95R", "SRSF2	P95L", "SH2B3	480_487del", "BCOR	L1646Pfs*6", "NRAS	Q61H"))
fish = layoutClones(fish)

fishPlot(fish,shape="spline",title.btm= hct_id,
         cex.title=1, vlines=c(timepoints), col.vline = "grey",
         vlab=c(timepoints), 
         bg.type = "solid",
         bg.col = "white"
)
drawLegend(fish,cex=1.5,xpos=-40)

#################

hct_id <- "hct_6"

timepoints <- c(-114, -16, 42, 84, 346, 514, 684)    

# -114	DNMT3A	E578X	0.459
# -114	NRAS	G12D	0.221
# -114	IDH1	R132C	0.221
# -16	DNMT3A	R882H	0.094
# -16	DNMT3A	E578X	0.051
# 42	DNMT3A	P904L	0.179
# 42	KRAS	T50I	0.199
# 42	RUNX1	I177V	0.116
# 84	DNMT3A	P904L	0.181
# 346	DNMT3A	P904L	0.156
# 346	KRAS	T50I	0.155
# 346	RUNX1	I177V	0.104
# 514	DNMT3A	P904L	0.17
# 514	KRAS	T50I	0.15
# 514	RUNX1	I177V	0.1
# 684	DNMT3A	P904L	0.129
# 684	KRAS	T50I	0.139
# 684	RUNX1	I177V	0.097

frac.table = matrix(
  c(45.9, 22.1, 22.1, 0.0001, 0.0001, 0.0001, 0.0001, #-114
    5.1, 0.0001, 0.0001, 9.4, 0.0001, 0.0001, 0.0001, #-16
    0.0001, 0.0001, 0.0001, 0.0001, 17.9, 19.9, 11.6, #42
    0.0001, 0.0001, 0.0001, 0.0001, 18.1, 0.0001, 0.0001, #84
    0.0001, 0.0001, 0.0001, 0.0001, 15.6, 15.5, 10.4, #346
    0.0001, 0.0001, 0.0001, 0.0001, 17, 15, 10, #514
    0.0001, 0.0001, 0.0001, 0.0001, 12.9, 13.9, 9.7),
  ncol=length(timepoints))
parents = c(0,0, 0, 0, 0, 0)
fish = createFishObject(frac.table,parents,timepoints=timepoints,
                        col = c("darkblue", "red", "yellow", "pink", "grey", "green", "white"),
                        clone.labels = c("DNMT3A	E578X", "NRAS	G12D", "IDH1	R132C", "DNMT3A	R882H", "DNMT3A	P904L", "KRAS	T50I", "RUNX1	I177V"))
fish = layoutClones(fish)

fishPlot(fish,shape="spline",title.btm= hct_id,
         cex.title=1, vlines=c(timepoints), col.vline = "grey",
         vlab=c(timepoints), 
         bg.type = "solid",
         bg.col = "white"
)
drawLegend(fish,cex=1.5,xpos=-40)

#################








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
library(tidyverse)
# install.packages("devtools")
library(devtools)
# install_github("chrisamiller/fishplot")
library(fishplot)

# ################################################################################# I ### Load data
# path <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", "CH_transfer_HCT")
# 
# WES_data <- 
#   readxl::read_xlsx(paste0(path, "/data/CICPT2144_cleaned vcf_01.21.22_reviewed w EHR and raw vcf_03.30.22_for Christelle.xlsx"))
# 
# # hct634851 <- read_csv("hct634851.csv")
# ################################################################################# II ### Data cleaning
# WES_data1 <- WES_data %>% 
#   # Create a HCT id
#   mutate(hct_id = paste0("hct_", dense_rank(BMT_date)), .before = 1) %>% 
#   # Keep only sequencing from clinical when have 2 platform data for the same date
#   arrange(hct_id, Sample_Type, Date_of_Collection, desc(NGS_Type)) %>% 
#   distinct(hct_id, Sample_Type, Date_of_Collection, GENE, VARIANT_P, .keep_all = TRUE)


# IDENTIFY CO-OCCURENCE AND MUTUAL EXCLUSIVITY in clones apparition
# In https://www.cbioportal.org/
# Select the cancer study (here acute myeloid and Myelodysplastic)
# Query by gene ex: ASXL1 RUNX1
# In the Mutual exclusivity tab, will show co-occurence or mutual exclusivity





# ################################################################################# II ### Fish plots

hct_id <- "hct_634851"

timepoints <- c(-38, 0, 81, 368, 1011)    

# Donor 0	ASXL1	G642fs	0.111
# -38	NA	NA	NA
# 81	ASXL1	G642fs	0.397
# 368	ASXL1	G642fs	0.481
# 368	RUNX1	P425fs	0.041
# 1011	ASXL1	G646Wfs*12	0.17    C

frac.table = matrix(
  c(0, 0, #-38
    0.0001, 0, #0
    39.7, 0.0001, #81
    48.1, 4.1, #368
    17, 0.0001),
  ncol=length(timepoints))
parents = c(0, 1)
fish = createFishObject(frac.table,parents,timepoints=timepoints,
                        col = c("red", "#ffff3f"),
                        clone.labels = c("ASXL1 G642fs", "RUNX1 P425fs"))
fish = layoutClones(fish)

png("hct_634851 recipient.png", width = 660, height = 480, bg = "transparent")
fishPlot(fish,shape="spline",title.btm= hct_id,
         cex.title=1, vlines=c(timepoints), col.vline = "gainsboro",
         vlab=c(timepoints), 
         pad.left = 0,
         bg.type = "solid",
         bg.col = "transparent"
)
drawLegend(fish,cex=1.0,xpos=fish@timepoints[1], nrow = 4)
dev.off()

timepoints <- c(-100, 0, 100)    
frac.table = matrix(
  c(11.1, 
    11.1,
    11.1),
  ncol=length(timepoints))
parents = c(0)
fish = createFishObject(frac.table,parents,timepoints=timepoints,
                        col = c("red"),
                        clone.labels = c("Donor-ASXL1 G642fs"))
fish = layoutClones(fish)


png("hct_634851 donor.png", width = 660, height = 480, bg = "transparent")
fishPlot(fish,shape="spline",title.btm= hct_id,
         cex.title=1, vlines=NULL,
         ramp.angle = 1,
         pad.left = 0.05,
         vlab=NULL, 
         bg.type = "solid",
         bg.col = "transparent"
)
drawLegend(fish,cex=1.0,xpos=fish@timepoints[1])
dev.off()
##################
hct_id <- "hct_657119"

timepoints <- c(-22, 0, 377, 516)    

# Donor 0	ASXL1	L956fs	0.023
# -22	NA	NA	NA
# 377	ASXL1	L956fs	0.116
# 516	ASXL1	L956fs	0.01

frac.table = matrix(
  c(0,
    0.0001,
    11.6, 
    1),
  ncol=length(timepoints))
parents = c(0)
fish = createFishObject(frac.table,parents,timepoints=timepoints,
                        col = c("red"),
                        clone.labels = c("ASXL1 L956fs"))
fish = layoutClones(fish)

png("hct_657119 recipient.png", width = 660, height = 480, bg = "transparent")
fishPlot(fish,shape="spline",title.btm= hct_id,
         cex.title=1, vlines=c(timepoints), col.vline = "gainsboro",
         pad.left = 0,
         vlab=c(timepoints), 
         bg.type = "solid",
         bg.col = "transparent"
)
drawLegend(fish,cex=1.0,xpos=fish@timepoints[1], nrow = 4)
dev.off()

timepoints <- c(-100, 0, 100)    
frac.table = matrix(
  c(2.3, 
    2.3,
    2.3),
  ncol=length(timepoints))
parents = c(0)
fish = createFishObject(frac.table,parents,timepoints=timepoints,
                        col = c("red"),
                        clone.labels = c("Donor-ASXL1 L956fs"))
fish = layoutClones(fish)

png("hct_657119 donor.png", width = 660, height = 480, bg = "transparent")
fishPlot(fish,shape="polygon",title.btm= hct_id,
         cex.title=1, vlines=NULL,
         ramp.angle = 1,
         pad.left = 0.05,
         vlab=NULL, 
         bg.type = "solid",
         bg.col = "transparent"
)
drawLegend(fish,cex=1.0,xpos=fish@timepoints[1])
dev.off()
#################

hct_id <- "hct_705854"

timepoints <- c(-132, -29, 0, 30, 115, 231, 337)    

# Donor 0	TET2	Q626X	0.86
# -132	NA	NA	NA        C
# -29	NF1	R1362X	0.364
# -29	CHEK2	R346C	0.297
# 30	TET2	Q626X	0.125   C
# 115	TET2	Q626X	0.111
# 231	TET2	Q626X	0.123   C
# 337	TET2	Q626X	0.124   C

frac.table = matrix(
  c(0.0001, 0.0001, 0,
    36.4, 29.7, 0, #-29
    0, 0, 0.0001, #0
    0, 0, 12.5,
    0, 0, 11.1,
    0, 0, 12.3,
    0, 0, 12.4),
  ncol=length(timepoints))
# CHEK2	NF1	0.695	0.695	Mutual exclusivity
parents = c(0,1, 0)
fish = createFishObject(frac.table,parents,timepoints=timepoints,
                        col = c("black", "grey90", 
                                "red"),
                        clone.labels = c("NF1 R1362X", "CHEK2 R346C", 
                                         "TET2 Q626X"))
fish = layoutClones(fish)

png("hct_705854 recipient.png", width = 660, height = 480, bg = "transparent")
fishPlot(fish,shape="spline",title.btm= hct_id,
         cex.title=1, vlines=c(timepoints), col.vline = "gainsboro",
         vlab=c(timepoints), 
         bg.type = "solid",
         bg.col = "transparent"
)
drawLegend(fish,cex=1.0,xpos=fish@timepoints[1], nrow = 4)
dev.off()

timepoints <- c(-100, 0, 100)    
frac.table = matrix(
  c(8.6, 
    8.6,
    8.6),
  ncol=length(timepoints))
parents = c(0)
fish = createFishObject(frac.table,parents,timepoints=timepoints,
                        col = c("red"),
                        clone.labels = c("Donor-TET2 Q626X"))
fish = layoutClones(fish)

png("hct_705854 donor.png", width = 660, height = 480, bg = "transparent")
fishPlot(fish,shape="polygon",title.btm= hct_id,
         cex.title=1, vlines=NULL,
         ramp.angle = 1,
         pad.left = 0.05,
         vlab=NULL, 
         bg.type = "solid",
         bg.col = "transparent"
)
drawLegend(fish,cex=1.0,xpos=fish@timepoints[1])
dev.off()
#################

hct_id <- "hct_715231"

timepoints <- c(-233, -23, 0,  86, 363, 727)    

# Donor 0	ASXL1	G462fs	0.062
# -233	IDH2	R140Q	0.4   C
# -233	SRSF2	P95R	0.399   C
# -233	SETBP1	D868N	0.109   C   lost
# -233	GNAS	R844H	0.378   C
# -23	ASXL1	H630fs	0.163   C
# -23	IDH2	R140Q	0.356   C
# -23	SRSF2	P95R	0.369   C
# -23	GNAS	R844H	0.366   C
# 86	ASXL1	G462fs	0.115
# 363	ASXL1	G462fs	0.268
# 727	ASXL1	G462fs	0.19    C

frac.table = matrix(
  c(0,       40, 39.9, 10.9, 37.8, 0.0001,
    0,       37, 36.9, 0.0001, 36.6, 16.3, # change 35.6 to 37
    0.0001,  0, 0, 0, 0, 0, #0
    11.5,    0, 0, 0, 0, 0, 
    26.8,    0, 0, 0, 0, 0, 
    19,      0, 0, 0, 0, 0),
  ncol=length(timepoints))
# ASXL1	SRSF2	<0.001	<0.001	 Co-occurrence
# IDH2	SRSF2	<0.001	<0.001	 Co-occurrence
# SRSF2	GNAS	<0.001	<0.001	 Co-occurrence
# ASXL1	IDH2	<0.001	<0.001	 Co-occurrence
# ASXL1	GNAS	0.342	0.685	 Co-occurrence
# IDH2	GNAS	0.541	0.878	 Co-occurrence
# IDH2	SERBP1	0.645	0.878	 Mutual exclusivity
# SRSF2	SERBP1	0.750	0.878	 Mutual exclusivity
# ASXL1	SERBP1	0.790	0.878	 Mutual exclusivity
# SERBP1	GNAS	0.970	0.970	 Mutual exclusivity
parents = c(0, 0, 2, 5, 3, 5)
fish = createFishObject(frac.table,parents,timepoints=timepoints,
                        col = c("red", "darkorchid", 
                                "green", "sienna1", 
                                "grey67", "blue"),
                        clone.labels = c("ASXL1 G462fs", "IDH2 R140Q", 
                                         "SRSF2 P95R", "SETBP1 D868N", 
                                         "GNAS R844H", "ASXL1 H630fs"),
                        fix.missing.clones=TRUE)
fish = layoutClones(fish)

png("hct_715231 recipient.png", width = 660, height = 480, bg = "transparent")
fishPlot(fish,shape="spline",title.btm= hct_id,
         cex.title=1, vlines=c(timepoints), col.vline = "gainsboro",
         vlab=c(timepoints), 
         bg.type = "solid",
         bg.col = "transparent"
)
drawLegend(fish,cex=1.0,xpos=fish@timepoints[1], nrow = 4)
dev.off()

timepoints <- c(-100, 0, 100)    
frac.table = matrix(
  c(6.2, 
    6.2,
    6.2),
  ncol=length(timepoints))
parents = c(0)
fish = createFishObject(frac.table,parents,timepoints=timepoints,
                        col = c("red"),
                        clone.labels = c("Donor-ASXL1 G462fs"))
fish = layoutClones(fish)

png("hct_715231 donor.png", width = 660, height = 480, bg = "transparent")
fishPlot(fish,shape="polygon",title.btm= hct_id,
         cex.title=1, vlines=NULL,
         ramp.angle = 1,
         pad.left = 0.05,
         vlab=NULL, 
         bg.type = "solid",
         bg.col = "transparent"
)
drawLegend(fish,cex=1.0,xpos=fish@timepoints[1])
dev.off()
#################

hct_id <- "hct_753399 represented \nat 0.5 of original size"

timepoints <- c(-140, -42, 0, 142, 292, 365, 463, 544)    

# Donor 0	SH2B3	480_487del	0.362
# -140	DNMT3A	R882C	0.351
# -140	DNMT3A		0.28
# -140	TET2	M1028Nfs*15	0.32
# -140	SRSF2	P95R	0.05
# -140	SRSF2	P95L	0.24
# -42	NA	NA	NA                  A
# 142	SH2B3	480_487del	0.204     A
# 142	DNMT3A	R882C	0.108
# 142	DNMT3A		0.121
# 142	TET2	M1028Nfs*15	0.1
# 142	SRSF2	P95L	0.1
# 292	DNMT3A	R882C	0.482
# 292	DNMT3A		0.493
# 292	TET2	M1028Nfs*15	0.45
# 292	SRSF2	P95L	0.506
# 292	BCOR	L1646Pfs*6	0.9
# 365	SH2B3	480_487del	0.077     A
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
# 544	SH2B3	480_487del	0.041     A
# 544	DNMT3A	R882C	0.311
# 544	DNMT3A		0.296
# 544	TET2	M1028Nfs*15	0.29
# 544	SRSF2	P95L	0.346
# 544	BCOR	L1646Pfs*6	0.41
# 544	NRAS	Q61H	0.25

frac.table = matrix(
  c(35.1, 28, "27.9", 5, 24,                     0,           0, 0, #-140
    0.0001, 0.0001, 0.0001, 0.00001, 0.0001,     0,           0, 0, #-42
    
    0, 0, 0, 0, 0,                               0,           0, 0, #0
    
    10.8, "10.7", 10, 0.0001, 10,                20.4,        0.0001, 0.0001, #142
    48.2, "48.1", 45, 0.0001, 50.6,            0.0001,      90, 0.0001, #292
    33.3, "33.1", 26, 0.0001, 37.5,            7.7,         0.0001, 0.0001, #365
    "18.1", "18", "17.9", 0.0001, 18,        0.0001,      13, 10, #463
    31.1, 29.6, 29, 0.0001, 34.6,              4.1,         41, 25),
  ncol=length(timepoints))
frac.table <- matrix(as.numeric(frac.table), ncol=length(timepoints))
frac.table <- 0.5 * frac.table

# TET2	SRSF2	<0.001	<0.001 Co-occurrence
# DNMT3A	SRSF2	<0.001	<0.001	Mutual exclusivity
# SRSF2	BCOR	<0.001	<0.001	Co-occurrence
# TET2	NRAS	<0.001	0.002	Mutual exclusivity
# DNMT3A	TET2	0.128	0.253	Mutual exclusivity
# BCOR	NRAS	0.160	0.253	Co-occurrence
# DNMT3A	BCOR	0.177	0.253	Co-occurrence
# DNMT3A	NRAS	0.241	0.301	Co-occurrence
# TET2	BCOR	0.297	0.331	Mutual exclusivity
# SRSF2	NRAS	0.500	0.500	Mutual exclusivity
parents = c(0,1, 2, 5, 0, 0, 0, 7)
fish = createFishObject(frac.table,parents,timepoints=timepoints, fix.missing.clones=TRUE,
                        col = c("darkblue", "slategray1", 
                                "darkorange1", "green", 
                                "forestgreen", "red", 
                                "purple4", "gold"),
                        clone.labels = c("DNMT3A R882C", "DNMT3A", 
                                         "TET2 M1028Nfs*15", "SRSF2 P95R", 
                                         "SRSF2 P95L", "SH2B3 480_487del", 
                                         "BCOR L1646Pfs*6", "NRAS Q61H"
                                         ))
fish = layoutClones(fish)

png("hct_753399 recipient.png", width = 660, height = 480, bg = "transparent")
fishPlot(fish,shape="spline",title.btm= hct_id,
         cex.title=1, vlines=c(timepoints), col.vline = "gainsboro",
         vlab=c(timepoints), 
         bg.type = "solid",
         bg.col = "transparent"
)
drawLegend(fish,cex=1.0,xpos=fish@timepoints[1], nrow = 4)
dev.off()

timepoints <- c(-100, 0, 100)    
frac.table = matrix(
  c(36.2, 
    36.2,
    36.2),
  ncol=length(timepoints))
frac.table <- 0.5*frac.table
parents = c(0)
fish = createFishObject(frac.table,parents,timepoints=timepoints,
                        col = c("red"),
                        clone.labels = c("Donor-SH2B3	480_487del"))
fish = layoutClones(fish)

png("hct_753399 donor.png", width = 660, height = 480, bg = "transparent")
fishPlot(fish,shape="polygon",title.btm= hct_id,
         cex.title=1, vlines=NULL,
         ramp.angle = 1,
         pad.left = 0.05,
         vlab=NULL, 
         bg.type = "solid",
         bg.col = "transparent"
)
drawLegend(fish,cex=1.0,xpos=fish@timepoints[1])
dev.off()
#################

hct_id <- "hct_753525"

timepoints <- c(-114, -16, 0, 42, 84, 346, 514, 684)    

# Donor 0	DNMT3A	P904L	0.116
# -114	DNMT3A	E578X	0.459
# -114	NRAS	G12D	0.221
# -114	IDH1	R132C	0.221
# -16	DNMT3A	R882H	0.094     A
# -16	DNMT3A	E578X	0.051     A
# 42	DNMT3A	P904L	0.179
# 42	KRAS	T50I	0.199
# 42	RUNX1	I177V	0.116
# 84	DNMT3A	P904L	0.181     A
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
  c(0, 45.9, 22.1, 0.0001, 22.1,        0, 0, #-114
    0, 9.5, 0.0001, 9.4, 0.0001,        0, 0, #-16 ###### Change 5.1 to 9.5
    
    0.0001, 0, 0, 0, 0,                      0.0001, 0.0001, # 0
    
    17.9, 0, 0, 0, 0,                      19.9, 11.6, #42
    18.1, 0, 0, 0, 0,                      0.0001, 0.0001, #84
    15.6, 0, 0, 0, 0,                      15.5, 10.4, #346
    17, 0, 0, 0, 0,                      15, 10, #514
    12.9, 0, 0, 0, 0,                      13.9, 9.7),
  ncol=length(timepoints))
# DNMT3A	IDH1	<0.001	<0.001	Co-occurrence
# KRAS	NRAS	<0.001	<0.001	Co-occurrence
# DNMT3A	RUNX1	<0.001	0.002	Mutual exclusivity
# IDH1	NRAS	0.008	0.019	Co-occurrence
# RUNX1	NRAS	0.057	0.095	Co-occurrence
# IDH1	RUNX1	0.057	0.095	Co-occurrence
# DNMT3A	KRAS	0.077	0.110	Mutual exclusivity
# IDH1	KRAS	0.103	0.129	Mutual exclusivity
# DNMT3A	NRAS	0.241	0.267	Co-occurrence
# KRAS	RUNX1	0.345	0.345	Co-occurrence
parents = c(0,0,2, 2, 3, 0, 6)
fish = createFishObject(frac.table,parents,timepoints=timepoints,
                        col = c("red", "slategray1", 
                                "gold", "darkblue", 
                                "blueviolet", "#70e000", 
                                "#ffff3f"),
                        clone.labels = c("DNMT3A P904L", "DNMT3A E578X", 
                                         "NRAS G12D", "DNMT3A R882H", 
                                         "IDH1 R132C", "KRAS T50I", 
                                         "RUNX1 I177V"))
fish = layoutClones(fish)

png("hct_753525 recipient.png", width = 660, height = 480, bg = "transparent")
fishPlot(fish,shape="spline",title.btm= hct_id,
         cex.title=1, vlines=c(timepoints), col.vline = "gainsboro",
         vlab=c(timepoints), 
         bg.type = "solid",
         bg.col = "transparent"
)
drawLegend(fish,cex=1.0,xpos=fish@timepoints[1], nrow = 4)
dev.off()

timepoints <- c(-100, 0, 100)    
frac.table = matrix(
  c(3.8, 
    3.8,
    3.8),
  ncol=length(timepoints))
parents = c(0)
fish = createFishObject(frac.table,parents,timepoints=timepoints,
                        col = c("red"),
                        clone.labels = c("DNMT3A P904L"))
fish = layoutClones(fish)

png("hct_753525 donor.png", width = 660, height = 480, bg = "transparent")
fishPlot(fish,shape="polygon",title.btm= hct_id,
         cex.title=1, vlines=NULL,
         ramp.angle = 1,
         pad.left = 0.05,
         vlab=NULL, 
         bg.type = "solid",
         bg.col = "transparent"
)
drawLegend(fish,cex=1.0,xpos=fish@timepoints[1])
dev.off()
#################

hct_id <- "hct_769464"

timepoints <- c(-114, -23, 0, 98, 339, 538, 696)    

# Donor 0	DNMT3A	R882C	0.038
# -114	U2AF1	Q157R	0.322
# -23	NA	NA	NA            A
# 98	TP53	V104M	0.038     A
# 98	DNMT3A	R882C	0.1     A
# 339	DNMT3A	R882C	0.151
# 538	DNMT3A	R882C	0.165
# 696	DNMT3A	R882C	0.147

frac.table = matrix(
  c(32.2, 0, 0,  #-114
    0.0001, 0.0001, 0.0001,  #-23
    0, 0, 0,  #0
    0, 3.8, 10,  #98
    0, 0.0001, 15.1,  #339
    0, 0.0001, 16.5,  #538
    0, 0.0001, 14.7),
  ncol=length(timepoints))
# DNMT3A	TP53	<0.001	<0.001	Mutual exclusivity
parents = c(0,0, 0)
fish = createFishObject(frac.table,parents,timepoints=timepoints,
                        col = c("orange", "blue", 
                                "red"),
                        clone.labels = c("U2AF1", "TP53 V104M", 
                                         "DNMT3A R882C"),
                        fix.missing.clones=TRUE)
fish = layoutClones(fish)

png("hct_769464 recipient.png", width = 660, height = 480, bg = "transparent")
fishPlot(fish,shape="spline",title.btm= hct_id,
         cex.title=1, vlines=c(timepoints), col.vline = "gainsboro",
         vlab=c(timepoints), 
         bg.type = "solid",
         bg.col = "transparent"
)
drawLegend(fish,cex=1.0,xpos=fish@timepoints[1], nrow = 4)
dev.off()

timepoints <- c(-100, 0, 100)    
frac.table = matrix(
  c(3.8, 
    3.8,
    3.8),
  ncol=length(timepoints))
parents = c(0)
fish = createFishObject(frac.table,parents,timepoints=timepoints,
                        col = c("red"),
                        clone.labels = c("DNMT3A R882C"))
fish = layoutClones(fish)

png("hct_769464 donor.png", width = 660, height = 480, bg = "transparent")
fishPlot(fish,shape="polygon",title.btm= hct_id,
         cex.title=1, vlines=NULL,
         ramp.angle = 1,
         pad.left = 0.05,
         vlab=NULL, 
         bg.type = "solid",
         bg.col = "transparent"
)
drawLegend(fish,cex=1.0,xpos=fish@timepoints[1])
dev.off()
#################







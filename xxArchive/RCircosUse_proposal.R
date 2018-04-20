setwd("C:/Users/Kirk/OneDrive/Documents/Dissertation")
library(RCircos)
library(data.table)
CytoBandIdeogram.Data <- read.table("CytoBandIdeogramExample.txt",header=T)
#Rescalling Ideogram data to fit C Elegans

MyCytoBand.Data<-CytoBandIdeogram.Data

chr4<-MyCytoBand.Data[MyCytoBand.Data$Chromosome=="chr4",]
new.chr0<-chr4
new.chr0$Chromosome<-"chr0"

chr1<-MyCytoBand.Data[MyCytoBand.Data$Chromosome=="chr1",]
new.chr4<-chr1
new.chr4$Chromosome<-"chr4"

MyCytoBand.Data<-MyCytoBand.Data[!(MyCytoBand.Data$Chromosome=="chr4"),]
MyCytoBand.Data<-rbind(new.chr0,MyCytoBand.Data,new.chr4)
MyCytoBand.Data<-MyCytoBand.Data[order(MyCytoBand.Data$Chromosome),]
MyCytoBand.Data$ChromStart<-MyCytoBand.Data$ChromStart/10
MyCytoBand.Data$ChromEnd<-MyCytoBand.Data$ChromEnd/10
MyCytoBand.Data$ChromEnd[148]<-14500000
MyCytoBand.Data$ChromEnd[193]<-18000000
MyCytoBand.Data$ChromEnd[227]<-20500000

tmp<-MyCytoBand.Data
tmp$Chromosome<-paste("SNP", tmp$Chromosome, sep=".")
MyCytoBand.Data$Chromosome<-paste("Transcript", MyCytoBand.Data$Chromosome, sep=".")
MyCytoBand.Data<-rbind(MyCytoBand.Data,tmp)
MyCytoBand.Data$Chromosome[MyCytoBand.Data$Chromosome=="Transcript.chr6"]<-"Transcript.chr10"
MyCytoBand.Data$Chromosome[MyCytoBand.Data$Chromosome=="SNP.chr6"]<-"SNP.chr10"



# Creating Link Data
output1<-read.csv("Proposal/iForm Output/iForm Output cQTL/iform.output.69.csv")
output2<-read.csv("Proposal/iForm Output/iform.output.59.csv")
output3<-read.csv("Proposal/iForm Output/iform.output.6.csv")
output4<-read.csv("Proposal/iForm Output/iform.output.77.csv")
output5<-read.csv("Proposal/iForm Output/iform.output.38.csv")
output6<-read.csv("Proposal/iForm Output/iform.output.15.csv")

# A_12_P118280 # chr 1  output.69 cQTL
# A_12_P107874 # chr 2  output.59
# A_12_P102922 # chr 3  output.6 need to find one with interactions
# A_12_P103547 # chr 4  output.77
# A_12_P106909 # chr 5  output.38
# A_12_P107771 # chr X  output.15



Link.Chr1<-Link.Data.Create(output1)
Link.Chr2<-Link.Data.Create(output2)
Link.Chr3<-Link.Data.Create(output3)
Link.Chr4<-Link.Data.Create(output4)
Link.Chr5<-Link.Data.Create(output5)
Link.Chr6<-Link.Data.Create(output6)

 ## Creating the RCircos Plot ###################

cyto.info <- MyCytoBand.Data
num.inside <- 5
num.outside <- 0
RCircos.Set.Core.Components(cyto.info, NULL, num.inside, num.outside)

jpegFileName <- sprintf("iForm_Network_Chr1_5.jpeg")
jpeg(jpegFileName, width = 2000, height = 2000, quality = 1000)
set.seed(1)
RCircos.Set.Plot.Area()
My.RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
My.RCircos.Link.Plot(Link.Chr1, track.num, TRUE)
My.RCircos.Link.Plot(Link.Chr2, track.num, TRUE)
My.RCircos.Link.Plot(Link.Chr3, track.num, TRUE)
My.RCircos.Link.Plot(Link.Chr4, track.num, TRUE)
My.RCircos.Link.Plot(Link.Chr5, track.num, TRUE)
My.RCircos.Link.Plot(Link.Chr6, track.num, TRUE)
dev.off()










y<-gene.expr[,"A_12_P103290"]

dat1<-data.frame(y,x)

fit1<-iForm.BIC.base2(dat1,"y")

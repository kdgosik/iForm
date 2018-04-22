library(RCircos)
My.RCircos.Link.Plot <- function (link.data, track.num, by.chromosome = FALSE) 
{
  RCircos.Pos <- RCircos.Get.Plot.Positions()
  RCircos.Pos <- RCircos.Pos[c(1:(dim(RCircos.Pos)[1]/2),(dim(RCircos.Pos)[1]):((dim(RCircos.Pos)[1]/2)+1)),]
  tmp.n <- round(dim(RCircos.Pos)[1]/16*7.5)
  RCircos.Pos <- rbind(RCircos.Pos[(tmp.n+1):dim(RCircos.Pos)[1],], RCircos.Pos[1:tmp.n,])
  RCircos.Par <- RCircos.Get.Plot.Parameters()
  link.data <- RCircos.Validate.Genomic.Data(link.data, plot.type = "link")
  one.track <- RCircos.Par$track.height + RCircos.Par$track.padding
  start <- RCircos.Par$track.in.start - (track.num - 1) * one.track
  base.positions <- RCircos.Pos * start
  data.points <- matrix(rep(0, nrow(link.data) * 2), ncol = 2)
  for (a.link in 1:nrow(link.data)) {
    data.points[a.link, 1] <- RCircos.Data.Point(link.data[a.link, 
                                                           1], link.data[a.link, 2])
    data.points[a.link, 2] <- RCircos.Data.Point(link.data[a.link, 
                                                           4], link.data[a.link, 5])
    if (data.points[a.link, 1] == 0 || data.points[a.link, 
                                                   2] == 0) {
      print("Error in chromosome locations ...")
      break
    }
  }
  points <- NULL
  for (a.link in 1:nrow(data.points)) {
    point.one <- data.points[a.link, 1]
    point.two <- data.points[a.link, 2]
    points <- c(points, point.one, point.two)
    if (point.one > point.two) {
      point.one <- data.points[a.link, 2]
      point.two <- data.points[a.link, 1]
    }
    P0 <- as.numeric(base.positions[point.one, ])
    P2 <- as.numeric(base.positions[point.two, ])
    links <- RCircos.Link.Line(P0, P2)
    
    if(link.data[a.link, 4] == "SNP.chr1"){
      lines(links$pos.x, links$pos.y, type = "l", col = "black", lwd=1)
    }
    if(link.data[a.link, 4] == "SNP.chr2"){
      lines(links$pos.x, links$pos.y, type = "l", col = "red", lwd=1)
    }
    if(link.data[a.link, 4] == "SNP.chr3"){
      lines(links$pos.x, links$pos.y, type = "l", col = "darkgreen", lwd=1)
    }
    if(link.data[a.link, 4] == "SNP.chr4"){
      lines(links$pos.x, links$pos.y, type = "l", col = "darkblue", lwd=1)
    }
    if(link.data[a.link, 4] == "SNP.chr5"){
      lines(links$pos.x, links$pos.y, type = "l", col = "purple", lwd=1)
    }
    if(link.data[a.link, 4] == "SNP.chr10"){
      lines(links$pos.x, links$pos.y, type = "l", col = "chocolate", lwd=1)
    }
    
  }
  
  text(-1.75, 0, label = "Transcripts", srt = 90, cex = 5)
  text(1.75, 0, label = "SNPs", srt = -90, cex = 5)
}






My.RCircos.Ribbon.Plot<-function (ribbon.data, track.num, by.chromosome = FALSE, twist = FALSE) 
{
  RCircos.Pos <- RCircos.Get.Plot.Positions()
  RCircos.Pos <- RCircos.Pos[c(1:(dim(RCircos.Pos)[1]/2),(dim(RCircos.Pos)[1]):((dim(RCircos.Pos)[1]/2)+1)),]
  tmp.n <- round(dim(RCircos.Pos)[1]/16*7.5)
  RCircos.Pos <- rbind(RCircos.Pos[(tmp.n+1):dim(RCircos.Pos)[1],], RCircos.Pos[1:tmp.n,])
  RCircos.Par <- RCircos.Get.Plot.Parameters()
  ribbon.data <- RCircos.Validate.Genomic.Data(ribbon.data, 
                                               plot.type = "link")
  one.track <- RCircos.Par$track.height + RCircos.Par$track.padding
  track.out <- RCircos.Par$track.in.start - (track.num - 1) * 
    one.track
  base.positions <- RCircos.Pos * track.out
  data.points <- matrix(rep(0, nrow(ribbon.data) * 4), ncol = 4)
  for (a.link in 1:nrow(ribbon.data)) {
    data.points[a.link, 1] <- RCircos.Data.Point(ribbon.data[a.link, 
                                                             1], ribbon.data[a.link, 2])
    data.points[a.link, 2] <- RCircos.Data.Point(ribbon.data[a.link, 
                                                             1], ribbon.data[a.link, 3])
    data.points[a.link, 3] <- RCircos.Data.Point(ribbon.data[a.link, 
                                                             4], ribbon.data[a.link, 5])
    data.points[a.link, 4] <- RCircos.Data.Point(ribbon.data[a.link, 
                                                             4], ribbon.data[a.link, 6])
    if (data.points[a.link, 1] == 0 || data.points[a.link, 
                                                   2] == 0 || data.points[a.link, 3] == 0 || data.points[a.link, 
                                                                                                         4] == 0) {
      stop("Error in chromosome locations ...")
    }
  }
  ribbon.colors <- RCircos.Get.Link.Colors(ribbon.data, by.chromosome)
  for (a.ribbon in 1:nrow(ribbon.data)) {
    start.one <- data.points[a.ribbon, 1]
    end.one <- data.points[a.ribbon, 2]
    if (twist == FALSE) {
      start.two <- data.points[a.ribbon, 3]
      end.two <- data.points[a.ribbon, 4]
    }
    else {
      start.two <- data.points[a.ribbon, 4]
      end.two <- data.points[a.ribbon, 3]
    }
    P0 <- as.numeric(base.positions[end.one, ])
    P2 <- as.numeric(base.positions[start.two, ])
    line.one <- RCircos.Link.Line(P0, P2)
    P0 <- as.numeric(base.positions[end.two, ])
    P2 <- as.numeric(base.positions[start.one, ])
    line.two <- RCircos.Link.Line(P0, P2)
    polygon.x <- c(base.positions[start.one:end.one, 1], 
                   line.one$pos.x, base.positions[start.two:end.two, 
                                                  1], line.two$pos.x)
    polygon.y <- c(base.positions[start.one:end.one, 2], 
                   line.one$pos.y, base.positions[start.two:end.two, 
                                                  2], line.two$pos.y)
    polygon(polygon.x, polygon.y, border = NA, col = ribbon.colors[a.ribbon])
  }
  text(-1.75, 0, label = "Transcripts", srt = 90, cex = 3)
  text(1.75, 0, label = "SNPs", srt = -90, cex = 3)
}




My.RCircos.Chromosome.Ideogram.Plot <- function () 
{
  RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
  RCircos.Pos <- RCircos.Get.Plot.Positions()
  RCircos.Pos <- RCircos.Pos[c(1:(dim(RCircos.Pos)[1]/2),(dim(RCircos.Pos)[1]):((dim(RCircos.Pos)[1]/2)+1)),]
  tmp.n <- round(dim(RCircos.Pos)[1]/16*7.5)
  RCircos.Pos <- rbind(RCircos.Pos[(tmp.n+1):dim(RCircos.Pos)[1],], RCircos.Pos[1:tmp.n,])
  RCircos.Par <- RCircos.Get.Plot.Parameters()
  right.side <- nrow(RCircos.Pos)/2
  outer.location <- RCircos.Par$chr.ideog.pos + RCircos.Par$chrom.width
  inner.location <- RCircos.Par$chr.ideog.pos
  chroms <- unique(RCircos.Cyto$Chromosome)
  for (a.chr in 1:length(chroms)) {
    the.chr <- RCircos.Cyto[RCircos.Cyto$Chromosome == chroms[a.chr], ]
    start <- the.chr$Location[1] - the.chr$Unit[1] + 1
    end <- the.chr$Location[nrow(the.chr)]
    mid <- round((end - start + 1)/2, digits = 0) + start
    
    if(gsub("[^0-9]","",chroms[a.chr])=="0") chr.color<-NA
    if(gsub("[^0-9]","",chroms[a.chr])=="1") chr.color<-"black"
    if(gsub("[^0-9]","",chroms[a.chr])=="2") chr.color<-"red"
    if(gsub("[^0-9]","",chroms[a.chr])=="3") chr.color<-"darkgreen"
    if(gsub("[^0-9]","",chroms[a.chr])=="4") chr.color<-"darkblue"
    if(gsub("[^0-9]","",chroms[a.chr])=="5") chr.color<-"purple"
    if(gsub("[^0-9]","",chroms[a.chr])=="10") chr.color<-"chocolate"    
    
    pos.x <- c(RCircos.Pos[start:end, 1] * outer.location, 
               RCircos.Pos[end:start, 1] * inner.location)
    pos.y <- c(RCircos.Pos[start:end, 2] * outer.location, 
               RCircos.Pos[end:start, 2] * inner.location)
    polygon(pos.x, pos.y, col = chr.color, border = NA)
    chr.name <- gsub("[^0-9]","",chroms[a.chr])
    text(RCircos.Pos[mid, 1] * (RCircos.Par$chr.name.pos-.1), 
         RCircos.Pos[mid, 2] * (RCircos.Par$chr.name.pos-.1), label = as.roman(chr.name), 
         srt = RCircos.Pos$degree[mid], cex = 4)
    
  }
  
}



Link.Data.Create<-function(output){
  
  transcript.location<-output[dim(output)[1],]
  output<-output[-dim(output)[1],]
  
  # chrX  4504351 - 4504292
  
  chr.trans<-substr(transcript.location[,3],4,10)
  
  if(chr.trans=="X") trans.Chromosome<-paste0("Transcript.chr",10)
  if(chr.trans=="I") trans.Chromosome<-paste0("Transcript.chr",1)
  if(chr.trans=="II") trans.Chromosome<-paste0("Transcript.chr",2)
  if(chr.trans=="III") trans.Chromosome<-paste0("Transcript.chr",3)
  if(chr.trans=="IV") trans.Chromosome<-paste0("Transcript.chr",4)
  if(chr.trans=="V") trans.Chromosome<-paste0("Transcript.chr",5)
  
  trans.chromStart<-round(transcript.location[,4],-5)
  trans.chromEnd<-trans.chromStart+100
  
  
  trans.Chromosome.1<-sapply(output[,1][!{output[,1]%like%"[.]"}],function(x){
    if(substr(x,2,2)=="_") {"SNP.chr10"}
    else{paste0("SNP.chr",substr(x,2,2))}
  })
  
  trans.chromStart.1<-round(as.numeric(sapply(output[,1][!{output[,1]%like%"[.]"}],function(x){
    if(substr(x,2,2)=="_") {gsub("[^0-9]","",x)}
    else{substr(x,4,15)}
  })),-5)
  
  trans.chromEnd.1<-trans.chromStart.1+100
  
  
  To<-matrix(rep(c(trans.Chromosome,trans.chromStart,trans.chromEnd),length(trans.chromStart.1)),ncol=3,byrow=T)
  
  My.Link<-cbind(To,trans.Chromosome.1,trans.chromStart.1,trans.chromEnd.1)
  My.Link<-data.frame(My.Link)
  colnames(My.Link)<-c("Chromosome","chromStart","chromEnd","Chromosome.1","chromStart.1","chromEnd.1")
  
  My.Link$chromStart<-as.numeric(paste(My.Link$chromStart))
  My.Link$chromStart.1<-as.numeric(paste(My.Link$chromStart.1))
  My.Link$chromEnd<-as.numeric(paste(My.Link$chromEnd))
  My.Link$chromEnd.1<-as.numeric(paste(My.Link$chromEnd.1))
  
  
  
  # Getting Interaction SNPs
  
  tmp1<-output[,1][{output[,1]%like%"[.]"}]
  tmp2<-strsplit(as.character(tmp1),".",fixed=T)
  tmp3<-do.call(rbind,tmp2)
  
  Chromosome<-sapply(tmp3[,1],function(x){
    if(substr(x,2,2)=="_") {"SNP.chr10"}
    else{paste0("SNP.chr",substr(x,2,2))}
  })
  
  chromStart<-round(as.numeric(sapply(tmp3[,1],function(x){
    if(substr(x,2,2)=="_") {gsub("[^0-9]","",x)}
    else{substr(x,4,15)}
  })),-5)
  
  chromEnd<-chromStart+100
  
  Chromosome.1<-sapply(tmp3[,2],function(x){
    if(substr(x,2,2)=="_") {"SNP.chr10"}
    else{paste0("SNP.chr",substr(x,2,2))}
  })
  
  chromStart.1<-round(as.numeric(sapply(tmp3[,2],function(x){
    if(substr(x,2,2)=="_") {gsub("[^0-9]","",x)}
    else{substr(x,4,15)}
  })),-5)
  
  chromEnd.1<-chromStart.1+100
  
  interaction.data<-data.frame(Chromosome,chromStart,chromEnd,Chromosome.1,chromStart.1,chromEnd.1)

  interaction.data$chromStart<-as.numeric(paste(interaction.data$chromStart))
  interaction.data$chromStart.1<-as.numeric(paste(interaction.data$chromStart.1))
  interaction.data$chromEnd<-as.numeric(paste(interaction.data$chromEnd))
  interaction.data$chromEnd.1<-as.numeric(paste(interaction.data$chromEnd.1))
  
  rbind(My.Link,interaction.data)
  
}

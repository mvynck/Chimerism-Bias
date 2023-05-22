# source functinos
# load data
setwd("/Users/mvynck/Library/CloudStorage/GoogleDrive-matthijs.vynck@gmail.com/Mijn Drive/Werk-AZSintJan/onderzoek_chimerisme/Chimerisme_biascorrection")
source("Code/biasFuns.R")
bias <- unname(as.vector(read.csv("Data/bias.csv")))[[1]]
samples <- list.files("./data/invitromixtures/")

# load libraries
library(ggplot2)
library(cowplot)


# example s

recipient <- read.table("Data/invitromixtures/Aref.txt", sep="\t")
donor <- read.table("Data/invitromixtures/Bref.txt", sep="\t")
Aid <- c(1, 2, 3, 4, 5, 7, 9, 11, 12, 13, 15, 16, 17, 19, 20, 21)
resA <- matrix(0, nrow = length(Aid), ncol = 24)
for(samp in 1:length(Aid)){
  followup <- read.table(paste0("./data/invitromixtures/",samples[Aid][samp]), sep="\t")
  resA[samp,] <- followUp(followup, recipient, donor, bias = bias)
}
recipient <- read.table("Data/invitromixtures/Bref.txt", sep="\t")
donor <- read.table("Data/invitromixtures/Aref.txt", sep="\t")
Bid <- c(1, 2, 3, 4, 5, 7, 9, 11, 12, 13, 15, 16, 17, 19, 20, 21)
resB <- matrix(0, nrow = length(Aid), ncol = 24)
for(samp in 1:length(Bid)){
  followup <- read.table(paste0("./data/invitromixtures/",samples[Bid][samp]), sep="\t")
  resB[samp,] <- followUp(followup, recipient, donor, bias = bias)
}

recipient <- read.table("Data/invitromixtures/Cref.txt", sep="\t")
donor <- read.table("Data/invitromixtures/Dref.txt", sep="\t")
Cid <- c(6, 8, 10, 14, 18, 22)
resC <- matrix(0, nrow = length(Cid), ncol = 24)
for(samp in 1:length(Cid)){
  followup <- read.table(paste0("./data/invitromixtures/",samples[Cid][samp]), sep="\t")
  resC[samp,] <- followUp(followup, recipient, donor, bias = bias)
}
recipient <- read.table("Data/invitromixtures/Dref.txt", sep="\t")
donor <- read.table("Data/invitromixtures/Cref.txt", sep="\t")
Did <- c(6, 8, 10, 14, 18, 22)
resD <- matrix(0, nrow = length(Cid), ncol = 24)
for(samp in 1:length(Did)){
  followup <- read.table(paste0("./data/invitromixtures/",samples[Did][samp]), sep="\t")
  resD[samp,] <- followUp(followup, recipient, donor, bias = bias)
}

res <- rbind(resA, resB, resC, resD)
res <- as.data.frame(res)
res$expected <- c(0.1, 0.1, 0.5, 0.5, 0.1, 0.5, 0, 10, 10, 10, 1, 1, 1, 5, 5, 5,
                  99.9, 99.9, 99.5, 99.5, 99.9, 99.5, 100, 90, 90, 90, 99, 99, 99, 95, 95, 95,
                  0.1, 0.5, 0, 10, 1, 5,
                  99.9, 99.5, 100, 90, 99, 95)/100
colnames(res) <- c("meanInf", "medianInf", "sdInf", "nInf",
                   "meanInfB", "medianInfB", "sdInfB", "nInfB",
                   "meanPotInf", "medianPotInf", "sdPotInf", "nPotInf",
                   "meanPotInfB", "medianPotInfB", "sdPotInfB", "nPotInfB",
                   "meanAll", "medianAll", "sdAll", "nAll",
                   "meanAllB", "medianAllB", "sdAllB", "nAllB", "Expected")
plot(res$Expected,res$meanInf-res$Expected)
plot(res$Expected,res$meanInfB-res$Expected)
plot(res$Expected,res$meanPotInf-res$Expected)
plot(res$Expected,res$meanPotInfB-res$Expected)






# nice plot for CD

recipient <- read.table("Data/invitromixtures/Cref.txt", sep="\t")
donor <- read.table("Data/invitromixtures/Dref.txt", sep="\t")
Cid <- c(6, 8, 10, 14, 18, 22)
resC <- matrix(0, nrow = length(Cid)*10, ncol = 2)
for(samp in 1:length(Cid)){
  followup <- read.table(paste0("./data/invitromixtures/",samples[Cid][samp]), sep="\t")
  resC[c((10*(samp-1)+1):(10*samp)),c(1:2)] <- followUp(followup, recipient, donor, bias = bias, reportMarker = TRUE, type = "inf")
}
resC <- as.data.frame(resC)
resC$Expected <- rep(c(0.1, 0.5, 0, 10, 1, 5), each = 10)/100
colnames(resC) <- c("InfM", "InfB", "Expected")
par(mfrow=c(1,2))
plot(resC$Expected, resC$InfM-resC$Expected, ylim=c(-0.018, 0.018))
plot(resC$Expected, resC$InfB-resC$Expected, ylim=c(-0.018, 0.018))
plotdfC <- data.frame(Difference=100*c(resC$InfM-resC$Expected, resC$InfB-resC$Expected),
                     Condition=factor(c(rep("No correction", nrow(resC)), rep("Bias-corrected", nrow(resC)))),
                     Expected=paste0("Expected: ",factor(100*rep(resC$Expected, 2)), "% HC"),
                     Group=paste0("C",c(1:nrow(resC),1:nrow(resC))),
                     Type=rep(c("Heterozygous", "Homozygous", "Heterozygous", "Heterozygous", "Heterozygous", "Homozygous","Heterozygous","Heterozygous","Heterozygous","Heterozygous"), 2*length(Cid)))
plotdfC$Condition <- relevel(plotdfC$Condition, "No correction")
ggplot(plotdfC, aes(x=Condition, y=Difference, color=Type, group=Group))+
  geom_point()+
  geom_line()+
  facet_wrap(~Expected)+
  theme_minimal()+
  ylab("Deviation from expected %HC (%)")


recipient <- read.table("Data/invitromixtures/Dref.txt", sep="\t")
donor <- read.table("Data/invitromixtures/Cref.txt", sep="\t")
Did <- c(6, 8, 10, 14, 18, 22)
resD <- matrix(0, nrow = length(Did)*4, ncol = 2)
for(samp in 1:length(Did)){
  followup <- read.table(paste0("./data/invitromixtures/",samples[Did][samp]), sep="\t")
  resD[c((4*(samp-1)+1):(4*samp)),c(1:2)] <- followUp(followup, recipient, donor, bias = bias, reportMarker = TRUE, type = "inf")
}
resD <- as.data.frame(resD)
resD$Expected <- rep(100-c(0.1, 0.5, 0, 10, 1, 5), each = 4)/100
colnames(resD) <- c("InfM", "InfB", "Expected")
par(mfrow=c(1,2))
plot(resD$Expected, resD$InfM-resD$Expected, ylim=c(-0.11, 0.11))
plot(resD$Expected, resD$InfB-resD$Expected, ylim=c(-0.11, 0.11))
plotdfD <- data.frame(Difference=100*c(resD$InfM-resD$Expected, resD$InfB-resD$Expected),
                      Condition=factor(c(rep("No correction", nrow(resD)), rep("Bias-corrected", nrow(resD)))),
                      Expected=paste0("Expected: ",factor(100*rep(resD$Expected, 2)), "% HC"),
                      Group=paste0("D",c(1:nrow(resD),1:nrow(resD))),
                      Type=rep(c("Heterozygous", "Homozygous", "Homozygous", "Heterozygous"), 2*length(Did)))
plotdfD$Condition <- relevel(plotdfD$Condition, "No correction")
ggplot(plotdfD, aes(x=Condition, y=Difference, color=Type, group=Group))+
  geom_point()+
  geom_line()+
  facet_wrap(~Expected)+
  theme_minimal()+
  ylab("Deviation from expected %HC (%)")

plotdfCD <- rbind(plotdfC, plotdfD)
plotdfCD$Expected <- factor(plotdfCD$Expected, levels=levels(factor(plotdfCD$Expected))[c(3, 1, 2, 4, 7, 5, 8, 9, 12, 10, 11, 6)])
plotdfCD$Marker <- "Type I"
ggplot(plotdfCD, aes(x=Condition, y=Difference, color=Type, group=Group))+
  geom_point()+
  geom_line()+
  facet_wrap(~Expected, nrow = 2)+
  theme_minimal()+
  xlab("")+
  ylab("Deviation from expected %HC (%)")


recipient <- read.table("Data/invitromixtures/Cref.txt", sep="\t")
donor <- read.table("Data/invitromixtures/Dref.txt", sep="\t")
resCB <- matrix(0, nrow = length(Cid)*2, ncol = 2)
for(samp in 1:length(Cid)){
  followup <- read.table(paste0("./data/invitromixtures/",samples[Cid][samp]), sep="\t")
  resCB[c((2*(samp-1)+1):(2*samp)),c(1:2)] <- followUp(followup, recipient, donor, bias = bias, reportMarker = TRUE, type = "potinf")
}
resCB <- as.data.frame(resCB)
resCB$Expected <- rep(c(0.1, 0.5, 0, 10, 1, 5), each = 2)/100
colnames(resCB) <- c("InfM", "InfB", "Expected")
par(mfrow=c(1,2))
plot(resCB$Expected, resCB$InfM-resCB$Expected, ylim=c(-0.11, 0.11))
plot(resCB$Expected, resCB$InfB-resCB$Expected, ylim=c(-0.11, 0.11))
plotdfCB <- data.frame(Difference=100*c(resCB$InfM-resCB$Expected, resCB$InfB-resCB$Expected),
                     Condition=factor(c(rep("No correction", nrow(resCB)), rep("Bias-corrected", nrow(resCB)))),
                     Expected=paste0("Expected: ",factor(100*rep(resCB$Expected, 2)), "% HC"),
                     Group=c(1:nrow(resCB),1:nrow(resCB)),
                     Type="Heterozygous")
plotdfCB$Condition <- relevel(plotdfCB$Condition, "No correction")
ggplot(plotdfCB, aes(x=Condition, y=Difference, color=Type, group=Group))+
  geom_point()+
  geom_line()+
  facet_wrap(~Expected, nrow = 2)+
  theme_minimal()+
  xlab("")+
  ylab("Deviation from expected %HC (%)")

recipient <- read.table("Data/invitromixtures/Dref.txt", sep="\t")
donor <- read.table("Data/invitromixtures/Cref.txt", sep="\t")
resDB <- matrix(0, nrow = length(Did)*8, ncol = 2)
for(samp in 1:length(Did)){
  followup <- read.table(paste0("./data/invitromixtures/",samples[Did][samp]), sep="\t")
  resDB[c((8*(samp-1)+1):(8*samp)),c(1:2)] <- followUp(followup, recipient, donor, bias = bias, reportMarker = TRUE, type = "potinf")
}
resDB <- as.data.frame(resDB)
resDB$Expected <- rep(100-c(0.1, 0.5, 0, 10, 1, 5), each = 8)/100
colnames(resDB) <- c("InfM", "InfB", "Expected")
par(mfrow=c(1,2))
plot(resDB$Expected, resDB$InfM-resDB$Expected, ylim=c(-0.11, 0.11))
plot(resDB$Expected, resDB$InfB-resDB$Expected, ylim=c(-0.11, 0.11))
plotdfDB <- data.frame(Difference=100*c(resDB$InfM-resDB$Expected, resDB$InfB-resDB$Expected),
                     Condition=factor(c(rep("No correction", nrow(resDB)), rep("Bias-corrected", nrow(resDB)))),
                     Expected=paste0("Expected: ", factor(100*rep(resDB$Expected, 2)),"% HC"),
                     Group=c(1:nrow(resDB),1:nrow(resDB)),
                     Type='Heterozygous')
plotdfDB$Condition <- relevel(plotdfDB$Condition, "No correction")
ggplot(plotdfDB, aes(x=Condition, y=Difference, color=Expected, group=Group))+
  geom_point()+
  geom_line()+
  theme_minimal()+
  ylab("Deviation from expected %HC (%)")

plotdfCDB <- rbind(plotdfCB, plotdfDB)
plotdfCDB$Expected <- factor(plotdfCDB$Expected, levels=levels(factor(plotdfCDB$Expected)))
plotdfCDB$Marker <- "Type II"
ggplot(plotdfCDB, aes(x=Condition, y=Difference, color=Type, group=Group))+
  geom_point()+
  geom_line()+
  facet_wrap(~Expected, nrow = 2)+
  theme_minimal()+
  xlab("")+
  ylab("Deviation from expected %HC (%)")

plotall <- rbind(plotdfCD, plotdfCDB)
colnames(plotall)[5:6] <- c("Constellation", "Type")
for(sampleid in 1:length(unique(plotall$Expected))){
  explev <- levels(plotall$Expected)[sampleid]
  nocorr <- mean(plotall$Difference[intersect(which(plotall$Expected == explev),which(plotall$Condition == "No correction"))])
  corr <- mean(plotall$Difference[intersect(which(plotall$Expected == explev),which(plotall$Condition == "Bias-corrected"))])
  plotall <- rbind(plotall, c(nocorr, "No correction", levels(plotall$Expected)[sampleid], paste0("mean", sampleid), "Not applicable", "Combined"))
  plotall <- rbind(plotall, c(corr, "Bias-corrected", levels(plotall$Expected)[sampleid], paste0("mean", sampleid), "Not applicable", "Combined"))
  plotall$Difference <- as.numeric(plotall$Difference)
}
ggplot(plotall, aes(x=Condition, y=(Difference), color=Type, shape=Constellation, group=Group))+
  geom_point()+
  geom_line()+
  facet_wrap(~Expected, nrow = 4)+
  theme_minimal()+
  theme(legend.position = "top")+
  xlab("")+coord_flip()+
  ylab("Deviation from expected %HC (%)")





# nice plot for AB

recipient <- read.table("Data/invitromixtures/Aref.txt", sep="\t")
donor <- read.table("Data/invitromixtures/Bref.txt", sep="\t")
Cid <- c(1:22)[-c(6, 8, 10, 14, 18, 22)]
resC <- matrix(0, nrow = length(Cid)*11, ncol = 2)
for(samp in 1:length(Cid)){
  followup <- read.table(paste0("./data/invitromixtures/",samples[Cid][samp]), sep="\t")
  resC[c((11*(samp-1)+1):(11*samp)),c(1:2)] <- followUp(followup, recipient, donor, bias = bias, reportMarker = TRUE, type = "inf")
}
resC <- as.data.frame(resC)
resC$Expected <- rep(c(0.1,0.1, 0.5, 0.5, 0.1, 0.5, 0, 10, 10, 10, 1,1,1, 5, 5, 5), each = 11)/100
resC$Replicate <- rep(samples[Cid], each = 11)
colnames(resC) <- c("InfM", "InfB", "Expected", "Replicate")
par(mfrow=c(1,2))
plot(resC$Expected, resC$InfM-resC$Expected, ylim=c(-0.028, 0.028))
plot(resC$Expected, resC$InfB-resC$Expected, ylim=c(-0.028, 0.028))
plotdfC <- data.frame(Difference=100*c(resC$InfM-resC$Expected, resC$InfB-resC$Expected),
                      Condition=factor(c(rep("No correction", nrow(resC)), rep("Bias-corrected", nrow(resC)))),
                      Expected=paste0("Expected: ",factor(100*rep(resC$Expected, 2)), "% HC"),
                      Group=paste0("C",c(1:nrow(resC),1:nrow(resC))),
                      Type=rep(c("Homozygous", "Heterozygous", "Homozygous", "Homozygous", "Heterozygous", "Homozygous", "Heterozygous", "Homozygous", "Heterozygous", "Homozygous", "Heterozygous"), 2*length(Cid)),
                      Replicate=rep(resC$Replicate, 2))
plotdfC$Condition <- relevel(plotdfC$Condition, "No correction")
ggplot(plotdfC, aes(x=Condition, y=Difference, color=Type, group=Group))+
  geom_point()+
  geom_line()+
  facet_wrap(~Expected)+
  theme_minimal()+
  ylab("Deviation from expected %HC (%)")


recipient <- read.table("Data/invitromixtures/Bref.txt", sep="\t")
donor <- read.table("Data/invitromixtures/Aref.txt", sep="\t")
Did <- c(1:22)[-c(6, 8, 10, 14, 18, 22)]
resD <- matrix(0, nrow = length(Did)*8, ncol = 2)
for(samp in 1:length(Did)){
  followup <- read.table(paste0("./data/invitromixtures/",samples[Did][samp]), sep="\t")
  resD[c((8*(samp-1)+1):(8*samp)),c(1:2)] <- followUp(followup, recipient, donor, bias = bias, reportMarker = TRUE, type = "inf")
}
resD <- as.data.frame(resD)
resD$Expected <- rep(100-c(0.1,0.1, 0.5, 0.5, 0.1, 0.5, 0, 10, 10, 10, 1,1,1, 5, 5, 5), each = 8)/100
resD$Replicate <- paste0("inv", rep(samples[Did], each = 8))
colnames(resD) <- c("InfM", "InfB", "Expected", "Replicate")
par(mfrow=c(1,2))
plot(resD$Expected, resD$InfM-resD$Expected, ylim=c(-0.11, 0.11))
plot(resD$Expected, resD$InfB-resD$Expected, ylim=c(-0.11, 0.11))
plotdfD <- data.frame(Difference=100*c(resD$InfM-resD$Expected, resD$InfB-resD$Expected),
                      Condition=factor(c(rep("No correction", nrow(resD)), rep("Bias-corrected", nrow(resD)))),
                      Expected=paste0("Expected: ",factor(100*rep(resD$Expected, 2)), "% HC"),
                      Group=paste0("D",c(1:nrow(resD),1:nrow(resD))),
                      Type=rep(c("Homozygous","Homozygous","Homozygous","Homozygous","Homozygous","Heterozygous", "Homozygous", "Heterozygous"), 2*length(Did)),
                      Replicate=rep(resD$Replicate, 2))
plotdfD$Condition <- relevel(plotdfD$Condition, "No correction")
ggplot(plotdfD, aes(x=Condition, y=Difference, color=Type, group=Group))+
  geom_point()+
  geom_line()+
  facet_wrap(~Expected)+
  theme_minimal()+
  ylab("Deviation from expected %HC (%)")

plotdfCD <- rbind(plotdfC, plotdfD)
plotdfCD$Expected <- factor(plotdfCD$Expected, levels=levels(factor(plotdfCD$Expected))[c(3, 1, 2, 4, 7, 5, 8, 9, 12, 10, 11, 6)])
plotdfCD$Marker <- "Type I"
ggplot(plotdfCD, aes(x=Condition, y=Difference, color=Type, group=Group))+
  geom_point()+
  geom_line()+
  facet_wrap(~Expected, nrow = 2)+
  theme_minimal()+
  xlab("")+
  ylab("Deviation from expected %HC (%)")


recipient <- read.table("Data/invitromixtures/Aref.txt", sep="\t")
donor <- read.table("Data/invitromixtures/Bref.txt", sep="\t")
resCB <- matrix(0, nrow = length(Cid)*2, ncol = 2)
for(samp in 1:length(Cid)){
  followup <- read.table(paste0("./data/invitromixtures/",samples[Cid][samp]), sep="\t")
  resCB[c((2*(samp-1)+1):(2*samp)),c(1:2)] <- followUp(followup, recipient, donor, bias = bias, reportMarker = TRUE, type = "potinf")
}
resCB <- as.data.frame(resCB)
resCB$Expected <- rep(c(0.1,0.1, 0.5, 0.5, 0.1, 0.5, 0, 10, 10, 10, 1,1,1, 5, 5, 5), each = 2)/100
resCB$Replicate <- rep(samples[Cid], each = 2)
colnames(resCB) <- c("InfM", "InfB", "Expected", "Replicate")
par(mfrow=c(1,2))
plot(resCB$Expected, resCB$InfM-resCB$Expected, ylim=c(-0.11, 0.11))
plot(resCB$Expected, resCB$InfB-resCB$Expected, ylim=c(-0.11, 0.11))
plotdfCB <- data.frame(Difference=100*c(resCB$InfM-resCB$Expected, resCB$InfB-resCB$Expected),
                       Condition=factor(c(rep("No correction", nrow(resCB)), rep("Bias-corrected", nrow(resCB)))),
                       Expected=paste0("Expected: ",factor(100*rep(resCB$Expected, 2)), "% HC"),
                       Group=c(1:nrow(resCB),1:nrow(resCB)),
                       Type="Heterozygous",
                       Replicate=rep(resCB$Replicate, 2))
plotdfCB$Condition <- relevel(plotdfCB$Condition, "No correction")
ggplot(plotdfCB, aes(x=Condition, y=Difference, color=Type, group=Group))+
  geom_point()+
  geom_line()+
  facet_wrap(~Expected, nrow = 2)+
  theme_minimal()+
  xlab("")+
  ylab("Deviation from expected %HC (%)")

recipient <- read.table("Data/invitromixtures/Bref.txt", sep="\t")
donor <- read.table("Data/invitromixtures/Aref.txt", sep="\t")
resDB <- matrix(0, nrow = length(Did)*5, ncol = 2)
for(samp in 1:length(Did)){
  followup <- read.table(paste0("./data/invitromixtures/",samples[Did][samp]), sep="\t")
  resDB[c((5*(samp-1)+1):(5*samp)),c(1:2)] <- followUp(followup, recipient, donor, bias = bias, reportMarker = TRUE, type = "potinf")
}
resDB <- as.data.frame(resDB)
resDB$Expected <- rep(100-c(0.1,0.1, 0.5, 0.5, 0.1, 0.5, 0, 10, 10, 10, 1,1,1, 5, 5, 5), each = 5)/100
resDB$Replicate <- paste0("inv", rep(samples[Did], each = 5))
colnames(resDB) <- c("InfM", "InfB", "Expected", "Replicate")
par(mfrow=c(1,2))
plot(resDB$Expected, resDB$InfM-resDB$Expected, ylim=c(-0.11, 0.11))
plot(resDB$Expected, resDB$InfB-resDB$Expected, ylim=c(-0.11, 0.11))
plotdfDB <- data.frame(Difference=100*c(resDB$InfM-resDB$Expected, resDB$InfB-resDB$Expected),
                       Condition=factor(c(rep("No correction", nrow(resDB)), rep("Bias-corrected", nrow(resDB)))),
                       Expected=paste0("Expected: ", factor(100*rep(resDB$Expected, 2)),"% HC"),
                       Group=c(1:nrow(resDB),1:nrow(resDB)),
                       Type='Heterozygous',
                       Replicate = rep(resDB$Replicate, 2))
plotdfDB$Condition <- relevel(plotdfDB$Condition, "No correction")
ggplot(plotdfDB, aes(x=Condition, y=Difference, color=Expected, group=Group))+
  geom_point()+
  geom_line()+
  theme_minimal()+
  ylab("Deviation from expected %HC (%)")

plotdfCDB <- rbind(plotdfCB, plotdfDB)
plotdfCDB$Expected <- factor(plotdfCDB$Expected, levels=levels(factor(plotdfCDB$Expected))[c(3,1,2,4,7,5,8,9,12,10,11,6)])
plotdfCDB$Marker <- "Type II"
ggplot(plotdfCDB, aes(x=Condition, y=Difference, color=Type, group=Group))+
  geom_point()+
  geom_line()+
  facet_wrap(~Expected, nrow = 2)+
  theme_minimal()+
  xlab("")+
  ylab("Deviation from expected %HC (%)")

plotall <- rbind(plotdfCD, plotdfCDB)
colnames(plotall)[c(5,7)] <- c("Constellation", "Type")
plotall$Expected<-as.character(plotall$Expected)
plotalltypeI <- plotall[plotall$Type=="Type I",]

#calculate sample averages
for(sampleid in 1:length(unique(plotall$Replicate))){
  samplelev <- unique(plotall$Replicate)[sampleid]
  nocorr <- mean(plotall$Difference[intersect(which(plotall$Replicate == samplelev),which(plotall$Condition == "No correction"))])
  corr <- mean(plotall$Difference[intersect(which(plotall$Replicate == samplelev),which(plotall$Condition == "Bias-corrected"))])
  nocorrtypeI <- mean(plotalltypeI$Difference[intersect(which(plotalltypeI$Replicate == samplelev),which(plotalltypeI$Condition == "No correction"))])
  corrtypeI <- mean(plotalltypeI$Difference[intersect(which(plotalltypeI$Replicate == samplelev),which(plotalltypeI$Condition == "Bias-corrected"))])
  plotall <- rbind(plotall, c(nocorr, "No correction", plotall$Expected[plotall$Replicate==unique(plotall$Replicate)[sampleid]][1], paste0("mean", sampleid), "Not applicable", samplelev, "Combined type I and II"))
  plotall <- rbind(plotall, c(corr, "Bias-corrected",  plotall$Expected[plotall$Replicate==unique(plotall$Replicate)[sampleid]][1], paste0("mean", sampleid), "Not applicable",  samplelev, "Combined type I and II"))
  plotall <- rbind(plotall, c(nocorrtypeI, "No correction", plotall$Expected[plotall$Replicate==unique(plotall$Replicate)[sampleid]][1], paste0("meantypeI", sampleid), "Not applicable", samplelev, "Combined type I"))
  plotall <- rbind(plotall, c(corrtypeI, "Bias-corrected",  plotall$Expected[plotall$Replicate==unique(plotall$Replicate)[sampleid]][1], paste0("meantypeI", sampleid), "Not applicable",  samplelev, "Combined type I"))
  plotall$Difference <- as.numeric(plotall$Difference)
}
plotall$Expected <- factor(plotall$Expected, levels=levels(factor(plotall$Expected))[c(3,1,2,4,7,5,8,9,12,10,11,6)])


# generate final plots
plotallmarkers <- plotall[plotall$Type %in% c("Type I", "Type II"),]
markerplot <- ggplot(plotallmarkers, aes(x=Condition, y=(Difference), color=Constellation, shape=Type, group=Group))+
  geom_point(size=1.3)+
  geom_line(size=0.4)+
  facet_wrap(~Expected, nrow = 4)+
  theme_minimal()+
  theme(legend.position = "top")+
  xlab("")+coord_flip()+
  ylab("Deviation from expected %HC (%)")+
  ylim(-9, 9)

wilcox.test(abs(plotallmarkers$Difference[plotallmarkers$Condition=="No correction"]),
            abs(plotallmarkers$Difference[plotallmarkers$Condition=="Bias-corrected"]),
            paired = TRUE)
markerbiastypei <- plotallmarkers[plotallmarkers$Type=="Type I",]
wilcox.test(abs(markerbiastypei$Difference[markerbiastypei$Condition=="No correction"]),
            abs(markerbiastypei$Difference[markerbiastypei$Condition=="Bias-corrected"]),
            paired = TRUE)
median(abs(markerbiastypei$Difference[markerbiastypei$Condition=="No correction"]))
median(abs(markerbiastypei$Difference[markerbiastypei$Condition=="Bias-corrected"]))

markerbiastypeii <- plotallmarkers[plotallmarkers$Type=="Type II",]
wilcox.test(abs(markerbiastypeii$Difference[markerbiastypeii$Condition=="No correction"]),
            abs(markerbiastypeii$Difference[markerbiastypeii$Condition=="Bias-corrected"]),
            paired = TRUE)
median(abs(markerbiastypeii$Difference[markerbiastypeii$Condition=="No correction"]))
median(abs(markerbiastypeii$Difference[markerbiastypeii$Condition=="Bias-corrected"]))


plotallsamples <- plotall[plotall$Type %in% c("Combined type I", "Combined type I and II"),]
sampleplot <- ggplot(plotallsamples, aes(x=Condition, y=(Difference), color=Type, group=Group))+
  geom_point(size=1.3)+
  geom_line(size=0.4)+
  facet_wrap(~Expected, nrow = 4)+
  theme_minimal()+
  theme(legend.position = "top")+
  xlab("")+coord_flip()+
  ylab("Deviation from expected %HC (%)")+
  ylim(-1.5, 1.5)

pdf("Figures/Figure3.pdf", height=9,width=11)
plot_grid(markerplot, sampleplot, nrow = 2, labels = c("A", "B"))
dev.off()

wilcox.test(abs(plotallsamples$Difference[plotallsamples$Condition=="No correction"]),
            abs(plotallsamples$Difference[plotallsamples$Condition=="Bias-corrected"]),
            paired = TRUE)
samplebiastypei <- plotallsamples[plotallsamples$Type=="Combined type I",]
wilcox.test(abs(samplebiastypei$Difference[samplebiastypei$Condition=="No correction"]),
            abs(samplebiastypei$Difference[samplebiastypei$Condition=="Bias-corrected"]),
            paired = TRUE)
median(abs(samplebiastypei$Difference[samplebiastypei$Condition=="No correction"]))
median(abs(samplebiastypei$Difference[samplebiastypei$Condition=="Bias-corrected"]))

samplebiastypeii <- plotallsamples[plotallsamples$Type=="Combined type I and II",]
wilcox.test(abs(samplebiastypeii$Difference[samplebiastypeii$Condition=="No correction"]),
            abs(samplebiastypeii$Difference[samplebiastypeii$Condition=="Bias-corrected"]),
            paired = TRUE)
median(abs(samplebiastypeii$Difference[samplebiastypeii$Condition=="No correction"]))
median(abs(samplebiastypeii$Difference[samplebiastypeii$Condition=="Bias-corrected"]))

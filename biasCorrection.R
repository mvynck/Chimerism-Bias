###########################################
### libraries                           ###
###########################################
library("ggplot2")
library("RColorBrewer")
library("dplyr")
library("cowplot")

###########################################
### source functinos                    ###
###########################################

setwd("~/SetAppropriateDir")
source("Code/biasFuns.R")
source("Code/dataBiasReliability.R")


###########################################
### load sample IDs                     ###
###########################################

# read
samples <- read.csv("Data/RecipientDonorFollowupSampleID.csv", sep=",")
colnames(samples) <- c("Recipient", "Donor", "Followup", "Run", "RunRecipient", "RunDonor")

# sample summary statistics
length(unique(samples$Recipient))
length(unique(samples$Donor))
length(samples$Followup)


###########################################
### load VAF data                       ###
###########################################

# initialize data frames to store VAFs
dat <- data.frame(matrix(ncol = 49, nrow = 0))
colnames(dat) <- c("file", paste0("MarkerRef", 1:24), paste0("MarkerAlt", 1:24))
VAFdat <- data.frame(matrix(ncol = 25, nrow = 0))
colnames(VAFdat) <- c("file", paste0("VAFRef", 1:24))
res <- matrix(0, nrow=nrow(samples), ncol=24)
read.depth.pre <- vector()
read.depth.post <- vector()
  
# load data into data frames
files <- list.files(path = "Data/Refdata")
for(sample in 1:nrow(samples)){
  recipient <- read.table(paste0("Data/Refdata/", samples$Recipient[sample], "RO.txt"), sep="\t")
  donor <- read.table(paste0("Data/Refdata/", samples$Donor[sample], "RD.txt"), sep="\t")
  
  followup <- read.table(paste0("Data/Monitoringdata/", samples$Followup[sample], ".txt"), sep="\t")  
  res[sample,] <- followUp(followup, recipient, donor, bias)
  
  # add to vector for read depth overview
  read.depth.pre <- c(read.depth.pre, recipient$V2, donor$V2)
  read.depth.post <- c(read.depth.post, followup$V2)
  
  #for bias estimation in this subset
  filename <- strsplit(samples$Recipient[sample], "\\.")[[1]][1]
  dat[sample*3-2,] <- c(filename, recipient[,3], recipient[,4])
  VAFdat[sample*3-2,] <- c(filename, recipient[,3]/(recipient[,3]+recipient[,4]))
  filename <- strsplit(samples$Donor[sample], "\\.")[[1]][1]
  dat[sample*3-1,] <- c(filename, donor[,3], donor[,4])
  VAFdat[sample*3-1,] <- c(filename, donor[,3]/(donor[,3]+donor[,4]))
  filename <- strsplit(samples$Followup[sample], "\\.")[[1]][1]
  dat[sample*3,] <- c(filename, followup[,3], followup[,4])
  VAFdat[sample*3,] <- c(filename, followup[,3]/(followup[,3]+followup[,4]))
}

# some reading depth summaries
ggplot(data.frame(read.depth.post), aes(x=read.depth.post))+geom_histogram(binwidth = 1000)
ggplot(data.frame(read.depth.pre), aes(x=read.depth.pre))+geom_histogram(binwidth = 25)


# screening data will be loaded multiple times when 
#    multiple follow up samples are available
# remove duplicates
VAFdat <- VAFdat[!duplicated(VAFdat$file),]
VAFdat.matrix <- matrix(as.numeric(as.matrix(VAFdat[,2:25])), ncol = 24, nrow = nrow(VAFdat))


###########################################
### estimate bias                       ###
###########################################

# note that in the manuscript bias is defined as the difference between the
#  VAF and 0.5, here we calculate the %ref - 0.5, thus the bias is obtained as
#  -(%ref - 0.5) = -(1-VAF-0.5) = (VAF - 0.5)
deviation <- c()
reliability <- c()
for(marker in 1:24){
  deviation[marker] <- -median(VAFdat.matrix[as.logical((VAFdat.matrix[, marker] < 0.60) * (VAFdat.matrix[, marker] > 0.40)), marker]-0.5)
}
median(abs(deviation))
sort(deviation)
bias <- deviation

###########################################
### reload data                         ###
### obtain bias-corrected VAFs          ###
### store stochastic effects for later  ###
###########################################

stochasticEffects <- data.frame(type=NULL, HC=NULL, Reads=NULL, Sample= NULL)
for(sample in 1:nrow(samples)){
  recipient <- read.table(paste0("Data/Refdata/", samples$Recipient[sample], "RO.txt"), sep="\t")
  donor <- read.table(paste0("Data/Refdata/", samples$Donor[sample], "RD.txt"), sep="\t")
  followup <- read.table(paste0("Data/Monitoringdata/", samples$Followup[sample], ".txt"), sep="\t")  
  res[sample,] <- followUp(followup, recipient, donor, bias)
  stochasticEffects <- rbind(stochasticEffects, stochasticityReadDepth(followup, recipient, donor, bias, samples$Followup[sample]))
  
  #for bias estimation in this subset
  filename <- strsplit(samples$Recipient[sample], "\\.")[[1]][1]
  dat[sample*3-2,] <- c(filename, recipient[,3], recipient[,4])
  VAFdat[sample*3-2,] <- c(filename, recipient[,3]/(recipient[,3]+recipient[,4]))
  filename <- strsplit(samples$Donor[sample], "\\.")[[1]][1]
  dat[sample*3-1,] <- c(filename, donor[,3], donor[,4])
  VAFdat[sample*3-1,] <- c(filename, donor[,3]/(donor[,3]+donor[,4]))
  filename <- strsplit(samples$Followup[sample], "\\.")[[1]][1]
  dat[sample*3,] <- c(filename, followup[,3], followup[,4])
  VAFdat[sample*3,] <- c(filename, followup[,3]/(followup[,3]+followup[,4]))
}


###########################################
### compare different measurement types ###
###########################################

res.sub <- res



# all samples
IDs=1:nrow(res.sub)

# get IDs for samples with at least 3 type I/II markers
IDb=which(res.sub[,4]>=3)
IDc=which(res.sub[,12]>=3)
IDss=intersect(IDs,IDc)
IDs=intersect(IDs,IDb)

# coefficient of variation
a=res.sub[IDs,3]/res.sub[IDs,1]
b=res.sub[IDs,7]/res.sub[IDs,5]


# final multipanel plot
df1 <- data.frame(typeIuncor=res.sub[IDs,5], y=b/a)
df2 <- data.frame(typeIuncor=res.sub[,5], y=res.sub[,5]-res.sub[,1])

df3 <- data.frame(typeIuncor=res.sub[,5], y=res.sub[,9]-res.sub[,5])
df4 <- data.frame(typeIuncor=res.sub[,5], y=res.sub[,13]-res.sub[,5])

df5 <- data.frame(typeIuncor=res.sub[,5], y=res.sub[,17]-res.sub[,5])
df6 <- data.frame(typeIuncor=res.sub[,5], y=res.sub[,21]-res.sub[,5])

p1 <- ggplot(df1, aes(typeIuncor, y))+
  geom_point()+
  xlab("%HC, type-I (corrected)")+
  ylab("Coef. of variation type-I (corr.) / (uncorr.)")+
  theme_minimal()+
  scale_x_continuous(trans='log10')
p2 <- ggplot(df2, aes(typeIuncor, y))+
  geom_point()+
  xlab("%HC, type-I (corrected)")+
  ylab("%HC, type-I (uncorr.) - type-I (corr.)")+
  theme_minimal()+
  scale_x_continuous(trans='log10')+
  ylim(-0.03, 0.03)
p3 <- ggplot(df3, aes(typeIuncor, y))+
  geom_point()+
  xlab("%HC, type-I (corrected)")+
  ylab("%HC, type-II (uncorr.) - type-I (corr.)")+
  theme_minimal()+
  scale_x_continuous(trans='log10')+
  ylim(-0.09, 0.09)
p4 <- ggplot(df4, aes(typeIuncor, y))+
  geom_point()+
  xlab("%HC, type-I (corrected)")+
  ylab("%HC, type-II (corr.) - type-I (corr.)")+
  theme_minimal()+
  scale_x_continuous(trans='log10')+
  ylim(-0.09, 0.09)
p5 <- ggplot(df5, aes(typeIuncor, y))+
  geom_point()+
  xlab("%HC, type-I (corrected)")+
  ylab("%HC, combined type-I/-II (uncorr.) - type-I (corr.)")+
  theme_minimal()+
  scale_x_continuous(trans='log10')+
  ylim(-0.03, 0.03)
p6 <- ggplot(df6, aes(typeIuncor, y))+
  geom_point()+
  xlab("%HC, type-I (corrected)")+
  ylab("%HC, combined type-I/-II (corr.) - type-I (corr.)")+
  theme_minimal()+
  scale_x_continuous(trans='log10')+
  ylim(-0.03, 0.03)


pdf("~/Desktop/Figure3.pdf", height=12.7,width=9)
plot_grid(p1, p2, p3, p4, p5, p6, nrow = 3, ncol = 2, labels="AUTO")
dev.off()
###########################################
### stochastic effects                  ###
###########################################

#subtract sample-wise average InfB from PotInfB
#retain only samples with low observed chimerism (no expected bias in informative markers)
stocTemp <- stochasticEffects
retainSamples <- stochasticEffects[stochasticEffects$type == "InfB",]
retainSamples <- unique(retainSamples$sample[retainSamples$HC < 0.01])
stochasticEffects <- stochasticEffects[stochasticEffects$sample %in% retainSamples,]
stochasticEffects$PotInfBaseline <- stochasticEffects$HC

for(sample in 1:length(unique(stochasticEffects$sample))){
  sampleLocs <- which(stochasticEffects$sample == unique(stochasticEffects$sample)[sample])
  potInfLocs <- which(stochasticEffects$type == "PotInfB")
  subset <- stochasticEffects[stochasticEffects$sample == unique(stochasticEffects$sample)[sample],]
  subsetInf <- subset[subset$type == "InfB",]
  stochasticEffects$PotInfBaseline[intersect(sampleLocs,potInfLocs)] <- stochasticEffects$HC[intersect(sampleLocs,potInfLocs)] - mean(subsetInf$HC)
}


par(mfrow=c(1,1))
summary(lm(abs(stochasticEffects$PotInfBaseline[stochasticEffects$type=="PotInfB"]
       )~stochasticEffects$Reads[stochasticEffects$type=="PotInfB"]))
plot(stochasticEffects$Reads[stochasticEffects$type=="PotInfB"], 100*abs(stochasticEffects$PotInfBaseline[stochasticEffects$type=="PotInfB"]),
     xlab = "Number of reads", ylab = "Observed absolute bias (%)")

# add expected 95% deviation, multiply by 4 as VAF is multipled by 2 to obtain HC,
#   Var(aX) = a**2*Var(X)
lines(seq(10000,90000,10),100*qnorm(0.975)*sqrt(4*0.5*0.5/seq(10000,90000,10)), type="l",col="orange")

ggptdf <- data.frame(x = stochasticEffects$Reads[stochasticEffects$type=="PotInfB"], y=100*abs(stochasticEffects$PotInfBaseline[stochasticEffects$type=="PotInfB"]))
gglinedf <- data.frame(x = seq(10000,90000,10), y = 100*qnorm(0.975)*sqrt(4*0.5*0.5/seq(10000,90000,10)))
pt1 <- ggplot(ggptdf, aes(x, y, color="black"))+geom_point()+
          theme_minimal()+
          xlab("Reading depth")+
          ylab("Absolute deviation (%)")+
          geom_line(data = gglinedf, aes(color = "blue"))+
          scale_colour_manual(name = c("point", "line"),
                              values = c("black", "blue"),
                              label = c("Observed", "Theoretical 95% quantile"),
                              guide = guide_legend(override.aes = list(
                                linetype = c("blank", "solid"),
                                shape = c(16, NA))))+
          theme(legend.position = "top", legend.title = element_blank())

# check the above approach by monte carlo sampling
yvals = vector()
xvals <- seq(10000,90000,5000)
for(i in 1:length(xvals)){
   n <- xvals[i]
   # 
   #   replicate rbinom: generate 10000 samples from a binomial with prob 0.5 (heterozygous marker)
   #   colSums / n: obtain 10000 averages
   #   -0.5: estimate bias
   #   abs: take absolute value as biases will be symmetric (more reliable tail quantile)
   #   quantile: obtain the 95% quantile (= two-tailed 2.5%, i.e. absolute value was taken)
   yvals[i] <- quantile(abs(colSums(replicate(10000, rbinom(n, size = 1, prob = 0.5)))/n-0.5)*2, 0.95)
}
lines(xvals, yvals*100, col = "blue") # OK, identical

# calculate % of observations below 95% threshold
stocSubset <- stochasticEffects[stochasticEffects$type=="PotInfB",]
below95 <- vector()
for(i in 1:nrow(stocSubset)){
  below95[i] <- qnorm(0.975)*sqrt(4*0.5*0.5/stocSubset$Reads[i])
}
mean(abs(stocSubset$PotInfBaseline) < below95) #65.2%
binom.test(sum(abs(stocSubset$PotInfBaseline) < below95), length(below95), 0.95) # p < 0.001
# 50% = typical bias due to binomial variability, function of sequencing depth
data.frame(Depth=seq(10000, 100000, 10000),Expected=qnorm(0.75)*sqrt(4*0.5*0.5/seq(10000, 100000, 10000)))

# typically, at least three markers required, stochastic (+remainder) effects may partially cancel out
# simulate by sampling in 5000 read bins
# simulate with increasing number of markers
stochasticEffectsSampling <- stocSubset
stochasticEffectsSampling$Reads <- 5000*round(stochasticEffectsSampling$Reads/5000, 0)
absSamplingNSims <- 10000
absBiasSampling <- matrix(0, nrow = absSamplingNSims, ncol = 11)
for(reads in seq(15000, 65000, 5000)){
  subset <- stochasticEffectsSampling$PotInfBaseline[stochasticEffectsSampling$Reads == reads]
  for(sim in 1:absSamplingNSims){
    absBiasSampling[sim, reads/5000 - 2] <- abs(mean(sample(subset, 3, replace = FALSE)))
  }
}
boxplot(absBiasSampling*100, names = seq(15000, 65000, 5000), xlab = "Reads", ylab = "Absolute deviation (%HC)")
abline(h=0.5, col = "blue")

stochasticEffectsSampling <- stocSubset
stochasticEffectsSampling$Reads <- 5000*round(stochasticEffectsSampling$Reads/5000, 0)
absSamplingNSims <- 10000
absBiasSampling2 <- matrix(0, nrow = absSamplingNSims, ncol = 10)
reads <- median(stochasticEffectsSampling$Reads) # 40000
for(markers in 1:10){
  subset <- stochasticEffectsSampling$PotInfBaseline[stochasticEffectsSampling$Reads == reads]
  for(sim in 1:absSamplingNSims){
    absBiasSampling2[sim, markers] <- abs(mean(sample(subset, markers, replace = FALSE)))
  }
}
boxplot(absBiasSampling2*100, names = 1:10, xlab = "Markers", ylab = "Absolute deviation (%HC)")

# theoretical
absBiasSampling3 <- matrix(0, nrow = absSamplingNSims, ncol = 10)
for(markers in 1:10){
  for(sim in 1:absSamplingNSims){
    genDevs <- rbinom(markers, reads, 0.5)/reads - 0.5
    absBiasSampling3[sim, markers] <- abs(mean(genDevs))
  }
}
boxplot(absBiasSampling3*100, names = 1:10, xlab = "Markers", ylab = "Absolute deviation (%HC)")

#combined figure
fig.df <- data.frame(obsDevs = 100*c(absBiasSampling2, absBiasSampling3), markers = factor(rep(rep(1:10, each = 10000), 2)), type = rep(c("Observed", "Theoretical"), each = 100000))
pt2 <- ggplot(fig.df, aes(x=markers, y=obsDevs, fill = type))+
        geom_boxplot()+
        theme_minimal()+
        xlab("Number of markers")+
        ylab("Absolute deviation (%)")+
        theme(legend.position="top", legend.title = element_blank())
  

plot_grid(pt1, pt2, ncol = 2, labels = c("A", "B"))

###########################################
### run specific biases                 ###
###########################################


# get non informative bias data from follow up samples
run.bias <- data.frame(matrix(0, ncol = 4, nrow = 0))
colnames(run.bias) <- c("marker", "VAF", "run", "sample")
for(sample in 1:nrow(samples)){
  recipient <- read.table(paste0("Data/Refdata/", samples$Recipient[sample], "RO.txt"), sep="\t")
  donor <- read.table(paste0("Data/Refdata/", samples$Donor[sample], "RD.txt"), sep="\t")
  followup <- read.table(paste0("Data/Monitoringdata/", samples$Followup[sample], ".txt"), sep="\t")  
  VAF <- followup$V3/(followup$V3+followup$V4)
  VAF.NI <- VAF[which(!1:24%in%sort(c(informative(recipient, donor), potentiallyInformative(recipient, donor))))]
  MarkerID.NI <- which(!1:24%in%sort(c(informative(recipient, donor), potentiallyInformative(recipient, donor))))
  run.bias <- rbind(run.bias, data.frame(marker=MarkerID.NI, VAF = VAF.NI, run = rep(samples$Run[sample], length(VAF.NI)), Sample = samples$Followup[sample]))
}
# remove homozygous markers
run.bias <- run.bias[!as.logical((run.bias$VAF < 0.4) + (run.bias$VAF > 0.6)),]
run.bias.fu <- run.bias


# next get non informative bias data from reference samples

ref.data <- read.csv("Data/Run-Ref.csv")
colnames(ref.data) <- c("run", "filename")
for(sample in 1:nrow(ref.data)){
  ref <- read.table(paste0("Data/Refdata/", ref.data$filename[sample]), sep="\t")
  VAF <- ref$V3/(ref$V3+ref$V4)
  VAF.het <- VAF[!(as.logical(c((VAF > 0.6) + (VAF < 0.4))))]
  MarkerID.het <- which(!as.logical(c((VAF > 0.6) + (VAF < 0.4))))
  run.bias <- rbind(run.bias, data.frame(marker=MarkerID.het, VAF = VAF.het, run = rep(ref.data$run[sample], length(VAF.het)), Sample = samples$Followup[sample]))
}
run.bias.pre <- run.bias[-c(1:nrow(run.bias.fu)),]


# reliability color
reliability.run.marker <- table(paste0(as.numeric(run.bias$run), ".",as.numeric(run.bias$marker)))
order <- paste0(rep(c(10, 11, 12, 13, 14, 15, 7, 8, 9), 24), ".", rep(1:24, each = 9))
hits <- rep(0, length(order))
for(element in 1:length(order)){
  hits[element] <- reliability.run.marker[which(names(reliability.run.marker) == order[element])]
  if(hits[element]>7){hits[element] <- 7}
}

plot(1:nrow(run.bias), run.bias$VAF, ylim = c(0.4, 0.6))
name.vec <- c(rep("", 4), "1", rep("", 4),
              rep("", 4), "2", rep("", 4),
              rep("", 4), "3", rep("", 4),
              rep("", 4), "4", rep("", 4),
              rep("", 4), "5", rep("", 4),
              rep("", 4), "6", rep("", 4),
              rep("", 4), "7", rep("", 4),
              rep("", 4), "8", rep("", 4),
              rep("", 4), "9", rep("", 4),
              rep("", 4), "10", rep("", 4),
              rep("", 4), "11", rep("", 4),
              rep("", 4), "12", rep("", 4),
              rep("", 4), "13", rep("", 4),
              rep("", 4), "14", rep("", 4),
              rep("", 4), "15", rep("", 4),
              rep("", 4), "16", rep("", 4),
              rep("", 4), "17", rep("", 4),
              rep("", 4), "18", rep("", 4),
              rep("", 4), "19", rep("", 4),
              rep("", 4), "20", rep("", 4),
              rep("", 4), "21", rep("", 4),
              rep("", 4), "22", rep("", 4),
              rep("", 4), "23", rep("", 4),
              rep("", 4), "24", rep("", 4))

table(table(paste0(run.bias$marker, "-",run.bias$run)))

run.bias$marker <- as.factor(run.bias$marker)

fit.lm1 <- lm(VAF~marker+run, data=run.bias)
anova(fit.lm1)


run.bias.for.plot <- run.bias
run.bias.for.plot$run <-
  recode(run.bias.for.plot$run, "7" = "a",
       "8" = "b",
       "9" = "c",
       "10" = "d",
       "11" = "e",
       "12" = "f",
       "13" = "g",
       "14" = "h",
       "15" = "i")
colnames(run.bias.for.plot) <- c("Marker", "VAF", "Run", "Sample")

pt1 <- ggplot(run.bias.for.plot, aes(x=Marker,y=100*VAF, color=Run)) +
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_dodge(width = 0.75), aes(group = Run), size = 1)+
  theme_bw()+
  theme(panel.border = element_blank(), legend.position = "none")+
  xlab("Marker ID")+
  ylab("Variant allele frequency (%)")+
  scale_colour_brewer(palette = "Set1")
run.specific.bias <- data.frame(Run = rep(7:15, each=24), Marker = rep(1:24, 9), Bias = 0)
count <- 1
for(run in 7:15){
  for(marker in 1:24){
    run.bias.temp <- run.bias[run.bias$marker==marker,]
    run.bias.temp <- run.bias.temp[run.bias.temp$run==run,]
    run.specific.bias$Bias[count] <- median(run.bias.temp$VAF) - 0.5
    count <- count + 1
  }
}
plot(1:216, run.specific.bias$Bias[order(run.specific.bias$Marker)])


# compare estimated bias in pre-transplant and post-transplant samples
#   reliability related to sequencing depth which is higher in post-
#   transplant samples?
bias.vec <- c()
n.vec <- c()
for(i in 1:24){
  bias.vec[i] <- mean(run.bias$VAF[run.bias$marker==i])-0.5
  n.vec[i] <- length(run.bias$VAF[run.bias$marker==i])
}
bias.vec.fu <- c()
n.vec.fu <- c()
for(i in 1:24){
  bias.vec.fu[i] <- mean(run.bias.fu$VAF[run.bias.fu$marker==i])-0.5
  n.vec.fu[i] <- length(run.bias.fu$VAF[run.bias.fu$marker==i])
}
bias.vec.pre <- c()
n.vec.pre <- c()
for(i in 1:24){
  bias.vec.pre[i] <- mean(run.bias.pre$VAF[run.bias.pre$marker==i])-0.5
  n.vec.pre[i] <- length(run.bias.pre$VAF[run.bias.pre$marker==i])
}
cor(bias.vec.fu, bias.vec.pre)
lm(bias.vec.fu ~ bias.vec.pre)
plot(bias.vec.fu, bias.vec.pre, xlab= "Bias in follow-up samples", ylab="Bias in pre-transplant samples")
abline(a = 0, b = 1, lty = 3)
text(0.02, 0.04, "r = 0.995")
legend("topleft", "Identity line", lty = 3)

# plot bias sorted from low to high with marker IDs

df.bias <- data.frame(bias=bias.vec*100,ID=c(1:24))
df.bias$ID <- factor(df.bias$ID, levels=c(1:24)[order(bias.vec)])
pt2 <- ggplot(df.bias, aes(x=ID, y=bias))+
  geom_bar(stat="identity")+
  theme_minimal()+
  scale_y_continuous("Bias (%)", breaks = -5.0:5.0)+xlab("Marker ID")+
  coord_flip()

# binomial variation does not explain
expected <- 0.5
sd <- sqrt(expected*(1-expected))
upper <- expected + qnorm(0.975)*sd/sqrt(c(5000:30000))
lower <- expected - qnorm(0.975)*sd/sqrt(c(5000:30000))
plot(5000:30000, upper, type="l", ylim=c(0.45,0.55), xlab = "number of reads")
lines(5000:30000, lower)

plot_grid(pt2, pt1, ncol = 2, rel_widths = c(2, 5), labels = c("A", "B"))


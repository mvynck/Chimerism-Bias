adjustMarker <- function(VAF, type, bias){
  # perform bias correction
  #  0  = ++ in the donor, -- in the recipient
  #  1  = -- in the donor, ++ in the recipient
  #  10 = ++ in the donor, +- in the recipient
  #  11 = -- in the donor, +- in the recipient
  #  20 = +- in the donor, ++ in the recipient
  #  21 = +- in the donor, -- in the recipient
  if(type==21){
    adjusted <- VAF/(0.50 + bias)
    return(1 - adjusted)
  }
  if(type==20){
    adjusted <- (VAF-1)/(-0.50+bias)
    return(1 - adjusted)
  }
  if(type==11){
    adjusted <- 2*VAF/(1 + 2*bias)
    return(adjusted)
  }
  if(type==10){
    adjusted <- 2*(VAF - 1)/(-1 + 2*bias)
    return(adjusted)
  }
  if(type==1){
    adjusted <- VAF/(1)
    return(adjusted)
  }
  if(type==0){
    adjusted <- (VAF - 1)/(-1)
    return(adjusted)
  }
}




informative <- function(recipient, donor){
  # calculate recipient and donor VAFs
  RecVAF <- recipient$V3/(recipient$V3+recipient$V4)
  DonVAF <- donor$V3/(donor$V3+donor$V4)
  # informative markers: homozygous in donor and not
  #   the same alleles in recipient
  #   homozygous is considered any VAF less than 10%
  #   or more than 90%
  infID <- c(setdiff(which(DonVAF < 0.10), which(RecVAF < 0.10)), 
             setdiff(which(DonVAF > 0.90), which(RecVAF > 0.90)))
  return(sort(infID))
}

informativeType <- function(recipient, donor){
  # get the marker type for informative markers
  # first get informative marker IDs    
  infMark <- informative(recipient, donor)
  
  # calculate the VAFs in donor and recipient
  # screening samples
  DonVAF <- donor$V3/(donor$V3+donor$V4)
  RecVAF <- recipient$V3/(recipient$V3+recipient$V4)
  
  # define a marker type
  #  0  = ++ in the donor, -- in the recipient
  #  1  = -- in the donor, ++ in the recipient
  #  10 = ++ in the donor, +- in the recipient
  #  11 = -- in the donor, +- in the recipient
  type <- c(DonVAF[infMark] < 0.10) + 10*c(RecVAF[infMark] > 0.10)*c(RecVAF[infMark] < 0.90)
  return(type)
}

potentiallyInformative <- function(recipient, donor){
  # calculate recipient and donor VAFs
  RecVAF <- recipient$V3/(recipient$V3+recipient$V4)
  DonVAF <- donor$V3/(donor$V3+donor$V4)
  # potentially informative markers: homozygous in recipient
  #   and heterozygous in donor
  potInfID <- c(intersect(intersect(which(DonVAF > 0.10), which(DonVAF < 0.90)), which(RecVAF < 0.10)), 
                intersect(intersect(which(DonVAF > 0.10), which(DonVAF < 0.90)), which(RecVAF > 0.90)))
  return(sort(potInfID))
}

potentialylInformativeType <- function(recipient, donor){
  potInfMark <- potentiallyInformative(recipient, donor)
  RecVAF <- recipient$V3/(recipient$V3+recipient$V4)
  # define a marker type
  #  20 = +- in the donor, ++ in the recipient
  #  21 = +- in the donor, -- in the recipient
  type <- c(RecVAF[potInfMark] < 0.10) + 20
  return(type)
}

hostPercentage <- function(VAF, type){
  # calculate host percentage
  # for different marker tyoes
  # without bias correction
  altHomoPerc <- VAF[type==1]
  refHomoPerc <- 1-VAF[type==0]
  altHeteroPerc <- VAF[type==11]*2
  refHeteroPerc <- (1-VAF[type==10])*2
  potInfAlt <- (0.5-VAF[type==21])*2
  potInfRef <- (VAF[type==20]-0.5)*2
  return(c(altHomoPerc,refHomoPerc, altHeteroPerc, refHeteroPerc, potInfAlt, potInfRef))
}



correctBias <- function(type = 0, HC = 0.50, bias = 0.10){
  #  0  = ++ in the donor, -- in the recipient
  #  1  = -- in the donor, ++ in the recipient
  #  10 = ++ in the donor, +- in the recipient
  #  11 = -- in the donor, +- in the recipient
  #  20 = +- in the donor, ++ in the recipient
  #  21 = +- in the donor, -- in the recipient
  
  b = bias
  #AA donor, aa recipient
  if(type == 0){
    corrected <- (HC-2*HC*b)/(-4*b*HC+2*b+1)
  }
  #aa donor, Aa recipient
  if(type == 11){
    corrected <- (2*b*HC+HC)/(2*b*(HC-1)+1)
  }
  #AA donor, Aa recipient
  if(type == 10){
    corrected <- (HC - 2*b*HC)/(-2*b*HC+2*b+1)
  }
  #aa donor, AA recipient
  if(type == 1){
    corrected <- (2*b*HC+HC)/(4*b*HC-2*b+1)
  }  
  #Aa donor, AA recipient
  if(type == 20){
    corrected <- (2*b + HC)/(2*b*HC +1)
  }  
  #Aa donor, aa recipient
  if(type == 21){
    corrected <- (2*b-HC)/(2*b*HC-1)
  }  
  return(corrected)
}

followUp <- function(followup, recipient, donor, bias){
  # calculate marker VAFs in the FU sample
  VAF <- followup$V3/(followup$V3+followup$V4)
  
  # get (pot.) inf. markers and their type
  infMark <- informative(recipient, donor)
  infType <- informativeType(recipient, donor)
  potInfMark <- potentiallyInformative(recipient, donor)
  potInfType <- potentialylInformativeType(recipient, donor)
  
  # get host percentages, with and without
  # bias correction
  # for informative markers
  if(length(infMark) > 0){
    infHC <- mapply(hostPercentage, VAF[infMark], infType)
    infHCBias <- c()
    for(marker in 1:length(infMark)){
      infHCBias[marker] <- correctBias(infType[marker], infHC[marker], bias = bias[infMark[marker]])
    }
    # if no inf. markers, return NA
  } else {
    infHC <- NA
    infHCBias <- NA
  }
  
  # and potentially informative markers
  if(length(potInfMark) > 0){
    potInfHC <- mapply(hostPercentage, VAF[potInfMark], potInfType)
    potInfHCBias <- c()
    for(marker in 1:length(potInfMark)){
      potInfHCBias[marker] <- correctBias(potInfType[marker], potInfHC[marker], bias = bias[potInfMark[marker]])
    }
    # if no pot. inf. markers, return NA
  } else {
    potInfHC <- NA
    potInfHCBias <- NA
  }
  
  # return 1x24 vector with summary statistics of host percentages
  # for inf. markers and pot. inf. markers, with and without bias
  # correction separately
  return(c(mean(infHC), median(infHC), sd(infHC), length(infHC),
           mean(infHCBias), median(infHCBias), sd(infHCBias), length(infHCBias),
           mean(potInfHC), median(potInfHC), sd(potInfHC), length(potInfHC),
           mean(potInfHCBias), median(potInfHCBias), sd(potInfHCBias), length(potInfHCBias),
           mean(c(infHC,potInfHC)), median(c(infHC, potInfHC)), sd(c(infHC, potInfHC)), length(c(infHC, potInfHC)),
           mean(c(infHCBias,potInfHCBias)), median(c(infHCBias, potInfHCBias)), sd(c(infHCBias, potInfHCBias)), length(c(infHCBias, potInfHCBias))))
}

stochasticityReadDepth <- function(followup, recipient, donor, bias, sample){
  # calculate marker VAFs in the FU sample
  VAF <- followup$V3/(followup$V3+followup$V4)
  
  # get (pot.) inf. markers and their type
  infMark <- informative(recipient, donor)
  infType <- informativeType(recipient, donor)
  potInfMark <- potentiallyInformative(recipient, donor)
  potInfType <- potentialylInformativeType(recipient, donor)
  
  # get host percentages, with and without
  # bias correction
  # for informative markers
  if(length(infMark) > 0){
    infHC <- mapply(hostPercentage, VAF[infMark], infType)
    infHCBias <- c()
    for(marker in 1:length(infMark)){
      infHCBias[marker] <- correctBias(infType[marker], infHC[marker], bias = bias[infMark[marker]])
    }
  # if no inf. markers, return NA
  } else {
    infHC <- NA
    infHCBias <- NA
  }
  
  # and potentially informative markers
  if(length(potInfMark) > 0){
    potInfHC <- mapply(hostPercentage, VAF[potInfMark], potInfType)
    potInfHCBias <- c()
    for(marker in 1:length(potInfMark)){
      potInfHCBias[marker] <- correctBias(potInfType[marker], potInfHC[marker], bias = bias[potInfMark[marker]])
    }
  # if no pot. inf. markers, return NA
  } else {
    potInfHC <- NA
    potInfHCBias <- NA
  }
  data.frame(sample=c(rep(sample, 2*length(infHC)+2*length(potInfHC))),
             type=c(rep("Inf", length(infHC)),
                    rep("InfB", length(infHCBias)), 
                    rep("potInf", length(potInfHC)), 
                    rep("PotInfB", length(potInfHCBias))),
             HC=c(infHC, infHCBias, potInfHC, potInfHCBias),
             Reads=c(followup$V3+followup$V4)[c(infMark, infMark, potInfMark, potInfMark)])

}

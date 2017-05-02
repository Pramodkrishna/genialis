#Step1
getwd()
setwd("/home/pramod/Documents/Kaggle/genialis/")
leu.data <- read.csv("leukemia (copy).csv",header = T, sep = ",",stringsAsFactors = F) #loading gene expression file.
d.all <- leu.data[,c(2:28,40:59)] 
d.all[d.all == 0] <- 1 #converting zero values to 1
d.log_all <- log(d.all,10)  # log(1) would return "0"; Normalization of data with log.  
d.log_all[is.na(d.log_all)] <- 0 
d.log_all <- round(d.log_all,2)
all.log <- cbind(leu.data$gene.patient,d.log_all) 
RM <- as.data.frame(rowMeans(all.log[,2:48])) #Rowmean 
RSD <- as.data.frame(transform(d.log_all,sd = apply(d.log_all,1,sd, na.rm = T))) #Sum of standard deviation
All.combined <- cbind(leu.data$gene.patient,RM$`rowMeans(all.log[, 2:48])`,RSD$sd) #Combining into One data frame 
All.data <- as.data.frame(All.combined,stringsAsFactors = F)
colnames(All.data) <- c("gene","mean.all","sd.all")
All.data$sd.all <- as.numeric(All.data$sd.all)
All.data$mean.all <- as.numeric(All.data$mean.all)
All.data$mean.all <- round(All.data$mean.all,2)
All.data$sd.all <- round(All.data$sd.all,2) #Final data for ALL class
head(All.data)
#AML data 
d.aml <- leu.data[,c(29:39,60:73)]
d.aml[d.aml == 0] <- 1        #converting zero values to 1
d.log_aml <- log(d.aml,10)    # log(1) would return "0"; Normalization of data with log.
d.log_aml[is.na(d.log_aml)] <- 0
d.log_aml <- round(d.log_aml,2)
aml.log <- cbind(leu.data$gene.patient,d.log_aml)
RM.aml <- as.data.frame(rowMeans(aml.log[,2:25])) #Row mean calculation.
RSD.aml <- as.data.frame(transform(aml.log,sd = apply(aml.log,1,sd, na.rm = T))) #Standard deviation calculation 
AML.combined <- cbind(leu.data$gene.patient,RM.aml$`rowMeans(aml.log[, 2:25])`,RSD.aml$sd)
AML.data <- as.data.frame(AML.combined,stringsAsFactors = F)
colnames(AML.data) <- c("gene","mean.aml","sd.aml")
AML.data$sd.aml <- as.numeric(AML.data$sd.aml)
AML.data$mean.aml <- as.numeric(AML.data$mean.aml)
AML.data$mean.aml <- round(AML.data$mean.aml,2)
AML.data$sd.aml <- round(AML.data$sd.aml,2)
head(AML.data)
#Merge data sets.
combined.data <- merge.data.frame( AML.data, All.data,by = "gene")
head(combined.data)
View(combined.data)
#############################Ranking based on Snr###############################
#Signal to Noise ratio. 
# Difference between means of each class divided by standard deviation of the classes for each gene.
for(i in 1:length(combined.data)){
  for(k in 1 : length(combined.data))
  combined.data$diff <- combined.data$mean.all -  combined.data$mean.aml
  combined.data$stdev <- combined.data$sd.aml + combined.data$sd.all
  combined.data$sn <- (combined.data$diff)/(combined.data$stdev)
  print(head(combined.data))
  }
 
#Ranking them according to their Siganl to Noise ratio.
combined.data$sn[combined.data$sn == "NaN"] <- 0
combined.data$sn <- round(combined.data$sn,2)
combined.data <- combined.data[order(combined.data$sn,decreasing = T),]
combined.data$rank <- rank(combined.data$sn,by = combined.data$sn)


View(combined.data)
#######################################
library(qdap)
pathway <- readLines("pathways.txt")                                       #loading pathway.txt file
pathlist <- strsplit(pathway,split = "\t")                                 
path.data <- do.call(rbind.data.frame,pathlist)
path.data <- data.frame(lapply(path.data, as.character), stringsAsFactors=FALSE)
colnames(path.data)[1] <- "pathway"
colnames(path.data)[2] <- "description"
x.df <- combined.data$gene


#Replacing the gene with values.
if(i in path.data$c..MAPK14....ace2Pathway....acetaminophenPathway....CHRNG....ARPC3... %n% x.df)
{
  path.data$c..MAPK14....ace2Pathway....acetaminophenPathway....CHRNG....ARPC3... <-  replace(path.data$c..MAPK14....ace2Pathway....acetaminophenPathway....CHRNG....ARPC3...,values = combined.data$sn) 
  print(path.data$c..MAPK14....ace2Pathway....acetaminophenPathway....CHRNG....ARPC3...)
}


if(i in path.data$c..TRAF2....COL4A1....PTGS2....TERT....ABI.2....GNGT1....GNAS... %n%  x.df)
{
  path.data$c..TRAF2....COL4A1....PTGS2....TERT....ABI.2....GNGT1....GNAS... <-  replace(path.data$c..TRAF2....COL4A1....PTGS2....TERT....ABI.2....GNGT1....GNAS...,values = combined.data$sn) 
  print(path.data$c..TRAF2....COL4A1....PTGS2....TERT....ABI.2....GNGT1....GNAS...)
}


View(path.data)


#Finding the sum of the genes with respect to the pathway, will yield the total value of pathway. 
# This will give us a robust way to rank the pathway, because higher the sum value of the pathway, means it contains genes which have greater expression value . 
  
#Was not able to clean the pathway file completely within the given time duration.
#Need to learn more of string manipulation.
  
  



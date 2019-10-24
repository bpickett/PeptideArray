library(reshape2)
library(ggplot2)

setwd("~/Desktop/PeptideArray")
mydata = read.csv(file="Summarized_Peptide_Results.tsv", header = TRUE, sep="\t")
#colnames(mydata) <- mydata[1, ] # the first row will be the header
mydata <- mydata[-1, ] #remove first label row
mydata <- mydata[-1, ] #remove second label row
mydata <- mydata[, -1] #remove first label column
rownames(mydata) <- mydata[,1] # the first column will be the rowName
colnames(mydata) <- gsub("^X", "",  colnames(mydata)) #remove leading X from all column headers
mydata <- mydata[, -1] #remove second label column
mydata[] <- lapply(mydata, function(x) as.numeric(as.character(x))) #convert factors to numeric
#dim(mydata)

##Normalize by peptide
NormZ <- as.data.frame(scale(mydata)) #converts values in data frame to Zscores ncol(mydata)
NormZ <- apply(NormZ,2,function(x) {ifelse(x < 0, 0, x)})#convert all negative numbers to "0", and uses data frame ##internal testing shows this increases OOB
NormZ1.p <- as.data.frame(NormZ) #converts values in data frame to Zscores ncol(mydata)
#ZColumns <- as.character(colnames(NormZ1.p))
#i <- as.vector(colSums(NormZ1.p >= 1.645))# find number of rows in each column with values >= x from mean. #1.645 is the 0.05 cutoff value for a one-tailed test

#Normalize by sample
#mydata.t <- t(mydata) #transposes matrix
#NormZ.t <- as.data.frame(scale(mydata.t)) #converts values in data frame to Zscores ncol(mydata.t)
#mydata.tN <- t(NormZ.t) #transposes matrix back
#NormZ1.s <- as.data.frame(apply(mydata.tN,2,function(x) {ifelse(x < 0, 0, x)}))#convert all negative numbers to "0", and uses data frame##internal testing shows this increases OOB
###mydata.tN1 <- as.data.frame(mydata.tN)
##ZColumns <- as.character(colnames(NormZ1.s))
##i <- as.vector(colSums(NormZ1.s >= 1.645))# find number of rows in each column with values >= x from mean

#Filter
#combined <- data.frame(ZColumns,i)
##cutoffValue <- round((nrow(NormZ)*.05)/1)#finds number representing 5% of samples (rounded to whole number)
#cutoffValue <- 1
##combined.subset <- subset(combined, combined$i>=cutoffValue)#finds columns above cutoff value (to be kept)
#rows2keep <- as.character(combined.subset$ZColumns)
##rows2keep

s1 = 1 ; e1 = 25 ; s2 = 516 ; e2 = 530#DENV1
s1 = 26 ; e1 = 76 ; s2 = 531 ; e2 = 549#DENV2
s1 = 77 ; e1 = 104 ; s2 = 550 ; e2 = 564#DENV3
s1 = 105 ; e1 = 151 ; s2 = 565 ; e2 = 594#DENV4
s1 = 256 ; e1 = 289 ; s2 = 685 ; e2 = 707#WNV
s1 = 491 ; e1 = 515 ; s2 = 0 ; e2 = 0#CHIKV
s1 = 421 ; e1 = 490 ; s2 = 797 ; e2 = 866 #ZIKV

###Choose 1(
##filter after peptide normalization
#NormZ.filtered.p <- NormZ1.p[colnames(NormZ1.p) %in% rows2keep]
#View(NormZ.filtered.p)
#write.table(NormZ.filtered.p, file = "DataMining-export-ZscoresFiltered.tsv", sep = "\t")
NormZ.filtered.p <- NormZ1.p[c(s1:e1,s2:e2)] #keep only columns with all data from specified virus species

#filter after sample normalization
#NormZ.filtered <- mydata.tN1[colnames(mydata.tN1) %in% rows2keep]

#NormZ.filtered.s <- NormZ1.s[c(s1:e1,s2:e2)] #keep only columns with all data from specified virus species
#write.table(NormZ.filtered.s, file = "zzDENV1-ZscoresUnfilteredSample.tsv", sep = "\t")
#write.table(data.frame("Species"=rownames(NormZ.filtered.s),NormZ.filtered.s), "zzZIKV-ZscoresUnfilteredSample.tsv", row.names=FALSE, sep = "\t")
write.table(data.frame("Species"=rownames(NormZ.filtered.p),NormZ.filtered.p), "zzZIKV-ZscoresUnfilteredPeptide.tsv", row.names=FALSE, sep = "\t")
#)

summary(NormZ.filtered.p)
#plot distribution of each column into one
#d <- melt(mydata[,1:50])
#d <- melt(NormZ[,1:50])
#d <- melt(NormZ.filtered[,1:50])
#d <- melt(NormZ.filtered[,1:ncol(NormZ.filtered)])
#ggplot(d,aes(x = value)) + 
#  facet_wrap(~variable,scales = "free_x") + 
#  geom_histogram(binwidth=.5)

#NormZ.dist <- as.matrix(dist(NormZ.filtered, method = "euclidean", diag = FALSE, upper = FALSE))
#plot(hclust(dist(NormZ.filtered),method="average"))
fit.p <- dist(NormZ.filtered.p, method="euclidean")#, labels = NormZ.filtered$column)#"average" = UPGMA 
fit <- hclust(dist(NormZ.filtered, method="euclidean"),method="average")#, labels = NormZ.filtered$column)#"average" = UPGMA 
fit <- hclust(dist(NormZ.filtered, method="euclidean"),method="ward.D2")#, labels = NormZ.filtered$column)#"average" = UPGMA 
fit <- hclust(dist(NormZ.filtered, method="maximum"),method="average")
fit <- hclust(dist(NormZ.filtered, method="maximum"),method="ward.D2")
fit <- hclust(dist(NormZ.filtered, method="manhattan"),method="average")
fit <- hclust(dist(NormZ.filtered, method="manhattan"),method="ward.D2")
fit <- hclust(dist(NormZ.filtered, method="canberra"),method="average")
fit <- hclust(dist(NormZ.filtered, method="canberra"),method="ward.D2")
fit <- hclust(dist(NormZ.filtered, method="minkowski"),method="average")
fit <- hclust(dist(NormZ.filtered, method="minkowski"),method="ward.D2")
plot(fit)
numClusters <- 4
groups <- cutree(fit, k=numClusters) # cut tree into k clusters
rect.hclust(fit, k=numClusters, border="red") # draw dendogram with red borders around the k clusters

#Visualize with MDS
MDS_results <- cmdscale(fit.p, eig = TRUE, k = 2)#MDS on distance
x <- MDS_results$points[, 1]
y <- MDS_results$points[, 2]
MDS_results1 <- do.call(rbind, Map(data.frame, X=x, Y=y))
wss <- (nrow(MDS_results1)-1)*sum(apply(MDS_results1,2,var))
j <- nrow(MDS_results1)-1
for (i in 2:j) wss[i] <- sum(kmeans(MDS_results1, centers=i)$withinss)#need to modify second integer to be (listLength-1)
k <- length(wss)
plot(1:k, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")#need to modify second integer to be (listLength-1)
# K-Means Cluster Analysis
fit_k <- kmeans(MDS_results1, 3) # kmeans on specified # of clusters
aggregate(MDS_results1,by=list(fit_k$cluster),FUN=mean)
MDSresults11 <- data.frame(MDS_results1, fit_k$cluster)
par(mfrow=c(1,1), mar=c(3, 3, 3, 7), xpd=TRUE)
#ggplot(MDS_results1, aes(x=MDS_results1$X, y=MDS_results1$Y)) 
#plot(MDSresults11$X, MDSresults11$Y, pch=21, bg=groups, col="black")
plot(MDSresults11$X, MDSresults11$Y, col=(fit_k$cluster +1), main="K-Means 1000 points", pch=20, cex=2)

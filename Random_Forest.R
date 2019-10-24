#taken from http://www.listendata.com/2014/11/random-forest-with-r.html
# Random Forest prediction of Kyphosis data
library(randomForest)
#library(inTrees)
#library(ROCR)
#library(tidyr)

setwd("~/Desktop/PeptideArray")
mydata = read.csv(file="zzZIKV-ZscoresUnfilteredPeptide.tsv", header = TRUE, sep="\t")
#mydata <- as.data.frame(apply(mydata,2,function(x) {ifelse(x < 0, 0, x)}))
#write.table(mydata, file = "ZIKV119-ZscoresFilteredSample0.tsv", sep = "\t")

#mydata = read.csv(file="CHIKV_Peptide_rfInput.tsv", header = TRUE, sep="\t")
#data(iris)
#nrow(mydata)
#summary(mydata)
#mydata <- iris
#mydata.na <- mydata
set.seed(111)
## artificially drop some data values.
#for (i in 1:4) mydata.na[sample(150, sample(20)), i] <- NA
#set.seed(222)
#mydata.imputed <- rfImpute(Species ~ ., mydata.na)
#set.seed(335)
#rf.myData <- varSelRF(mydata, cl, ntree = 100000, ntreeIterat = 10000, vars.drop.frac = 0.2)
numTry <- sqrt(ncol(mydata))
mydata.rf <- randomForest(Species ~ ., data=mydata, ntree=100000, importance = TRUE, do.trace = 10000, proximity = TRUE, mtry = numTry)#, do.trace = 10000, mtry = 10) 
importance.df <- importance(mydata.rf, scale = TRUE)#per Strobl and Zeileis (2008), who reported that a scaling approach can lead to undesirable results.
print(mydata.rf)
print(importance.df)
varImpPlot(mydata.rf)
#used <- varUsed(mydata.rf, by.tree=FALSE, count=TRUE)
#combinedVars <- data.frame(used,importance.df)
#mtry <- tuneRF(mydata.rf[-5],mydata.rf$Species, ntreeTry=100000,
#               stepFactor=1.5,improve=0.01, trace=TRUE, plot=TRUE)
#best.m <- mtry[mtry[, 2] == min(mtry[, 2]), 1]
#print(mtry)
#print(best.m)

# #Evaluate the performance of the random forest for classification.
#pred2=predict(mydata.rf,type = "prob")
prediction <- predict(mydata.rf,mydata,type="prob")
votes <- 100000*(predict(mydata.rf,mydata,type="vote",norm.votes=TRUE))#converts prediction (percentages) into number of votes
rfDecision <- getTree(randomForest(mydata[,-5], mydata[,5], ntree=10), 3, labelVar=TRUE)
print(rfDecision)
View(rfDecision)
#data("ROCR.simple")
#View(ROCR.simple)
# #prediction is ROCR function
#perf = prediction(pred2[,2], mydata$Species)
# 
# #performance in terms of true and false positive rates
# #1. Area under curve
#auc = performance(perf, "auc")
# #2. True Positive and Negative Rate
# pred3 = performance(perf, "tpr","fpr")
# #3. Plot the ROC curve
# plot(pred3,main="ROC Curve for Random Forest",col=2,lwd=2)
# abline(a=0,b=1,lwd=2,lty=2,col="gray")
# 
# #Plot a sample tree
# cforest(Species~., data=mydata, controls=cforest_control(mtry=2, mincriterion=0))

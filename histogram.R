fcq1 <- read.table("Maize_RPKM_nogenelen.txt_SP.txt_random-100000-99-0.662.txt",header=T,sep="\t")

#hist(log(fcq1$Fruits, base=10))
#hist(sqrt(fcq1$Fruits))

## subset data
randompccs <- fcq1$random_pcc.99.percentile.0.662
## omit NAs
randompccs <- na.omit(randompccs)
#randompccs <- log10(randompccs)
## get 5% and 95% thresholds
outputQuantile <- quantile(randompccs, c(.01,.99)) #duration is a column, and then you get the 5% and 95% of the distribution
## get upper and lower percentiles
dfQuantile <- cbind(outputQuantile)
dfQuantile <- data.frame(dfQuantile)
upperperc <- dfQuantile$outputQuantile[[2]]
lowerperc <- dfQuantile$outputQuantile[[1]]
## put data in dataframe indicating data above upperpercent
dat <- data.frame(x=randompccs, above=randompccs>=upperperc) #make 95th percentile different color
library(ggplot2)
## make histogram plot with thresholds
qplot(x,data=dat,geom="histogram",fill=above, binwidth=0.01, xlab="Spearman's rank", main="Random expectation- 100000 draws, 99 percentile")
dat2 <- data.frame( x=randompccs, below=randompccs<=lowerperc) #make percentile different color
qplot(x,data=dat2,geom="histogram",fill=below, binwidth=0.01, xlab="PCC", main="Random expectation")


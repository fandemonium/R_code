library(plyr)
library(BayesFactor)

args <- commandArgs(TRUE)

si <- readRDS(args[1])

temp<-si[, -c(1, 3:21, 110, 114, 116:120, 123)]
for (i in colnames(temp)){
	sink(args[2], append = TRUE)
	tryCatch({
	print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
	test<-temp[, c(i, "myear", "foam.type", "id")]
	colnames(test)[1]<-"to_test"
	test<-test[! is.na(test$to_test), ]
	bf1<-anovaBF(to_test ~ foam.type + myear + id, data=test, whichRandom=c("myear", "id"))
	average<-ddply(test, .(foam.type), summarize, avg=mean(to_test))
	average<-average[order(average$avg), ]
	samples = posterior(bf1, iterations = 10000)
	consistent = (samples[, paste("foam.type", average[3,1], sep="-")] > samples[, paste("foam.type", average[2,1], sep="-")]) & (samples[, paste("foam.type", average[2,1], sep="-")] > samples[, paste("foam.type", average[1,1], sep="-")])
	N_consistent = sum(consistent)
	bf_restriction_against_full = (N_consistent / 10000) / (1 / 6)
	bf_full_against_null = as.vector(bf1)
	bf_restriction_against_null = bf_restriction_against_full * bf_full_against_null
	print(bf_restriction_against_null)
	}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
	sink()
}


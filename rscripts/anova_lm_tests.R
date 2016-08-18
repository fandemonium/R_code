### residual normality test ####
shapiro_tb<-data.frame()
for (i in names(fs[, c(7:18, 21)])){
 tryCatch({
 a1<-lm(log(fs[, i]+1) ~ fs$category)
 a1.stdres<-rstandard(a1)
 shapiro_test<-shapiro.test(a1.stdres)
 stat_param<-i
 if (shapiro_test$p.value > 0.05){
 	test<-data.frame(stat_param, shapiro_test$p.value)
 	shapiro_tb<-rbind(shapiro_tb, test)
 	}
  }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}	
### run anova on normal residual samples ###
anova_tb<-data.frame()
for (i in names(fs[, colnames(fs) %in% shapiro_tb$stat_param])){
 a1<-anova(lm(log(fs[, i]+1) ~ fs$category))
 id <- i
 p <- a1$"Pr(>F)"[1]
 if (p <= "0.05"){
 	test<-data.frame(id, p)
	anova_tb<-rbind(anova_tb, test)
	}
}
### creat box plot for all significant lm's ####
temp.df<-data.frame(fs[, c("category")], fs[, colnames(fs) %in% anova_tb$id])
colnames(temp.df)[1]<-"category"
temp.df<-melt(temp.df, id.vars="category")
ggplot(temp.df, aes(category, value)) + geom_boxplot()+ facet_grid(variable ~., scales="free")+theme_bw()
### post-hoc analysis ###
posthoc_tb<-data.frame()
for (i in names(fs[, colnames(fs) %in% anova_tb$id])){
	tryCatch({
	m1<-aov(log(fs[, i]+1) ~ fs$category)
	posthoc<-TukeyHSD(m1)
	test<-data.frame(posthoc$fs)
	sig<-test[test[,4]<=0.05,]
	sig$category<-i
	posthoc_tb<-rbind(posthoc_tb, sig)
	}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
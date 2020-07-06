# chow test 
# from pachevalier/chow_test.R (https://gist.github.com/pachevalier/6131725)
# see stackoverflow (https://stats.stackexchange.com/questions/66250/compare-2-regression-lines-in-r) for example. 

mc  <- lm(formula = SBR ~ Age, data = ch)
m1  <-  lm(formula = SBR ~ Age, data = subset(ch, Sex == "M"))
m2  <-  lm(formula = SBR ~ Age, data = subset(ch, Sex == "F"))
sc  <- sum(mc$residuals^2)
s1  <- sum(m1$residuals^2)
s2  <- sum(m2$residuals^2)
k  <- 2
# Test statistic
fstat  <- (sc - (s1 + s2)) / k / (s1 + s2) * (length(mc$residuals) - 2*k)  
fstat

# Rejection region
qf(.95,k, length(mc$residuals) - 2*k)

# Pvalue
pf(fstat,k, length(mc$residuals) - 2*k)

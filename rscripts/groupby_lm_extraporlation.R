##### example table 
####assay_name	assay_date	Time	temp	A1	A2	A3	A4
####FullPlate969_test	20200319	0	29.9	0.379	0.37	0.339	0.331
####FullPlate969_test	20200319	30	30	0.339	0.31	0.286	0.281
####FullPlate969_test	20200319	60	30	0.314	0.28	0.268	0.251


library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(reshape2)

# read in the sample data you sent me in the first email
data <- read.csv("~/Downloads/test.csv")
# melt
data.melt <- melt(data, id=c("assay_name", "assay_date", "Time", "temp"))

# get the slopes 
tbl <- data.melt %>% 
	group_by(variable) %>% #groupby plate id
	nest() %>% 
	mutate(model = map(data, ~lm(value ~ Time, data = .x) %>% tidy)) %>% #linear regression and look at the stats per term
	unnest(model) %>% 
	filter(term == 'Time') %>% #get time related stats only
	select(estimate) %>% #get slope only
	data.frame() #personal preference

# r2 is evaluated at the full model levell. so create a different table
r2 <- data.melt %>% 
	group_by(variable) %>% 
	nest() %>% 
	mutate(model = map(data, ~lm(value ~ Time, data = .x) %>% glance)) %>% #linear regression and look at the stats at full model level. 
	unnest(model) %>% 
	select(adj.r.squared) %>% #get r2
	data.frame() #personal perference

final <- cbind(tbl, r2)

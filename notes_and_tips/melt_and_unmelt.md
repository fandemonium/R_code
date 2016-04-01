library(reshape2)

#to convert data from wide formate data.frame to long format:
#reshape2:melt allows you to melt without naming id columns, which is nice when you every column is a variable.
melt(wide_data_frame)

#to convert from long to wide:
#if all of you factors have equal number of samples, then you can just do:
all.wide<-dcast(all, tree ~ type, value.var="Td")
#otherwise, R will have no idea what to do with the extra values of some factors. Then you will have to aggregate those values:
all.wide<-dcast(all, tree ~ type, value.var="Td", fun.aggregate=mean, na.rm=T)

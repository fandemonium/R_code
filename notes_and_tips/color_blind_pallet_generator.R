library(RColorBrewer)

# create a function that generates color palettes based on "Dark2"
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))

# determine the number of different colors needed, eg. 10
colors = getPallette(10)

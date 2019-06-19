#usage: Rscript get_neon_raw_tables.R "DP1.20086.001"

library(httr)
library(jsonlite)
library(dplyr, quietly=T)
library(downloader)
library(geoNEON)
library(neonUtilities)

args <- commandArgs(TRUE)

product_id <- args[1]
# product_id = "DP1.20086.001"

req <- GET(paste0("http://data.nenonscience.org/api/v0/products/", product_id))

#req.content <- content(req, as="parsed")

avail <- fromJSON(content(req, as="text"). simplifyDataFrame=T, flatten=T)

# list all urls contain all sites and dates
data.urls <- unlist(avail$data$siteCodes$availableDataUrls)

to_download <- data.frame()
for (i in seq(1, length(data.urls))){
      data <- GET(data.urls[i])
      data.file <- fromJSON(content(data, as="text"))
      temp <- data.file$data$files %>% filter(grepl("rawdata", ignore.case=T, name)) %>% filter(grepl("expanded", name)) %>% select(name, url)
      to_download <- rbind(to_download, temp)
}

out <- paste0("curl -o ", to_download$name, " ", to_download$url)

write.table(out, paste0(product_id, "_curl_command.sh"), quote=F, row.names=F, col.names=F)

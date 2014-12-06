getwd()
setwd("~/Dropbox/Public/DomConsAirfare CSV.csv")
data <- read.table("DomConsAirfare CSV",header=TRUE,sep=",")
data
(cl <- kmeans(data, 2))
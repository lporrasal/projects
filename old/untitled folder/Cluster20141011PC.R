getwd()
setwd("C:/Users/salilmehta/Desktop")
data <- read.table("Cluster.csv",header=TRUE,sep=",")
data
(cl <- kmeans(data, 2))
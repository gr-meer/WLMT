library(data.table)

names6m <- read.table(list_of_samples_of_interest, stringsAsFactors = FALSE, header = FALSE)

print(system.time({  data0 <- fread(paste0("zcat ", names6m[1,1]),stringsAsFactors = FALSE) }))
samplename <- substr(names6m[1,1],8,str_length(names6m[1,1])-10)
data0 <- data0[,c("index","coverage")]

colnames(data0)
colnames(data0)[2:ncol(data0)] <-seq(2:ncol(data0))

for(i in 2:nrow(names6m))
{
  data1 <- fread(paste0("zcat ", names6m[i,1]), header = TRUE, stringsAsFactors = FALSE)
  data1 <- data.frame(data1)
  data1 <- data1[,c("index","coverage")]
  colnames(data1)[2]<- paste0("coverage",i+1000)
  data0 <- merge(data0, data1, all = FALSE, by="index")
  samplename <- substr(names6m[i,1],8,str_length(names6m[i,1])-10)
  print(paste0("File ", i, " finished processing"))
  print(Sys.time())
}

d0_covered_5_90 <- data0[apply(data0[,-1],1,function(x){sum(x>=5)/(ncol(data0)-1)>=0.90}),]

write.csv(list_of_sites_d0_covered_5_90,list_of_sites_with_good_coverage)
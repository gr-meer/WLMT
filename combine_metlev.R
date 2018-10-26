library(data.table)

#list_of_samples_of_interest contains paths of the *.metcov.gz files

names <- read.table(list_of_samples_of_interest, stringsAsFactors = FALSE, header = FALSE)
#all_suitable_sites - a dataframe where column #2 is a list of all sites which pass coverage requirement
selected_sites <- read.csv( all_suitable_sites, header = TRUE, stringsAsFactors = FALSE)
selected_sites <- selected_sites[,2]
selected_sites <- data.frame(selected_sites)
colnames(selected_sites)[1] <- "index"

data0 <- selected_sites
for(i in 1:nrow(names))
{
  data1 <- fread(paste0("zcat ", names[i,1]), header = TRUE, stringsAsFactors = FALSE)
  data1 <- data.frame(data1)
  data1 <- data1[,c("index","metlev")]
  colnames(data1)[2]<- paste0("metlev",i)
  data0 <- merge(data0, data1, all = FALSE, by="index")
  samplename <- substr(names[i,1],8,str_length(names[i,1])-10)
  print(paste0("File ", i, " finished processing"))
  print(Sys.time())
}

#autosomes_suitable_sites - a dataframe where column #1 is a list of only the sites which pass coverage requirement and also located on autosomes
autosomes_sites <- read.csv(autosomes_suitable_sites, header = FALSE, stringsAsFactors = FALSE)
autosomes_sites <- data.frame(autosomes_sites)
colnames(autosomes_sites)[1] <- "index"
everything <- data0
everything <- merge(autosomes_sites, data0, all = FALSE, by="index")
head(everything)

#everything <- everything[-c(1)]
everything <- t(everything)

#resulting_table - name of the final dataframe which can be used to construct or to apply Whole Lifespan Multi-Tissue Clock
write.table(everything, resulting_table,col.names = FALSE,row.names = FALSE,sep = '\t')

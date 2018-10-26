# This script converts files produced by Bismark v0.15.0 to a format with following columns:
# 1. Chormosome name
# 2. Unique CpG site index which is based on chromosome number and the position on this chromosome
# 3. Coverage
# 4. Methylation level

# Example: 
# chr18,chr18_3101849,133,83.4586466165414


library(data.table)
library(stringr)

#locations_of_samples_of_interest is a table where the first row contains name of a sample 
# and the second - path to the corresponding *.bismark.cov.gz file produced by Bismark v0.15.0

names <- read.table(locations_of_samples_of_interest, stringsAsFactors = FALSE, header = TRUE)
for(i in 1:nrow(names))
{
  print(i)
  print(names[i,1])
  print(Sys.time())
  if(str_sub(names[i,2], start= -3) != '.gz')
  {
    print(system.time({  data0 <- fread(names[i,2],header = FALSE, stringsAsFactors = FALSE) }))
  }
  else
  {print(system.time({  data0 <- fread(paste0("zcat ", names[i,2])) }))}
  print(system.time({  
    index <- paste0(data0[[1]],"_",data0[[2]])
    coverage <- data0[[5]] + data0[[6]]
    data1 <- data.frame("chr" = data0[[1]], "index" = index, "coverage" = coverage, "metlev" = data0[[4]], stringsAsFactors = FALSE)
    outfilename <- paste0("metcov/",names[i,1],".metcov")
    fwrite(data1,file = outfilename,col.names = TRUE)
    system(paste0("pigz -p 8 ",outfilename))
  }))
}




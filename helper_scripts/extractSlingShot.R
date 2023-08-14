library(SingleCellExperiment)
library(slingshot)
library(RColorBrewer)

data = readRDS('path/to/SingleCellExperiment/rds/object')

curve = slingCurves(data)

lineage_counter = 1

for (x in curve){
  
  lineage_title = paste("Lineage",as.character(lineage_counter), sep = "_")
  file_title = paste(lineage_title,".csv", sep = "")
  
  lineage_counter = lineage_counter + 1
  
  lineage = x$s
  
  write.csv(lineage,file = file_title, row.names = F)
  
}


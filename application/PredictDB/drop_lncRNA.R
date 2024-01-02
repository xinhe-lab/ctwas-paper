library(RSQLite)
library(tools)

weight_dir <- "/project2/compbio/predictdb/mashr_models/"
output_dir <- "/project2/mstephens/wcrouse/predictdb_nolnc/"

weights <- list.files(weight_dir)
weights <- weights[grep(".db", weights)]
weights <- paste0(weight_dir, weights)

for (i in 1:length(weights)){
  print(i)
  
  weight <- weights[i]
  weight_stem <- rev(unlist(strsplit(file_path_sans_ext(weight), "/")))[1]
  
  sqlite <- RSQLite::dbDriver("SQLite")
  db = RSQLite::dbConnect(sqlite, weight)
  query <- function(...) RSQLite::dbGetQuery(db, ...)
  weights_table <- query("select * from weights")
  extra_table <- query("select * from extra")
  dbDisconnect(db)
  
  #subset to protein coding genes only
  extra_table <-  extra_table[extra_table$gene_type=="protein_coding",,drop=F]
  weights_table <- weights_table[weights_table$gene %in% extra_table$gene,]
  
  weight_info = read.table(gzfile(paste0(file_path_sans_ext(weight), ".txt.gz")), header = T)
  weight_info <- weight_info[weight_info$GENE %in% extra_table$gene,]
  
  #write file
  if (!file.exists(paste0(output_dir, weight_stem, "_nolnc.db"))){
    db <- dbConnect(RSQLite::SQLite(), paste0(output_dir, weight_stem, "_nolnc.db"))
    dbWriteTable(db, "extra", extra_table)
    dbWriteTable(db, "weights", weights_table)
    dbDisconnect(db)
    
    weight_info_gz <- gzfile(paste0(output_dir, weight_stem, "_nolnc.txt.gz"), "w")
    write.table(weight_info, weight_info_gz, sep=" ", quote=F, row.names=F, col.names=T)
    close(weight_info_gz)
  }
}

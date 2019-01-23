# load package
library(rtracklayer)
library(RIPSeeker)
library(chipseq)

# load function
prepareChipseq <- function(reads){
  frag_len = median(estimate.mean.fraglen(reads))
  cat(paste0('median size is ', round(frag_len)))
  reads_extended = resize(reads, width = frag_len)
  return(trim(reads_extended))
}

# load bed file
path <- "/mnt/c/chipseq/test3"
file_list <- list.files(path)
for(file_name in file_list){
    file_split_1 <- strsplit(file_name, split = "_", fixed = TRUE)[[1]]
    file_split_2 <- strsplit(file_name, split = ".", fixed = TRUE)[[1]]
    if(length(file_split_1) == 3 & file_split_1[3] == "chr6.bed"){
        file_import <- import.bed(file.path(path, file_name))
        file_export <- prepareChipseq(file_import)
        exportGRanges(file_export, paste0(path, "/", file_split_2[1], ".txt"), exportFormat = "txt")
    }
}
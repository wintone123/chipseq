# load package
library(rtracklayer)
library(RIPSeeker)
library(chipseq)

# load function
prepareChipseq <- function(reads){
  frag_len = median(estimate.mean.fraglen(reads))
  reads_extended = resize(reads, width = frag_len)
  return(trim(reads_extended))
}

# load bed file
path <- "/mnt/c/chipseq/test3"
file_list <- list.files(path)
for(file_name in file_list){
    file_split <- strsplit(file_name, split = "_", fixed = TRUE)[[1]]
    if(length(file_split) == 2 & file_split[2] == "best.bed"){
        print(paste0("extending ", file_name, "......"))
        file_import <- import.bed(file.path(path, file_name))
        file_import <- file_import[file_import@seqnames %in% paste0("chr", c(1:19, "X", "Y"))]
        file_export <- prepareChipseq(file_import)
        exportGRanges(file_export, paste0(path, "/", file_split[1], "_best_ext.txt"), exportFormat = "txt")
    }
}
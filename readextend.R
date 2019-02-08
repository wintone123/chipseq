# load package
library(rtracklayer)
library(RIPSeeker)
library(chipseq)

# load function
PrepareChipseq <- function(reads){
  frag_len = median(estimate.mean.fraglen(reads))
  reads_extended = resize(reads, width = frag_len)
  return(trim(reads_extended))
}

# load bed file
path <- "/mnt/c/chipseq/test3"
dir.create(file.path(path, "ext_split_files"))
file_list <- list.files(file.path(path, "split_files"))
for (file_name in file_list) {
    file_split <- strsplit(file_name, split = ".", fixed = TRUE)[[1]]
    print(paste0("extending ", file_name, "......"))
    file_import <- import.bed(file.path(path, "split_files", file_name))
    file_export <- prepareChipseq(file_import)
    exportGRanges(file_export, paste0(path, "/ext_split_files/", file_split[1], "_ext.txt"), exportFormat = "txt")
}
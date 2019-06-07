# load library
library(Biostrings)

# load necessary data

path <- "/mnt/c/chipseq/test2"  # input file path
load_list <- list.files(file.path(path, "split_files"))
load_list <- load_list[which(load_list != "output")]
dir.create(file.path(path,"output"), showWarnings = FALSE)  # create output folder

# genome_list
genome_list <- c(10:19, "X","Y")
print("=============Let's Go!=============")
# analysis process
for(chrom_1 in genome_list){
  print(chrom_1)
  chrom <- paste0("chr", chrom_1)
  print(chrom)
  read_list <- vector()
  for(file_name in load_list){
    if(strsplit(file_name, split = "_", fixed = TRUE)[[1]][3] == paste0(chrom, ".bed"))
      base_name <- strsplit(file_name, split = ".", fixed = TRUE)[[1]][1]
      read_list <- c(read_list, base_name)
    }
  print(read_list)
}

# load packages
library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(chipseq)
library(biomaRt)
library(dplyr)

print("==========Let's go!==========")
print("loading data......")
# load function
BinChipseq <- function(reads, bins) {
  mcols(bins)$score = countOverlaps(bins, reads)
  return(bins)
}
TileCollection <- function (collection, bin_num) {
    for (i in 1:nrow(collection)) {
        start <- collection$start_position[i]
        end <- collection$end_position[i]
        seq_list <- seq(start, end, length.out = bin_num + 1)
        start_list <- vector(length = bin_num)
        end_list <- vector(length =  bin_num)
        for (j in 1:(bin_num + 1)) {
            if (j == 1) {
                start_list[j] <- seq_list[j]
            } else if (j == bin_num + 1) {
                end_list[j-1] <- seq_list[j]
            } else {
                start_list[j] <- seq_list[j] + 1
                end_list[j-1] <- seq_list[j]
            }
        }
        if (i == 1) {
            data_df <- data.frame(seqname = rep(collection$chromosome_name[i], bin_num),
                                  start = start_list, 
                                  end = end_list,
                                  strand = rep(collection$strand[i], bin_num))
        } else {
            data_df_temp <- data.frame(seqname = rep(collection$chromosome_name[i], bin_num),
                                       start = start_list, 
                                       end = end_list,
                                       strand = rep(collection$strand[i], bin_num))
            data_df <- rbind(data_df, data_df_temp)
        }
    }
    GRanges_temp <- GRanges(seqnames = Rle(data_df$seqname),
                            ranges = IRanges(start = data_df$start, end = data_df$end),
                            strand = Rle(data_df$strand))
    return(GRanges_temp)
}
PrepareChipseq <- function (reads) {
  frag_len = median(estimate.mean.fraglen(reads))
  reads_extended = resize(reads, width = frag_len)
  return(trim(reads_extended))
}

# parameter set
bed_path <- "/mnt/c/chipseq/test2/split_files"
path <- "/mnt/c/chipseq/test7/"
cpg_file <-"cpg_mm10_2.bed"
chrom_list <- c(1:19,"X","Y")
bin_num <- 60

# load file
cpg_mm10 <- import.bed(file.path(path, cpg_file))
bed_file_list <- list.files(bed_path)

# main
print("+++++++++loop start+++++++++")
for (num in 1:length(chrom_list)){
    chrom <- chrom_list[num]
    chrom_name <- paste0("chr", chrom)
    print(paste0("analysis ", chrom_name, "......"))

    # load bed file
    for (file in bed_file_list) {
        name_split <- strsplit(file, split = "_", fixed = TRUE)[[1]]
        if (name_split[1] == "rep1" & name_split[3] == paste0(chrom_name, ".bed")) {
            chrom_bed <- PrepareChipseq(import.bed(file.path(bed_path, file)))
        }
    }
    length_bed <- length(chrom_bed) / 1000000

    # import gene info
    cpg <- cpg_mm10[cpg_mm10@seqnames == chrom_name]
    start(cpg) <- start(cpg) - width(cpg)
    end(cpg) <- end(cpg) + width(cpg)
    cpg_df <- data.frame(chromosome_name = cpg@seqnames,
                         start_position = start(cpg),
                         end_position = end(cpg),
                         strand = cpg@strand,
                         stringsAsFactors = FALSE)
    cpg_bins <- TileCollection(cpg_df, bin_num)
    
    # overlap
    bin_cpg <- BinChipseq(chrom_bed, cpg_bins)
    data_df <- data.frame(chrom = rep(NA, nrow(cpg_df)*bin_num),
                          position = rep(NA, nrow(cpg_df)*bin_num),
                          Hits = rep(NA, nrow(cpg_df)*bin_num))
    for (i in 1:nrow(cpg_df)) {
        data_df[(bin_num*(i-1)+1):(bin_num*i),]$chrom <- rep(chrom_name, bin_num)
        data_df[(bin_num*(i-1)+1):(bin_num*i),]$position <- c(1:bin_num)
        data_df[(bin_num*(i-1)+1):(bin_num*i),]$Hits <- bin_cpg$score[((i-1)*bin_num+1):(i*bin_num)] / length_bed
    }

    if (num == 1) {
        score <- data_df
    } else {
        score <- rbind(score, data_df)
    }

    # clean memory
    # for (name in bed_list) {
    #     rm(name)
    # }
    # for (name in peak_list) {
    #     rm(name)
    # }
}

# output csv
print(paste0("writing h3k27_on_cpg_2......"))
write.csv(score, file = file.path(path, "h3k27_on_cpg_2.csv"))
print("===========Finish!===========")
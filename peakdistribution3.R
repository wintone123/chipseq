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
PrepareChipseq <- function(reads){
  frag_len = median(estimate.mean.fraglen(reads))
  reads_extended = resize(reads, width = frag_len)
  return(trim(reads_extended))
}
BinChipseq <- function(reads, bins){
  mcols(bins)$score = countOverlaps(bins, reads)
  return(bins)
}
ImportNp <- function(path){
  temp <- read.delim2(path, quote = "/", head = FALSE)
  GRanges_temp <- GRanges(seqnames = Rle(temp[,1]),
                          ranges = IRanges(start = temp[,2], end = temp[,3]),
                          strand = Rle("*"),
                          name = as.character(temp[,4]),
                          score = as.numeric(temp[,5]))
  return(GRanges_temp)
}
OverlapPeaks <- function(GRanges_1, GRanges_2){
    df_1 <- data.frame(start(GRanges_1), end(GRanges_1))
    df_2 <- data.frame(start(GRanges_2), end(GRanges_2))
    start_out <- vector()
    end_out <- vector()
    loop_s <- 1
    for (i in 1:nrow(df_1)) {
        start_1 <- df_1[i,1]
        end_1 <- df_1[i,2]
        for (j in loop_s:nrow(df_2)) {
            start_2 <- df_2[j,1]
            end_2 <- df_2[j,2]
            if (end_2 <= start_1) {
                loop_s <- j
            } else if (start_2 >= end_1) {
                break
            } else {
                if (start_2 <= start_1) {
                    start_out <- append(start_out, start_1)
                    if (end_2 <= end_1 & end_2 > start_1) {
                        end_out <- append(end_out, end_2)
                    } else if (end_2 > end_1) {
                        end_out <- append(end_out, end_1)
                    }
                } else if (start_2 > start_1 & start_2 < end_1) {
                    start_out <- append(start_out, start_2)
                    if (end_2 <= end_1) {
                        end_out <- append(end_out, end_2)
                    } else {
                        end_out <- append(end_out, end_1)
                    }
                }
            }
        }
    }
    new_GRanges <- GRanges(seqnames = Rle(GRanges_1@seqnames[1]),
                           ranges = IRanges(start = start_out, end = end_out),
                           strand = Rle("*"))
    return(new_GRanges)
}

# parameter set
path <- "/mnt/c/chipseq/test2"
bed_file_list <- list.files(file.path(path, "split_files"))
peak_file_list <- list.files(file.path(path, "narrowpeak_split_files"))
chrom_list <- c(1:19, "X", "Y")
up_position <- 1000
down_position <- 1000
bin_num <- 30
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "mmusculus_gene_ensembl",  
                host = "jul2018.archive.ensembl.org")  
ds <- useDataset("mmusculus_gene_ensembl", mart = mart)

# main
print("+++++++++loop start+++++++++")
for (num in 1:length(chrom_list)){
    chrom <- chrom_list[num]
    chrom_name <- paste0("chr", chrom)
    print(paste0("analysis ", chrom_name, "......"))
    bed_list <- vector()
    peak_list <- vector()

    # load bed file
    for (file in bed_file_list) {
        name_split <- strsplit(file, split = "_", fixed = TRUE)[[1]]
        if (name_split[3] == paste0(chrom_name, ".bed")) {
            assign(name_split[1], PrepareChipseq(import.bed(file.path(path, "split_files", file))))
            if (name_split[1] != "control") {
                bed_list <- append(bed_list, name_split[1])
            }
        }
    }

    # load peak file
    for (file in peak_file_list) {
        name_split <- strsplit(file, split = "_", fixed = TRUE)[[1]]
        if (name_split[4] == paste0(chrom_name, ".narrowPeak")) {
            assign(paste0(name_split[1], "_peak"), ImportNp(file.path(path, "narrowpeak_split_files", file)))
            peak_list <- append(peak_list, paste0(name_split[1], "_peak"))
        }
    }
    for (i in 2:length(peak_list)) {
        if (i <= 2) {
            overlap_peak <- OverlapPeaks(eval(parse(text = peak_list[1])), eval(parse(text = peak_list[2])))
        } else {
            overlap_peak <- OverlapPeaks(overlap_peak, eval(parse(text = peak_list[i])))
        }
    }

    # import gene info
    gene_info <- getBM(attributes = c("ensembl_gene_id","external_gene_name",
                                      "chromosome_name", "start_position",
                                      "end_position", "strand", "gene_biotype"), 
                       mart = ds, filter = "chromosome_name", values = chrom)
    gene_info <- distinct(gene_info, external_gene_name, .keep_all = TRUE)
    gene_info_G <- GRanges(seqnames = Rle(chrom_name),
                           ranges = IRanges(start = gene_info$start_position, end = gene_info$end_position),
                           strand = Rle("*"))
    gene_info_fil <- gene_info[unique(queryHits(findOverlaps(gene_info_G, overlap_peak))),]
    for (i in 1:nrow(gene_info_fil)) {
        start <- gene_info_fil$start_position[i] - up_position
        end <- gene_info_fil$start_position[i] + down_position
        seq_list <- seq(start, end, length.out = bin_num+1)
        start_list <- vector(length = length(seq_list)-1)
        end_list <- vector(length =  length(seq_list)-1)
        for (j in 1:length(seq_list)){
            if (j == 1) {
                start_list[j] <- seq_list[j]
            } else if (j == length(seq_list)) {
                end_list[j-1] <- seq_list[j]
            } else {
                start_list[j] <- seq_list[j] + 1
                end_list[j-1] <- seq_list[j]
            }
        }
        if (i == 1) {
            gene_info_ext <- data.frame(name = rep(gene_info_fil$external_gene_name[i], bin_num), start = start_list, end = end_list, strand = rep(gene_info_fil$strand[i], bin_num))
        } else {
            gene_info_ext_temp <- data.frame(name = rep(gene_info_fil$external_gene_name[i], bin_num), start = start_list, end = end_list, strand = rep(gene_info_fil$strand[i], bin_num))
            gene_info_ext <- rbind(gene_info_ext, gene_info_ext_temp)
        }
    }
    gene_info_ext_G <- GRanges(seqnames = Rle(chrom_name),
                               ranges = IRanges(start = gene_info_ext$start, end = gene_info_ext$end),
                               strand = Rle(gene_info_ext$strand),
                               name = Rle(gene_info_ext$name))
    
    # overlap
    for (i in 1:length(bed_list)) {
        if (i == 1) {
            bin_test <- BinChipseq(eval(parse(text = bed_list[i])), gene_info_ext_G)
        } else {
            bin_test_temp <- BinChipseq(eval(parse(text = bed_list[i])), gene_info_ext_G)
            bin_test$score <- bin_test$score + bin_test_temp$score
        }
    }
    bin_test$score <- bin_test$score/length(bed_list)
    bin_control <- BinChipseq(eval(parse(text = "control")), gene_info_ext_G)
    data_df <- data.frame(name = rep(NA, nrow(gene_info_fil)*bin_num*2),
                          chrom = rep(NA, nrow(gene_info_fil)*bin_num*2),
                          item = rep(NA, nrow(gene_info_fil)*bin_num*2),
                          position = rep(NA, nrow(gene_info_fil)*bin_num*2),
                          Hits = rep(NA, nrow(gene_info_fil)*bin_num*2))
    for (i in 1:nrow(gene_info_fil)) {
        print(gene_info_fil$external_gene_name[i])
        data_df[(bin_num*2*i-bin_num*2+1):(bin_num*2*i),]$name <- rep(gene_info_fil$external_gene_name[i], bin_num*2)
        data_df[(bin_num*2*i-bin_num*2+1):(bin_num*2*i),]$chrom <- rep(chrom_name, bin_num*2)
        data_df[(bin_num*2*i-bin_num*2+1):(bin_num*2*i-bin_num),]$item <- rep("test", bin_num)
        data_df[(bin_num*2*i-bin_num+1):(bin_num*2*i),]$item <- rep("control", bin_num)
        data_df[(bin_num*2*i-bin_num*2+1):(bin_num*2*i-bin_num),]$position <- c(1:bin_num)
        data_df[(bin_num*2*i-bin_num+1):(bin_num*2*i),]$position <- c(1:bin_num)
        if (as.logical(bin_test@strand[(i-1)*bin_num+1] == "+")) {
            data_df[(bin_num*2*i-bin_num*2+1):(bin_num*2*i-bin_num),]$Hits <- bin_test$score[((i-1)*bin_num+1):(i*bin_num)]
        } else {
            data_df[(bin_num*2*i-bin_num*2+1):(bin_num*2*i-bin_num),]$Hits <- rev(bin_test$score[((i-1)*bin_num+1):(i*bin_num)])
        }
        if (as.logical(bin_control@strand[(i-1)*bin_num+1] == "+")) {
            data_df[(bin_num*2*i-bin_num+1):(bin_num*2*i),]$Hits <- bin_control$score[((i-1)*bin_num+1):(i*bin_num)]
        } else {
            data_df[(bin_num*2*i-bin_num+1):(bin_num*2*i),]$Hits <- rev(bin_control$score[((i-1)*bin_num+1):(i*bin_num)])
        }
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
print(paste0("writing TSS_10kb......"))
write.csv(score, file = file.path(path, "TSS_10kb.csv"))
print("===========Finish!===========")
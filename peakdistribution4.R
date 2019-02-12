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
BinChipseq <- function(reads, bins){
  mcols(bins)$score = countOverlaps(bins, reads)
  return(bins)
}
TileCollection <- function(collection, bin_num) {
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
                                  strand = rep(collection$strand[i], bin_num),
                                  name = rep(collection$external_gene_name[i], bin_num))
        } else {
            data_df_temp <- data.frame(seqname = rep(collection$chromosome_name[i], bin_num),
                                       start = start_list, 
                                       end = end_list,
                                       strand = rep(collection$strand[i], bin_num),
                                       name = rep(collection$external_gene_name[i], bin_num))
            data_df <- rbind(data_df, data_df_temp)
        }
    }
    GRanges_temp <- GRanges(seqnames = Rle(data_df$seqname),
                            ranges = IRanges(start = data_df$start, end = data_df$end),
                            strand = Rle(data_df$strand),
                            name = Rle(data_df$name))
    return(GRanges_temp)
}

# parameter set
path <- "/mnt/c/chipseq/test7/"
cpg_file <-"cpg_mm10.bed"
chrom_list <- c(1:19, "X", "Y")
up_position <- 5000
down_position <- 5000
bin_num <- 50
cpg_n <- 1
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "mmusculus_gene_ensembl",  
                host = "jul2018.archive.ensembl.org")  
ds <- useDataset("mmusculus_gene_ensembl", mart = mart)

# load file
cpg_mm10 <- import.bed(file.path(path, cpg_file))

# main
print("+++++++++loop start+++++++++")
for (num in 1:length(chrom_list)){
    chrom <- chrom_list[num]
    chrom_name <- paste0("chr", chrom)
    cpg <- cpg_mm10[cpg_mm10@seqnames == chrom_name]
    print(paste0("analysis ", chrom_name, "......"))

    # import gene info
    gene_info <- getBM(attributes = c("ensembl_gene_id","external_gene_name",
                                      "chromosome_name", "start_position",
                                      "end_position", "strand", "gene_biotype"), 
                       mart = ds, filter = "chromosome_name", values = chrom)
    gene_info <- distinct(gene_info, external_gene_name, .keep_all = TRUE)
    gene_info$TSS <- rep(0, nrow(gene_info))
    for (i in 1:nrow(gene_info)) {
        if (gene_info$strand[i] == 1) {
            gene_info$TSS[i] = gene_info$start_position[i]
        } else {
            gene_info$TSS[i] = gene_info$end_position[i]
        }
    }
    gene_info_G <- GRanges(seqnames = Rle(chrom_name),
                           ranges = IRanges(start = gene_info$TSS-up_position, end = gene_info$end_position+down_position),
                           strand = Rle("*"),
                           name = gene_info$external_gene_name)
    gene_info_G <- BinChipseq(cpg, gene_info_G)
    gene_info_fil_G <- gene_info_G[gene_info_G$score >= cpg_n]
    gene_info_fil <- data.frame(chromosome_name = gene_info_fil_G@seqnames,
                                start_position = start(gene_info_fil_G),
                                end_position = end(gene_info_fil_G),
                                strand = gene_info_fil_G@strand,
                                external_gene_name = gene_info_fil_G$name,
                                stringsAsFactors = FALSE)
    gene_info_fil_bins <- TileCollection(gene_info_fil, bin_num)
    
    # overlap
    bin_cpg <- BinChipseq(cpg, gene_info_fil_bins)
    data_df <- data.frame(name = rep(NA, nrow(gene_info_fil)*bin_num),
                          chrom = rep(NA, nrow(gene_info_fil)*bin_num),
                          position = rep(NA, nrow(gene_info_fil)*bin_num),
                          Hits = rep(NA, nrow(gene_info_fil)*bin_num))
    for (i in 1:nrow(gene_info_fil)) {
        print(gene_info_fil$external_gene_name[i])
        data_df[(bin_num*(i-1)+1):(bin_num*i),]$name <- rep(gene_info_fil$external_gene_name[i], bin_num)
        data_df[(bin_num*(i-1)+1):(bin_num*i),]$chrom <- rep(chrom_name, bin_num)
        data_df[(bin_num*(i-1)+1):(bin_num*i),]$position <- c(1:bin_num)
        if (as.logical(bin_cpg@strand[(i-1)*bin_num+1] == "+")) {
            data_df[(bin_num*(i-1)+1):(bin_num*i),]$Hits <- bin_cpg$score[((i-1)*bin_num+1):(i*bin_num)]
        } else {
            data_df[(bin_num*(i-1)+1):(bin_num*i),]$Hits <- rev(bin_cpg$score[((i-1)*bin_num+1):(i*bin_num)])
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
print(paste0("writing cpg_on_TSS_5kb......"))
write.csv(score, file = file.path(path, "cpg_on_TSS_10kb.csv"))
print("===========Finish!===========")
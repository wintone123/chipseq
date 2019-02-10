# load library
library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(chipseq)
library(biomaRt)
library(BSgenome.Mmusculus.UCSC.mm10)  
library(Gviz)
library(tidyverse)

# load function
prepareChipseq <- function(reads){
  frag_len = median(estimate.mean.fraglen(reads))
  cat(paste0('median size is ', round(frag_len)))
  reads_extended = resize(reads, width = frag_len)
  return(trim(reads_extended))
}
BinChipseq <- function(reads, bins){
  mcols(bins)$score = countOverlaps(bins, reads)
  return(bins)
}
import_np <- function(path){
  temp <- read.delim2(path, quote = "/", head = FALSE)
  GRanges_temp <- GRanges(seqnames = Rle(temp[,1]),
                          ranges = IRanges(start = temp[,2], end = temp[,3]),
                          strand = Rle("*"),
                          name = as.character(temp[,4]),
                          score = as.numeric(temp[,5]))
  return(GRanges_temp)
}
overlap_peaks <- function(GRanges_1, GRanges_2){
    df_1 <- data.frame(start(GRanges_1), end(GRanges_1))
    df_2 <- data.frame(start(GRanges_2), end(GRanges_2))
    start_out <- vector()
    end_out <- vector()
    loop_s <- 1
    for(i in c(1:nrow(df_1))){
        start_1 <- df_1[i,1]
        end_1 <- df_1[i,2]
        for(j in c(loop_s:nrow(df_2))){
            start_2 <- df_2[j,1]
            end_2 <- df_2[j,2]
            if (end_2 <= start_1){
                loop_s <- j
            } else if (start_2 >= end_1) {
                break
            } else{
                if(start_2 <= start_1){
                    start_out <- append(start_out, start_1)
                    if(end_2 <= end_1 & end_2 > start_1){
                        end_out <- append(end_out, end_2)
                    } else if(end_2 > end_1){
                        end_out <- append(end_out, end_1)
                    }
                } else if(start_2 > start_1 & start_2 < end_1){
                    start_out <- append(start_out, start_2)
                    if(end_2 <= end_1){
                        end_out <- append(end_out, end_2)
                    } else{
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

# load peak file
rep1 <- import.bed("c:/chipseq/test4/rep1_best_chr6.bed")
rep2 <- import.bed("c:/chipseq/test4/rep2_best_chr6.bed")
control <- import.bed("c:/chipseq/test4/control_best_chr6.bed")
rep1_peak <- import_np("c:/chipseq/test4/rep1_best_peaks_chr6.narrowpeak")
rep2_peak <- import_np("c:/chipseq/test4/rep2_best_peaks_chr6.narrowpeak")
rep1_peak_1 <- import_np("c:/chipseq/test4/rep1_best_peaks_chr6_1.narrowpeak")
rep2_peak_1 <- import_np("c:/chipseq/test4/rep2_best_peaks_chr6_1.narrowpeak")
rep1_ext_peak <- import_np("c:/chipseq/test3/rep1_best_chr6_ext_peaks.narrowPeak")
rep2_ext_peak <- import_np("c:/chipseq/test3/rep1_best_chr6_ext_peaks.narrowPeak")


# pic
rep1_peak_track <- AnnotationTrack(rep1_peak, chromosome = "chr6", genome = "mm10", name = "rep1", col = "blue", fill = "blue", stacking = "dense")
rep1_ext_peak_track <- AnnotationTrack(rep1_ext_peak, chromosome = "chr6", genome = "mm10", name = "rep1_ext", col = "blue", fill = "blue", stacking = "dense")
rep1_track <- DataTrack(rep1_on_chr6, type = "histogram", ylim = c(0,150), genome = "mm10", size = 2, name = "rep1")
control_track <- DataTrack(control_on_chr6, type = "histogram", ylim = c(0,150), genome = "mm10", size = 2, name = "input")
enriched_track <- AnnotationTrack(enriched, chromosome = "chr6", genome = "mm10", name = "rep1", col = "yellow", fill = "yellow", stacking = "dense")
axis_track <- GenomeAxisTrack()
genome_track <- IdeogramTrack(chromosome = "chr6", genome = "mm10")
plotTracks(c(genome_track, axis_track, rep1_track, rep1_peak_track, rep1_peak_1_track, rep2_track, rep2_peak_track, rep2_peak_1_track),
           from = 122707000, to = 122715000)

# load element
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "mmusculus_gene_ensembl",  # human (hsapiens_gene_ensembl)
                host = "jul2018.archive.ensembl.org")  # ensembl release 93 (july 2018)
ds <- useDataset("mmusculus_gene_ensembl", mart = mart)
nanog_info_df <- getBM(attributes = c("ensembl_gene_id","external_gene_name", "chromosome_name", "strand",
                                      "start_position", "end_position", "exon_chrom_start", "exon_chrom_end",
                                      "5_utr_start", "5_utr_end", "3_utr_start", "3_utr_end"), 
                       mart = ds, filter = "external_gene_name", values = "nanog")
egs <- getBM(attributes = c("ensembl_gene_id","external_gene_name",
                            "chromosome_name", "start_position",
                            "end_position", "strand", "gene_biotype"), mart = ds, 
                            filter = "chromosome_name", 
                            values = "6")

# peaks in tss +-1000bp
egs$TSS <- rep(0, nrow(egs))
for (i in 1:nrow(egs)) {
    if (egs$strand[i] == 1) {
        egs$TSS[i] = egs$start_position[i]
    } else {
        egs$TSS[i] = egs$end_position[i]
    }
}
tss_chr6 <- GRanges(seqnames = Rle("chr6"),
                    ranges = IRanges(start = egs$TSS-200, end = egs$TSS+200),
                    strand = Rle("*"))
overlap_peaks <- overlap_peaks(rep1_peak_1, rep2_peak_1)
egs_fil <- egs[unique(queryHits(findOverlaps(tss_chr6, overlap_peaks))),]
for (i in 1:nrow(egs_fil)) {
  if (i == 1) {
    if (egs_fil$strand[i] == 1) {
      start_temp <- egs_fil$start_position[i] + seq(-1000, 900, length.out = 20)
    } else {
      start_temp <- egs_fil$end_position[i] + seq(-900, 1000, length.out = 20)
    }
    strand_temp <- rep(egs_fil$strand[i], 20)
    tss <- data.frame(start = start_temp, strand = strand_temp)
  } else {
    if(egs_fil$strand[i] == 1){
      start_temp <- egs_fil$start_position[i] + seq(-1000, 900, length.out = 20)
    } else {
      start_temp <- egs_fil$end_position[i] + seq(-900, 1000, length.out = 20)
    }
    strand_temp <- rep(egs_fil$strand[i], 20)
    tss_temp <- data.frame(start = start_temp, strand = strand_temp)
    tss <- rbind(tss, tss_temp)
  }
}
tss_bins <- GRanges(seqnames = Rle("chr6"),
                    ranges = IRanges(start = tss$start, width = 100),
                    strand = Rle(tss$strand))
rep1_on_tss <- BinChipseq(rep1, tss_bins)
rep2_on_tss <- BinChipseq(rep2, tss_bins)
control_on_tss <- BinChipseq(control, tss_bins)
total_on_tss <- tss_bins
total_on_tss$score <- rep1_on_tss$score + rep2_on_tss$score
score <- data.frame(position = c(1:20), score = rep(0, 20), control = rep(0, 20))
for (i in 1:nrow(egs_fil)) {
  if (as.logical(total_on_tss@strand[(i-1)*20+1] == "+")) {
    score_temp <- total_on_tss$score[((i-1)*20+1):(i*20)]
  } else {
    score_temp <- rev(total_on_tss$score[((i-1)*20+1):(i*20)])
  }
  score$score <- score$score + score_temp
} 
for (i in 1:nrow(egs_fil)) {
  if (as.logical(control_on_tss@strand[(i-1)*20+1] == "+")) {
    score_temp <- control_on_tss$score[((i-1)*20+1):(i*20)]
  } else {
    score_temp <- rev(control_on_tss$score[((i-1)*20+1):(i*20)])
  }
  score$control <- score$control + score_temp
} 

# peak in -100% to 200% area of gene
gene_chr6 <- GRanges(seqnames = Rle("chr6"),
                     ranges = IRanges(start = egs$start_position, end = egs$end_position),
                     strand = Rle("*"))
egs_fil <- egs[unique(queryHits(findOverlaps(gene_chr6, overlap_peaks))),]
for(i in c(1:nrow(egs_fil))){
  start <- 2 * egs_fil$start_position[i] - egs_fil$end_position[i]
  end <- 2 * egs_fil$end_position[i] - egs_fil$start_position[i]
  seq_list <- seq(start, end, length.out = 271)
  start_list <- vector()
  end_list <- vector()
  for (j in 1:length(seq_list)) {
      if (j == 1) {
        start_list <- append(start_list, seq_list[j])
      } else if(j == length(seq_list)) {
        end_list <- append(end_list, seq_list[j])
      } else {
        start_list <- append(start_list, seq_list[j]+1)
        end_list <- append(end_list, seq_list[j])
      }
    }
  if (i ==1 ) {
	  egs_ext <- data.frame(start = start_list, end = end_list, strand = rep(egs_fil$strand[i], 270))
  } else {
	  egs_ext_temp <- data.frame(start = start_list, end = end_list, strand = rep(egs_fil$strand[i], 270))
	  egs_ext <- rbind(egs_ext, egs_ext_temp)
  }
}
egs_ext_G <- GRanges(seqnames = Rle("chr6"),
					           ranges = IRanges(start = egs_ext$start, end = egs_ext$end),
					           strand = Rle(egs_ext$strand))
rep1_on_egs_ext <- BinChipseq(rep1, egs_ext_G)
control_on_egs_ext <- BinChipseq(control, egs_ext_G)
score <- data.frame(position = c(1:270), score = rep(0, 270), control = rep(0, 270))
for (i in c(1:nrow(egs_fil))) {
  if (as.logical(rep1_on_egs_ext@strand[(i-1)*270+1] == "+")) {
    score_temp <- rep1_on_egs_ext$score[((i-1)*270+1):(i*270)]
  } else {
    score_temp <- rev(rep1_on_egs_ext$score[((i-1)*270+1):(i*270)])
  }
  score$score <- score$score + score_temp
} 
for (i in 1:nrow(egs_fil)) {
  if (as.logical(control_on_egs_ext@strand[(i-1)*270+1] == "+")) {
    score_temp <- control_on_egs_ext$score[((i-1)*270+1):(i*270)]
  } else {
    score_temp <- rev(control_on_egs_ext$score[((i-1)*270+1):(i*270)])
  }
  score$control <- score$control + score_temp
} 
score <- gather(score, c(score, control), key = "item", value = "hits")

# ratio calculation
ratio_list <- vector(length = length(ratio_bins))
control_score <- control_bins$score
rep1_score <- rep1_bins$score
for (i in 1:length(ratio_bins)) {
  if (control_score[i] !=0 & rep1_score[i] != 0) {
    ratio_list[i] <- rep1_score[i] / control_score[i]
  } else {
    ratio_list[i] <- 0
  }
}
ratio_bins$score <- ratio_list

# -10kb start end +10kb
gene_chr6 <- GRanges(seqnames = Rle("chr6"),
                     ranges = IRanges(start = egs$start_position, end = egs$end_position),
                     strand = Rle("*"))
egs_fil <- egs[unique(queryHits(findOverlaps(gene_chr6, overlap_peaks))),]
for(i in c(1:nrow(egs_fil))){
  start <- egs_fil$start_position[i]
  end <- egs_fil$end_position[i]
  before_start <- start - 10000
  after_end <- end + 10000
  seq_list <- vector()
  seq_list <- append(seq_list, seq(before_start, start, length.out = 51))
  seq_list <- append(seq_list, seq(start, end, length.out = 101))
  seq_list <- append(seq_list, seq(end, after_end, length.out = 51))
  seq_list <- unique(seq_list)
  start_list <- vector()
  end_list <- vector()
  for (j in 1:length(seq_list)) {
      if (j == 1) {
        start_list <- append(start_list, seq_list[j])
      } else if(j == length(seq_list)) {
        end_list <- append(end_list, seq_list[j])
      } else {
        start_list <- append(start_list, seq_list[j]+1)
        end_list <- append(end_list, seq_list[j])
      }
    }
  if (i ==1 ) {
	  egs_ext <- data.frame(start = start_list, end = end_list, strand = rep(egs_fil$strand[i], 270))
  } else {
	  egs_ext_temp <- data.frame(start = start_list, end = end_list, strand = rep(egs_fil$strand[i], 270))
	  egs_ext <- rbind(egs_ext, egs_ext_temp)
  }
}
egs_ext_G <- GRanges(seqnames = Rle("chr6"),
					           ranges = IRanges(start = egs_ext$start, end = egs_ext$end),
					           strand = Rle(egs_ext$strand))
rep1_on_egs_ext <- BinChipseq(rep1, egs_ext_G)
control_on_egs_ext <- BinChipseq(control, egs_ext_G)
score <- data.frame(position = c(1:270), score = rep(0, 270), control = rep(0, 270))
for (i in c(1:nrow(egs_fil))) {
  if (as.logical(rep1_on_egs_ext@strand[(i-1)*270+1] == "+")) {
    score_temp <- rep1_on_egs_ext$score[((i-1)*270+1):(i*270)]
  } else {
    score_temp <- rev(rep1_on_egs_ext$score[((i-1)*270+1):(i*270)])
  }
  score$score <- score$score + score_temp
} 
for (i in 1:nrow(egs_fil)) {
  if (as.logical(control_on_egs_ext@strand[(i-1)*270+1] == "+")) {
    score_temp <- control_on_egs_ext$score[((i-1)*270+1):(i*270)]
  } else {
    score_temp <- rev(control_on_egs_ext$score[((i-1)*270+1):(i*270)])
  }
  score$control <- score$control + score_temp
} 
score <- gather(score, c(score, control), key = "item", value = "hits")

# heatmap
tss_fil <- filter(TSS_10kb, item == "test") %>% spread(key = position, value = Hits) 
tss_fil2 <- select(tss_fil, c(4:33)) 
tss_mat <- as.matrix(tss_fil2)
rownames(tss_mat) <- tss_fil$name
pheatmap(tss_mat, cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = FALSE, show_colnames = FALSE, width = 3)
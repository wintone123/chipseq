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

# load peak file
rep1 <- import.bed("c:/chipseq/test4/rep1_best_chr6.bed")
rep2 <- import.bed("c:/chipseq/test4/rep2_best_chr6.bed")
control <- import.bed("c:/chipseq/test3/control_best_chr6.bed")
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
                            "end_position", "strand"), mart = ds, 
                            filter = "chromosome_name", 
                            values = "6")

# peaks in tss +-1000bp
egs$TSS <- rep(0, nrow(egs))
for(i in c(1:nrow(egs))){
    if(egs$strand[i] == 1){
        egs$TSS[i] = egs$start_position[i]
    } else{
        egs$TSS[i] = egs$end_position[i]
    }
}
tss_chr6 <- GRanges(seqnames = Rle("chr6"),
                    ranges = IRanges(start = egs$TSS-200, end = egs$TSS+200),
                    strand = Rle("*"))
egs_fil <- egs[unique(queryHits(findOverlaps(tss_chr6, rep1_peak_1))),]
for(i in c(1:nrow(egs_fil))){
  if(i == 1){
    if(egs_fil$strand[i] == 1){
      start_temp <- egs_fil$start_position[i] + seq(-1000, 900, length.out = 20)
    }else {
      start_temp <- egs_fil$end_position[i] + seq(-900, 1000, length.out = 20)
    }
    strand_temp <- rep(egs_fil$strand[i], 20)
    tss <- data.frame(start = start_temp, strand = strand_temp)
  }else{
    if(egs_fil$strand[i] == 1){
      start_temp <- egs_fil$start_position[i] + seq(-1000, 900, length.out = 20)
    }else {
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
score <- data.frame(position = c(1:20), score = rep(0, 20))
for(i in c(1:nrow(egs_fil))){
  if(as.logical(rep1_on_tss@strand[(i-1)*20+1] == "+")){
    score_temp <- rep1_on_tss$score[c(((i-1)*20+1):(i*20))]
  } else{
    score_temp <- rev(rep1_on_tss$score[c(((i-1)*20+1):(i*20))])
  }
  score$score <- score$score + score_temp
}
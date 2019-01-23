# load library
library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(chipseq)
library(biomaRt)
library(BSgenome.Mmusculus.UCSC.mm10)  
library(Gviz)

# load function
prepareChipseq <- function(reads){
  frag_len = median(estimate.mean.fraglen(reads))
  cat(paste0('median size is ', round(frag_len)))
  reads_extended = resize(reads, width = frag_len)
  return(trim(reads_extended))
}
BinChipseq <- function(reads, bins, name){
  mcols(bins)$score = countOverlaps(bins, reads)
  return(bins)
}

# load peak file
rep1 <- import.bed("c:/chipseq/test3/rep1_best_chr6.bed")
rep2 <- import.bed("c:/chipseq/test3/rep2_best_chr6.bed")
control <- import.bed("c:/chipseq/test3/control_best_chr6.bed")
rep1_peak <- import_np("c:/chipseq/test3/rep1_best_peaks_chr6.narrowpeak")
rep2_peak <- import_np("c:/chipseq/test3/rep2_best_peaks_chr6.narrowpeak")
rep1_ext_peak <- import_np("c:/chipseq/test3/rep1_best_chr6_ext_peaks.narrowPeak")
rep2_ext_peak <- import_np("c:/chipseq/test3/rep1_best_chr6_ext_peaks.narrowPeak")


# pic
rep1_peak_track <- AnnotationTrack(rep1_peak, chromosome = "chr6", genome = "mm10", name = "rep1", col = "blue", fill = "blue", stacking = "dense")
rep1_ext_peak_track <- AnnotationTrack(rep1_ext_peak, chromosome = "chr6", genome = "mm10", name = "rep1_ext", col = "blue", fill = "blue", stacking = "dense")
rep1_track <- DataTrack
enriched_track <- AnnotationTrack(enriched, chromosome = "chr6", genome = "mm10", name = "rep1", col = "yellow", fill = "yellow", stacking = "dense")
axis_track <- GenomeAxisTrack()
genome_track <- IdeogramTrack(chromosome = "chr6", genome = "mm10")
plotTracks(c(genome_track, axis_track, rep1_peak_track, rep1_ext_peak_track),
           from = 122600000, to = 122700000)
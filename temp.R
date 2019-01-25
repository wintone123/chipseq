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
rep1_track <- DataTrack(rep1_on_chr6, type = "histogram", ylim = c(0,150), genome = "mm10", size = 2, name = "rep1")
control_track <- DataTrack(control_on_chr6, type = "histogram", ylim = c(0,150), genome = "mm10", size = 2, name = "input")
enriched_track <- AnnotationTrack(enriched, chromosome = "chr6", genome = "mm10", name = "rep1", col = "yellow", fill = "yellow", stacking = "dense")
axis_track <- GenomeAxisTrack()
genome_track <- IdeogramTrack(chromosome = "chr6", genome = "mm10")
plotTracks(c(genome_track, axis_track, control_track, rep1_track, rep1_peak_track, rep1_ext_peak_track),
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

# load library
library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(chipseq)
library(biomaRt)
library(BSgenome.Mmusculus.UCSC.mm9)  # human (BSgenome.Hsapiens.UCSC.hg38)
library(Gviz)

# si for mm9 
genome <- BSgenome.Mmusculus.UCSC.mm9  # human (BSgenome.Hsapiens.UCSC.hg38)
si <- seqinfo(genome)
si <- si[paste0('chr', c(1:19, 'X', 'Y'))]

# load interestined chromosome and position (mouse mm9 genome)
chromosome <- "chr6"  # chromosome (chr 1~19 or X, Y)
start <- 122600000  # start point
end <- 122900000  # end point

# load bed file and read extension
prepareChipseq <- function(reads){
  frag_len = median(estimate.mean.fraglen(reads))
  cat(paste0('median size is ', round(frag_len)))
  reads_extended = resize(reads, width = frag_len)
  return(trim(reads_extended))
}
path = "d:/chipseq/test1"  # bed file path 
load_list <- list.files(path)
read_list <- vector()
for(file_name in load_list){
  if(substr(file_name, start = nchar(file_name)-2, stop = nchar(file_name)) == "bed"){
    if(substr(file_name, start = nchar(file_name)-8, stop = nchar(file_name)-4) != "peaks"){
      base_name <- strsplit(file_name, split = ".", fixed = TRUE)[[1]][1]
      assign(base_name, prepareChipseq(import.bed(file.path(path, file_name))))
      read_list <- c(read_list, base_name)
    }
  }
}

# bm for mm9 
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "mmusculus_gene_ensembl",  # human (hsapiens_gene_ensembl)
                host = "jul2018.archive.ensembl.org")  # ensembl release 93 (july 2018)
fm <- Gviz:::.getBMFeatureMap()
fm["symbol"] <- "external_gene_name"
bm <- BiomartGeneRegionTrack(chromosome = chromosome, genome = "mm9",  # human(hg38)
                             start = start, end = end,
                             biomart = mart, filter = list("with_refseq_mrna" = TRUE),
                             size = 4, name = "Refseq", utr5 = "red3", utr3 = "red3",
                             protein_coding = "black", col.line = NULL, cex = 7,
                             collapseTranscript = "longest", featureMap = fm)

# export bm plot
AT <- GenomeAxisTrack()
plotTracks(c(bm,AT), from = start, to = end,
           transcriptAnnotation = "symbol", window = "auto",
           cex.title = 1, fontsize = 10)

# bin generation 
binsize <- 200  # binsize (200 bp)
bins <- tileGenome(si[chromosome], tilewidth = binsize,
                   cut.last.tile.in.chrom = TRUE)

# normalization factor calculation
length_list <- vector()
for(base_name in read_list){
  length_list <- c(length_list, length(eval(parse(text = base_name))))
}
fold_list <- length_list / min(length_list)

# binning with normalization
BinChipseq <- function(reads, bins, name){
  index <- match(name, read_list)
  mcols(bins)$score = countOverlaps(bins, reads) / fold_list[index]
  return(bins)
}
bin_track_list <- vector()
for(base_name in read_list){
  bin_name <- paste0(base_name, "_200bins")
  assign(bin_name, BinChipseq(eval(parse(text = base_name)), bins, base_name))
  bin_track_name <- paste0(bin_name, "_track")
  if(bin_track_name == "control_track"){
    assign(bin_track_name, DataTrack(eval(parse(text = bin_name)), strand = "*", genome = "mm9",  # human (hg38) 
    col.histogram = "gray", fill.histrogram = "black",
    name = base_name, col.axis = "black", cex.axis = 0.4, ylim = c(0,150))) 
  }
  else{
    assign(bin_track_name, DataTrack(eval(parse(text = bin_name)), strand = "*", genome = "mm9",  # human (hg38)
    col.histogram = "gray", fill.histrogram = "steelblue",
    name = base_name, col.axis = "steelblue", cex.axis = 0.4, ylim = c(0,150))) 
  }
  bin_track_list <- c(bin_track_list, bin_track_name)
}

# export bins to bm plot
bin_track_temp <- vector()
for(bin_track_name in bin_track_list){
  bin_track_temp <- c(bin_track_temp, eval(parse(text = bin_track_name)))
}
plotTracks(c(bin_track_temp, bm, AT), from = start, 
           to = end, transcriptAnnotation = "symbol", window = "auto", 
           type = "histogram", cex.title = 0.7, fontsize = 10)

# load peaks
peak_list <- vector()
peak_track_list <- vector()
for(file_name in load_list){
  if(substr(file_name, start = nchar(file_name)-2, stop = nchar(file_name)) == "bed"){
    if(substr(file_name, start = nchar(file_name)-8, stop = nchar(file_name)-4) == "peaks"){
      peak_name <- strsplit(file_name, split = ".", fixed = TRUE)[[1]][1]
      assign(peak_name, import.bed(file.path(path, file_name)))
      peak_list <- c(peak_list, peak_name)
      peak_track_name <- paste0(peak_name, "_track")
      assign(peak_track_name, AnnotationTrack(eval(parse(text = peak_name)), genome = "mm9",  # human (hg38)
             name = peak_name, chromosome = chromosome, shape = "box", fill = "blue3", size = 2))
      peak_track_list <- c(peak_track_list, peak_track_name)
    }  
  }
}


# export peaks to bm plot
bin_peak_track_temp <- vector()
for(track_name in sort(c(bin_track_list, peak_track_list))){
  bin_peak_track_temp <- c(bin_peak_track_temp, eval(parse(text = track_name)))
}
plotTracks(c(bin_peak_track_temp, bm, AT), from = start, to = end, 
           transcriptAnnotation = "symbol", window = "auto",
           type = "histogram", cex.title = 0.7, fontsize = 10)

# peaks enrichment
peak_temp <- vector()
for(peak_name in peak_list){
  peak_temp <- c(peak_temp, eval(parse(text = peak_name)))
}
enriched_regions = Reduce(subsetByOverlaps, peak_temp)
enriched_regions_track <- AnnotationTrack(enriched_regions, genome = "mm9", name = "enriched regions",  # human (hg38)
                                          chromosome = chromosome, shape = "box", fill = "green3", size = 2)

# export peaks enrichment to bm plot
bin_peak_track_temp <- vector()
enrich_list <- sort(c(bin_track_list, peak_track_list))
for(name in enrich_list){
  bin_peak_track_temp <- c(bin_peak_track_temp, eval(parse(text = name)))
}
plotTracks(c(bin_peak_track_temp, enriched_regions_track, bm, AT), from = start, to = end, 
           transcriptAnnotation = "symbol", window = "auto",
           type = "histogram", cex.title = 0.7, fontsize = 10)

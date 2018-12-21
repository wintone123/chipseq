# load library
library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(biomaRt)
library(chipseq)
library(BSgenome.Mmusculus.UCSC.mm10)  # human (BSgenome.Hsapiens.UCSC.hg38)
library(Gviz)

# load necessary data
genome <- BSgenome.Mmusculus.UCSC.mm10
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "mmusculus_gene_ensembl",  # human (hsapiens_gene_ensembl)
                host = "jul2018.archive.ensembl.org")  # ensembl release 93 (july 2018)


# load file
prepareChipseq <- function(reads){
  frag_len = median(estimate.mean.fraglen(reads))
  reads_extended = resize(reads, width = frag_len)
  return(trim(reads_extended))
}
path <- "/mnt/c/chipseq/test2/split_files"
dir.create(file.path(path,"output"), showWarnings = FALSE)  # create output folder
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

# geonme_list
geonme_list <- c(1:19, "X","Y")

# analysis process
for(chromsome in geonme_list){
    chromsome <- paste0("chr",chromsome)

    # create data export file (1)
    colname <- vector()
    export <- promoter_peaks_export <- data.frame(t(unlist(read_list)), stringsAsFactors = FALSE)
    export_temp <- vector()
    for(name in read_list){
        export_temp <- append(export_temp, length(eval(parse(text = name))))
        print(paste(name, "has", length(eval(parse(text = name))), "peaks", sep = " "))
    }
        export[nrow(export)+1,] <- export_temp
        rownames(export)[nrow(export)] <- "total peaks"

    # egs isolation
    ds <- useDataset("mmusculus_gene_ensembl", mart = mart)  # human (hsapiens_gene_ensembl)
    egs <- getBM(attributes = c("ensembl_gene_id","external_gene_name",
                                "chromosome_name", "start_position",
                                "end_position", "strand"), mart = ds, 
                                filter = "chromosome_name", 
                                values = strsplit(chromosome, split = "r")[[1]][2])
    egs_regions <- GRanges(seqnames = Rle(paste0("chr", egs$chromosome_name)),
                        ranges = IRanges(start = egs$start_position, end = egs$end_position),
                        strand = Rle(rep("*", nrow(egs))),
                        gene = egs$external_gene_name)

    # peaks in egs regions
    export_temp <- vector()
    for(name in read_list){
        egs_peak_sum <- sum(countOverlaps(egs_regions, eval(parse(text = name))))
        print(paste(name, "has", egs_peak_sum, "peaks in egs regions", sep = " "))
        export_temp <- append(export_temp, egs_peak_sum)
    }
    export[nrow(export)+1,] <- export_temp
    rownames(export)[nrow(export)] <- "egs regions"

    # create data export file (2)
    colnames(export) <- export[1,]
    export <- export[c(2:nrow(export)),]
    write.csv(export, paste0(path, "/output/", chromosome, "_export.csv"))
}
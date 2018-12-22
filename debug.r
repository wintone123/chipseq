# load library
library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(biomaRt)
library(chipseq)
library(BSgenome.Mmusculus.UCSC.mm10)  
library(Gviz)

# load necessary data
genome <- BSgenome.Mmusculus.UCSC.mm10
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "mmusculus_gene_ensembl",  # human (hsapiens_gene_ensembl)
                host = "jul2018.archive.ensembl.org")  # ensembl release 93 (july 2018)
prepareChipseq <- function(reads){
  frag_len = median(estimate.mean.fraglen(reads))
  reads_extended = resize(reads, width = frag_len)
  return(trim(reads_extended))
}
path <- "/mnt/c/chipseq/test2/split_files"  # input file path
load_list <- list.files(path)
dir.create(file.path(path,"output"), showWarnings = FALSE)  # create output folder

# genome_list
genome_list <- c(1:19, "X","Y")
print("=============Let's Go!=============")
# analysis process
for(chrom_1 in genome_list){
  chrom <- paste0("chr", chrom_1)

  # load file
  read_list <- vector()
  for(file_name in load_list){
    if(substr(file_name, start = nchar(file_name)-2, stop = nchar(file_name)) == "bed"){
      if(substr(file_name, start = nchar(file_name)-7, stop = nchar(file_name)-4) == chrom){
        base_name <- strsplit(file_name, split = ".", fixed = TRUE)[[1]][1]
        assign(base_name, prepareChipseq(import.bed(file.path(path, file_name))))
        read_list <- c(read_list, base_name)
      }
    }
  }

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
  print("+++++++++++++++++++++++++++++++++++")

  # egs isolation
  ds <- useDataset("mmusculus_gene_ensembl", mart = mart)  
  egs <- getBM(attributes = c("ensembl_gene_id","external_gene_name",
                              "chromosome_name", "start_position",
                              "end_position", "strand"), mart = ds, 
                              filter = "chromosome_name", 
                              values = chrom_1)
  egs_regions <- GRanges(seqnames = Rle(chrom),
                         ranges = IRanges(start = egs$start_position, end = egs$end_position),
                         strand = Rle("*"),
                         gene = egs$external_gene_name)

  # egs overlapping regions
  egs_overlap <- data.frame("chromosome_name","start","end", stringsAsFactors=FALSE)
  egs <- egs[order(egs$start_position),]
  for(i in c(1: nrow(egs))){
    if(i == 1){
      start <- egs[i,]$start_position
      end <- egs[i,]$end_position
      egs_overlap_list <- c(start, end)
    } else{
      if(egs[i,]$start_position > end){
        start <- egs[i,]$start_position
        end <- egs[i,]$end_position
        egs_overlap_list <- append(egs_overlap_list, c(start, end))
      } else{
        if(egs[i,]$end_position >= end){
          end <- egs[i,]$end_position
          egs_overlap_list[length(egs_overlap_list)] <- end
        }
      }
    }
  }
  for(i in c(1:(length(egs_overlap_list)/2))){
    egs_overlap[nrow(egs_overlap)+1,] <- c(chrom, egs_overlap_list[i*2-1], egs_overlap_list[i*2])
  }
  colnames(egs_overlap) <- egs_overlap[1,]
  egs_overlap <- egs_overlap[c(2:nrow(egs_overlap)),]
  egs_overlap$start <- as.numeric(egs_overlap$start)
  egs_overlap$end <- as.numeric(egs_overlap$end)
  egs_overlap_regions <- GRanges(seqnames = Rle(chrom),
                                ranges = IRanges(start = egs_overlap$start, end = egs_overlap$end),
                                strand = Rle("*"))

  # intragenic regions isolation from egs overlap (gene distance > 200kb)
  egs_overlap <- egs_overlap[order(egs_overlap$start),]
  gd_list <- c(1)
  for(i in c(1:nrow(egs_overlap))){
    if(i == 1){
      start <- egs_overlap[i,]$star-100000
      end <- egs_overlap[i,]$end+100000
      gd_list <- append(gd_list, c(start, end))
    } else{
      if(egs_overlap[i,]$start - end > 200000){
        start <- egs_overlap[i,]$start-100000
        end <- egs_overlap[i,]$end+100000
        gd_list <- append(gd_list, c(start, end))
      } else{
        end <- egs_overlap[i,]$end+100000
        gd_list[length(gd_list)] <- end
      }
    }
  }
  if(chrom_1 == "X"){
    n <- 20
  } else if(chrom_1 == "Y"){
    n <- 21
  } else{
    n <- as.numeric(chrom_1)
  }
  gd_list <- append(gd_list, seqlengths(genome)[n])
  gd_temp <- data.frame("chromosome_name", "start", "end", stringsAsFactors=FALSE)
  if(gd_list[2] - gd_list[1] > 100000){
    a <- 1
  } else{
    a<- 2
  }
  if(gd_list[length(gd_list)] - gd_list[length(gd_list)-1] > 100000){
    b <- length(gd_list)/2
  } else{
    b <- length(gd_list)/2-1
  }
  for(i in c(a:b)){
    gd_temp[nrow(gd_temp)+1,] <- c(chrom, gd_list[i*2-1], gd_list[i*2])
  }
  colnames(gd_temp) <- gd_temp[1,]
  gd_temp <- gd_temp[c(2:nrow(gd_temp)),]
  gd_temp$start <- as.numeric(gd_temp$start)
  gd_temp$end <- as.numeric(gd_temp$end)
  gd_regions <- GRanges(seqnames = Rle(chrom), 
                        ranges = IRanges(start = gd_temp$start, end = gd_temp$end))

  # peaks in intragenic regions (gene distance > 200kb)
  export_temp <- vector()
  for(name in read_list){
    intragene_peak_sum <- sum(countOverlaps(gd_regions, eval(parse(text = name))))
    print(paste(name, "has", intragene_peak_sum, "peaks in intragenic regions", sep = " "))
    export_temp <- append(export_temp, intragene_peak_sum)
  }
  export[nrow(export)+1,] <- export_temp
  rownames(export)[nrow(export)] <- "intragenic regions"
  print("+++++++++++++++++++++++++++++++++++")

  # create data export file (2)
  colnames(export) <- export[1,]
  export <- export[c(2:nrow(export)),]
  write.csv(export, paste0(path, "/output/", chrom, "_export.csv"))
  print(paste0(path, "/output/", chrom, "_export.csv"))

  # remove files
  for(name in read_list){
    remove(name)
  }
  print("===============Done!===============")
}
# load library
library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(biomaRt)
library(chipseq)
library(BSgenome.Mmusculus.UCSC.mm10)  # human (BSgenome.Hsapiens.UCSC.hg38)
library(Gviz)

# si for mm9 
genome <- BSgenome.Mmusculus.UCSC.mm10  # human (BSgenome.Hsapiens.UCSC.hg38)

# load interestined chromosome and position (mouse mm9 genome)
chromosome <- "chr6"  # chromosome (chr 1~19 or X, Y)

# load bed file and read extension
prepareChipseq <- function(reads){
  frag_len = median(estimate.mean.fraglen(reads))
  reads_extended = resize(reads, width = frag_len)
  return(trim(reads_extended))
}
path = "/mnt/c/chipseq/test1"  # bed file path 
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

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "mmusculus_gene_ensembl",  # human (hsapiens_gene_ensembl)
                host = "jul2018.archive.ensembl.org")  # ensembl release 93 (july 2018)

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

# egs overlapping regions
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
egs_overlap <- data.frame("chromosome_name","start","end", stringsAsFactors=FALSE)
for(i in c(1:(length(egs_overlap_list)/2))){
  egs_overlap[nrow(egs_overlap)+1,] <- c(chromosome, egs_overlap_list[i*2-1], egs_overlap_list[i*2])
}
colnames(egs_overlap) <- egs_overlap[1,]
egs_overlap <- egs_overlap[c(2:nrow(egs_overlap)),]
egs_overlap$start <- as.numeric(egs_overlap$start)
egs_overlap$end <- as.numeric(egs_overlap$end)
egs_overlap_regions <- GRanges(seqnames = Rle(egs_overlap$chromosome_name),
                               ranges = IRanges(start = egs_overlap$start, end = egs_overlap$end),
                               strand = Rle(rep("*", nrow(egs_overlap))))

# peaks in egs overlapp regions
export_temp <- vector()
for(name in read_list){
  egs_overlap_peak_sum <- sum(countOverlaps(egs_overlap_regions, eval(parse(text = name))))
  print(paste(name, "has", egs_overlap_peak_sum, "peaks in egs overlap regions", sep = " "))
  export_temp <- append(export_temp, egs_overlap_peak_sum)
}
export[nrow(export)+1,] <- export_temp
rownames(export)[nrow(export)] <- "egs overlap regions"

# set promoter (+/- 200 bp around TSS)
egs$TSS <- ifelse(egs$strand == "1", egs$start_position, egs$end_position)
promoter_regions <- GRanges(seqnames = Rle(paste0("chr", egs$chromosome_name)),
                            ranges = IRanges(start = egs$TSS - 200, end = egs$TSS + 200),
                            strand = Rle(rep("*", nrow(egs))),
                            gene = egs$external_gene_name)

# peaks in promoter regions (+/- 200 bp around TSS)
export_temp <- vector()
for(name in read_list){
  promoter_peak_sum <- sum(countOverlaps(promoter_regions, eval(parse(text = name))))
  print(paste(name, "has", promoter_peak_sum, "peaks in promoter regions", sep = " "))
  export_temp <- append(export_temp, promoter_peak_sum)
}
export[nrow(export)+1,] <- export_temp
rownames(export)[nrow(export)] <- "promoter regions"

# exon isolation
exon <- getBM(attributes = c("ensembl_gene_id","external_gene_name",
                            "chromosome_name", "exon_chrom_start",
                            "exon_chrom_end", "strand"), mart = ds,
                            filter = "chromosome_name", 
                            values = strsplit(chromosome, split = "r")[[1]][2])
exon_regions <- GRanges(seqnames = Rle(paste0("chr", exon$chromosome_name)),
                        ranges = IRanges(start = exon$exon_chrom_start, end = exon$exon_chrom_end),
                        srtand = Rle(rep("*", nrow(exon))),
                        gene = exon$external_gene_name)

# peaks in exon regions
export_temp <- vector()
for(name in read_list){
  exon_peak_sum <- sum(countOverlaps(exon_regions, eval(parse(text = name))))
  print(paste(name, "has", exon_peak_sum, " peaks in exon regions", sep = " "))
  export_temp <- append(export_temp, exon_peak_sum)
}
export[nrow(export)+1,] <- export_temp
rownames(export)[nrow(export)] <- "exon regions"

# exon extend regions and intron from exon extend regions
exon_extend <- data.frame("external_gene_name", "exon_start", "exon_end", stringsAsFactors = FALSE)
intron_from_exon_extend <- data.frame("external_gene_name", "intron_start", "intron_end", stringsAsFactors = FALSE)
for(name in unique(exon$external_gene_name)){
  name_index_list <- which(exon$external_gene_name == name)
  name_temp <- exon[name_index_list,]
  name_temp <- name_temp[order(name_temp$exon_chrom_start),]
  for(i in c(1:nrow(name_temp))){
    if(i == 1){
      start <- name_temp[i,4]
      end <- name_temp[i,5]
      exon_extend_list <- c(start, end)
    } else{
      if(name_temp[i,4] > end){
        start <- name_temp[i,4]
        end <- name_temp[i,5]
        exon_extend_list <- append(exon_extend_list, c(start, end))
      } else{
        if(name_temp[i,5] >= end){
          end <- name_temp[i,5]
          exon_extend_list[length(exon_extend_list)] <- end
        }
      }
    }
  }
  for(i in c(1:(length(exon_extend_list)/2))){
    exon_extend[nrow(exon_extend)+1,] <- c(name, exon_extend_list[i*2-1], exon_extend_list[i*2])
  }
  if(length(exon_extend_list) >= 4){
    for(i in c(1:(length(exon_extend_list)/2-1))){
      intron_from_exon_extend[nrow(intron_from_exon_extend)+1,] <- c(name, exon_extend_list[i*2]+1, exon_extend_list[i*2+1]-1)
    }
  }
}
colnames(exon_extend) <- exon_extend[1,]
exon_extend <- exon_extend[c(2:nrow(exon_extend)),]
exon_extend$exon_start <- as.numeric(exon_extend$exon_start)
exon_extend$exon_end <- as.numeric(exon_extend$exon_end)
exon_extend_regions <- GRanges(seqnames = Rle(chromosome),
                                ranges = IRanges(start = exon_extend$exon_start, end = exon_extend$exon_end),
                                strand = Rle(rep("*", nrow(intron_from_exon_extend))))
colnames(intron_from_exon_extend) <- intron_from_exon_extend[1,]
intron_from_exon_extend <- intron_from_exon_extend[c(2:nrow(intron_from_exon_extend)),]
intron_from_exon_extend$intron_start <- as.numeric(intron_from_exon_extend$intron_start)
intron_from_exon_extend$intron_end <- as.numeric(intron_from_exon_extend$intron_end)
intron_from_exon_extend_regions <- GRanges(seqnames = Rle(chromosome),
                                    ranges = IRanges(start = intron_from_exon_extend$intron_start, 
                                                     end = intron_from_exon_extend$intron_end),
                                    strand = Rle(rep("*", nrow(intron_from_exon_extend))))

# peaks in extend gene exon regions
export_temp <- vector()
for(name in read_list){
  exon_extend_peak_sum <- sum(countOverlaps(exon_extend_regions, eval(parse(text = name))))
  print(paste(name, "has", exon_extend_peak_sum, "peaks in exon extend regions", sep = " "))
  export_temp <- append(export_temp, exon_extend_peak_sum)
}
export[nrow(export)+1,] <- export_temp
rownames(export)[nrow(export)] <- "exon extend regions"

# peaks in gene intron regions from extend exon
export_temp <- vector()
for(name in read_list){
  intron_peak_sum <- sum(countOverlaps(intron_from_exon_extend_regions, eval(parse(text = name))))
  print(paste(name, "has", intron_peak_sum, "peaks in intron from exon extend regions", sep = " "))
  export_temp <- append(export_temp, intron_peak_sum)
}
export[nrow(export)+1,] <- export_temp
rownames(export)[nrow(export)] <- "intron from exon extend regions"

# intron extend regions and exon from intron extend regions
intron_extend <- data.frame("external_gene_name", "intron_start", "intron_end", stringsAsFactors = FALSE)
exon_from_intron_extend <- data.frame("external_gene_name", "exon_start", "exon_end", stringsAsFactors = FALSE)
for(name in unique(exon$external_gene_name)){
  name_index_list <- which(exon$external_gene_name == name)
  name_temp <- exon[name_index_list,]
  name_temp <- name_temp[order(name_temp$exon_chrom_start),]
  for(i in c(1:nrow(name_temp))){
    if(i == 1){
      start <- name_temp[i,4]
      end <- name_temp[i,5]
      intron_extend_list <- c(start, end)
    }
    else{
      if(name_temp[i,4] > end){
        start <- name_temp[i,4]
        end <- name_temp[i,5]
        intron_extend_list <- append(intron_extend_list, c(start, end))
      } else{
        if(name_temp[i,5] <= end){
          start <- name_temp[i,4]
          end <- name_temp[i,5]
          intron_extend_list[length(intron_extend_list)-1] <- start
          intron_extend_list[length(intron_extend_list)] <- end
        } else{
          start <- name_temp[i,4]
          intron_extend_list[length(intron_extend_list)-1] <- start
        }
      }
    }
  }
  for(i in c(1:(length(intron_extend_list)/2))){
    exon_from_intron_extend[nrow(exon_from_intron_extend)+1,] <- c(name, intron_extend_list[i*2-1], intron_extend_list[i*2])
  }
  if(length(intron_extend_list) >= 4){
    for(i in c(1:(length(intron_extend_list)/2-1))){
      intron_extend[nrow(intron_extend)+1,] <- c(name, intron_extend_list[i*2]+1, intron_extend_list[i*2+1]-1)
    }
  }
}
colnames(intron_extend) <- intron_extend[1,]
intron_extend <- intron_extend[c(2:nrow(intron_extend)),]
intron_extend$intron_start <- as.numeric(intron_extend$intron_start)
intron_extend$intron_end <- as.numeric(intron_extend$intron_end)
intron_extend_regions <- GRanges(seqnames = Rle(chromosome),
                                ranges = IRanges(start = intron_extend$intron_start, end = intron_extend$intron_end),
                                strand = Rle("*"))
colnames(exon_from_intron_extend) <- exon_from_intron_extend[1,]
exon_from_intron_extend <- exon_from_intron_extend[c(2:nrow(exon_from_intron_extend)),]
exon_from_intron_extend$exon_start <- as.numeric(exon_from_intron_extend$exon_start)
exon_from_intron_extend$exon_end <- as.numeric(exon_from_intron_extend$exon_end)
exon_from_intron_extend_regions <- GRanges(seqnames = Rle(chromosome),
                                    ranges = IRanges(start = exon_from_intron_extend$exon_start, 
                                                     end = exon_from_intron_extend$exon_end),
                                    strand = Rle("*"))

# peaks in gene intron extend regions 
export_temp <- vector()
for(name in read_list){
  intron_peak_sum <- sum(countOverlaps(intron_extend_regions, eval(parse(text = name))))
  print(paste(name, "has", intron_peak_sum, "peaks in intron extend regions", sep = " "))
  export_temp <- append(export_temp, intron_peak_sum)
}
export[nrow(export)+1,] <- export_temp
rownames(export)[nrow(export)] <- "intron extend regions"

# peaks in gene exon regions from intron extend
export_temp <- vector()
for(name in read_list){
  exon_peak_sum <- sum(countOverlaps(exon_from_intron_extend_regions, eval(parse(text = name))))
  print(paste(name, "has", exon_peak_sum, "peaks in exon from intron extend regions", sep = " "))
  export_temp <- append(export_temp, exon_peak_sum)
}
export[nrow(export)+1,] <- export_temp
rownames(export)[nrow(export)] <- "exon from intron extend regions"

# 5' UTR isolation and setting
utr_5 <- na.omit(getBM(attributes = c("ensembl_gene_id","external_gene_name",
                            "chromosome_name", "5_utr_start",
                            "5_utr_end", "strand"), mart = ds,
                            filter = "chromosome_name", 
                            values = strsplit(chromosome, split = "r")[[1]][2]))
utr_5_regions <- GRanges(seqnames = Rle(chromosome),
                         ranges = IRanges(start = utr_5$'5_utr_start', end = utr_5$'5_utr_end'),
                         strand = Rle("*"),
                         gene = utr_5$external_gene_name)

# 5' UTR extend
utr_5_extend <- data.frame("external_gene_name", "start", "end", stringsAsFactors = FALSE)
for(name in unique(utr_5$external_gene_name)){
  name_index_list <- which(utr_5$external_gene_name == name)
  name_temp <- utr_5[name_index_list,]
  name_temp <- name_temp[order(name_temp$"5_utr_start"),]
  for(i in c(1, length(name_index_list))){
    if(i == 1){
      start <- name_temp[i,]$"5_utr_start"
      end <- name_temp[i,]$"5_utr_end"
      utr_5_list <- c(start, end)
    } else{
      if(name_temp[i,]$"5_utr_start" <= start){
        if(name_temp[i,]$"5_utr_end" >= end){
          start <- name_temp[i,]$"5_utr_start"
          end <- name_temp[i,]$"5_utr_end"
          utr_5_list <- c(start, end)
        } else{
          start <- name_temp[i,]$"5_utr_start"
          utr_5_list[1] <- start
        }
      } else{
        if(name_temp[i,]$"5_utr_end" >= end){
          start <- name_temp[i,]$"5_utr_start"
          end <- name_temp[i,]$"5_utr_end"
          utr_5_list[2] <- end
        }
      }
    }
  }
  utr_5_extend[nrow(utr_5_extend)+1,] <- c(name, utr_5_list)
}
colnames(utr_5_extend) <- utr_5_extend[1,]
utr_5_extend <- utr_5_extend[c(2:nrow(utr_5_extend)),]
utr_5_extend$start <- as.numeric(utr_5_extend$start)
utr_5_extend$end <- as.numeric(utr_5_extend$end)
utr_5_extend_regions <- GRanges(seqnames = Rle(chromosome),
                                ranges = IRanges(start = utr_5_extend$start, 
                                                 end = utr_5_extend$end),
                                strand = Rle("*"),
                                gene = Rle(utr_5_extend$external_gene_name))

# peaks in 5' UTR extend regions
export_temp <- vector()
for(name in read_list){
  utr_5_peak_sum <- sum(countOverlaps(utr_5_extend_regions, eval(parse(text = name))))
  print(paste(name, "has", utr_5_peak_sum, "peaks in 5' UTR extend regions", sep = " "))
  export_temp <- append(export_temp, utr_5_peak_sum)
}
export[nrow(export)+1,] <- export_temp
rownames(export)[nrow(export)] <- "5' UTR extend regions"

# 5' UTR overlap
utr_5_overlap <- data.frame("chromosome_name", "start", "end", stringsAsFactors = FALSE)
utr_5_extend <- utr_5_extend[order(utr_5_extend$start),]
for(i in c(1:nrow(utr_5_extend))){
  if(i == 1){
    start <- utr_5_extend[i,]$start
    end <- utr_5_extend[i,]$end
    utr_5_overlap_list <- c(start, end)
  } else{
    if(utr_5_extend[i,]$start > end){
      start <- utr_5_extend[i,]$start
      end <- utr_5_extend[i,]$end
      utr_5_overlap_list <- append(utr_5_overlap_list, c(start, end))
    } else{
      if(utr_5_extend[i,]$end >= end){
        end <- utr_5_extend[i,]$end
        utr_5_overlap_list[length(utr_5_overlap_list)] <- end
      }
    }
  }
}
for(i in c(1:(length(utr_5_overlap_list)/2))){
  utr_5_overlap[nrow(utr_5_overlap)+1,] <- c(chromosome, utr_5_overlap_list[i*2-1], utr_5_overlap_list[i*2])
}
colnames(utr_5_overlap) <- utr_5_overlap[1,]
utr_5_overlap <- utr_5_overlap[c(2:nrow(utr_5_overlap)),]
utr_5_overlap$start <- as.numeric(utr_5_overlap$start)
utr_5_overlap$end <- as.numeric(utr_5_overlap$end)
utr_5_overlap_regions <- GRanges(seqnames = Rle(chromosome),
                                ranges = IRanges(start = utr_5_overlap$start, 
                                                 end = utr_5_overlap$end),
                                strand = Rle("*"))

# peaks in 5' UTR overlap regions
export_temp <- vector()
for(name in read_list){
  utr_5_peak_sum <- sum(countOverlaps(utr_5_overlap_regions, eval(parse(text = name))))
  print(paste(name, "has", utr_5_peak_sum, "peaks in 5' UTR overlap regions", sep = " "))
  export_temp <- append(export_temp, utr_5_peak_sum)
}
export[nrow(export)+1,] <- export_temp
rownames(export)[nrow(export)] <- "5' UTR overlap regions"

# 3' UTR isolation and setting
utr_3 <- na.omit(getBM(attributes = c("ensembl_gene_id","external_gene_name",
                              "chromosome_name", "3_utr_start",
                              "3_utr_end", "strand"), mart = ds,
                              filter = "chromosome_name", 
                              values = strsplit(chromosome, split = "r")[[1]][2]))
utr_3_regions <- GRanges(seqnames = Rle(chromosome),
                         ranges = IRanges(start = utr_3$'3_utr_start', end = utr_3$'3_utr_end'),
                         strand = Rle("*"),
                         gene = utr_3$external_gene_name)

# 3' UTR extend
utr_3_extend <- data.frame("external_gene_name", "start", "end", stringsAsFactors = FALSE)
for(name in unique(utr_3$external_gene_name)){
  name_index_list <- which(utr_3$external_gene_name == name)
  name_temp <- utr_3[name_index_list,]
  name_temp <- name_temp[order(name_temp$"3_utr_start"),]
  for(i in c(1, length(name_index_list))){
    if(i == 1){
      start <- name_temp[i,]$"3_utr_start"
      end <- name_temp[i,]$"3_utr_end"
      utr_3_list <- c(start, end)
    } else{
      if(name_temp[i,]$"3_utr_start" <= start){
        if(name_temp[i,]$"3_utr_end" >= end){
          start <- name_temp[i,]$"3_utr_start"
          end <- name_temp[i,]$"3_utr_end"
          utr_3_list <- c(start, end)
        } else{
          start <- name_temp[i,]$"3_utr_start"
          utr_3_list[1] <- start
        }
      } else{
        if(name_temp[i,]$"3_utr_end" >= end){
          start <- name_temp[i,]$"3_utr_start"
          end <- name_temp[i,]$"3_utr_end"
          utr_3_list[2] <- end
        }
      }
    }
  }
  utr_3_extend[nrow(utr_3_extend)+1,] <- c(name, utr_3_list)
}
colnames(utr_3_extend) <- utr_3_extend[1,]
utr_3_extend <- utr_3_extend[c(2:nrow(utr_3_extend)),]
utr_3_extend$start <- as.numeric(utr_3_extend$start)
utr_3_extend$end <- as.numeric(utr_3_extend$end)
utr_3_extend_regions <- GRanges(seqnames = Rle(chromosome),
                                ranges = IRanges(start = utr_3_extend$start, 
                                                 end = utr_3_extend$end),
                                strand = Rle("*"),
                                gene = Rle(utr_3_extend$external_gene_name))

# peaks in 3' extend UTR regions
export_temp <- vector()
for(name in read_list){
  utr_3_peak_sum <- sum(countOverlaps(utr_3_extend_regions, eval(parse(text = name))))
  print(paste(name, "has", utr_3_peak_sum, "peaks in 3' UTR extend regions", sep = " "))
  export_temp <- append(export_temp, utr_3_peak_sum)
}
export[nrow(export)+1,] <- export_temp
rownames(export)[nrow(export)] <- "3' UTR extend regions"

# 3' UTR overlap
utr_3_overlap <- data.frame("chromosome_name", "start", "end", stringsAsFactors = FALSE)
utr_3_extend <- utr_3_extend[order(utr_3_extend$start),]
for(i in c(1:nrow(utr_3_extend))){
  if(i == 1){
    start <- utr_3_extend[i,]$start
    end <- utr_3_extend[i,]$end
    utr_3_overlap_list <- c(start, end)
  } else{
    if(utr_3_extend[i,]$start > end){
      start <- utr_3_extend[i,]$start
      end <- utr_3_extend[i,]$end
      utr_3_overlap_list <- append(utr_3_overlap_list, c(start, end))
    } else{
      if(utr_3_extend[i,]$end >= end){
        end <- utr_3_extend[i,]$end
        utr_3_overlap_list[length(utr_3_overlap_list)] <- end
      }
    }
  }
}
for(i in c(1:(length(utr_3_overlap_list)/2))){
  utr_3_overlap[nrow(utr_3_overlap)+1,] <- c(chromosome, utr_3_overlap_list[i*2-1], utr_3_overlap_list[i*2])
}
colnames(utr_3_overlap) <- utr_3_overlap[1,]
utr_3_overlap <- utr_3_overlap[c(2:nrow(utr_3_overlap)),]
utr_3_overlap$start <- as.numeric(utr_3_overlap$start)
utr_3_overlap$end <- as.numeric(utr_3_overlap$end)
utr_3_overlap_regions <- GRanges(seqnames = Rle(chromosome),
                                ranges = IRanges(start = utr_3_overlap$start, 
                                                 end = utr_3_overlap$end),
                                strand = Rle("*"))

# peaks in 3' overlap UTR regions
export_temp <- vector()
for(name in read_list){
  utr_3_peak_sum <- sum(countOverlaps(utr_3_overlap_regions, eval(parse(text = name))))
  print(paste(name, "has", utr_3_peak_sum, "peaks in 3' UTR overlap regions", sep = " "))
  export_temp <- append(export_temp, utr_3_peak_sum)
}
export[nrow(export)+1,] <- export_temp
rownames(export)[nrow(export)] <- "3' UTR overlap regions"

# intragenic regions isolation from egs overlap (gene distance > 200kb)
egs_overlap <- egs_overlap[order(egs_overlap$start),]
gd_list <- c(1)
for(i in c(1:nrow(egs_overlap))){
  if(i == 1){
    start <- egs_overlap[i,]$star-100000
    end <- egs_overlap[i,]$end+100000
    gd_list <- append(gd_list, c(start, end))
  }
  else{
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
gd_list <- append(gd_list, seqlengths(genome)[as.numeric(strsplit(chromosome, split = "r")[[1]][2])])
gd_temp <- data.frame("chromosome_name", "gd_start", "gd_end", stringsAsFactors=FALSE)
for(i in c(ifelse(gd_list[2] - gd_list[1] > 100000, 1, 2):ifelse(gd_list[length(gd_list)] - gd_list[length(gd_list)-1] > 100000, length(gd_list)/2, length(gd_list)/2-1))){
  gd_temp[nrow(gd_temp)+1,] <- c(chromosome, gd_list[i*2-1], gd_list[i*2])
}
colnames(gd_temp) <- gd_temp[1,]
gd_temp <- gd_temp[c(2:nrow(gd_temp)),]
gd_temp$gd_start <- as.numeric(gd_temp$gd_start)
gd_temp$gd_end <- as.numeric(gd_temp$gd_end)
gd_regions <- GRanges(seqnames = Rle(chromosome), 
                      ranges = IRanges(start = gd_temp$gd_start, end = gd_temp$gd_end))

# peaks in intragenic regions (gene distance > 200kb)
export_temp <- vector()
for(name in read_list){
  intragene_peak_sum <- sum(countOverlaps(gd_regions, eval(parse(text = name))))
  print(paste(name, "has", intragene_peak_sum, "peaks in intragenic regions", sep = " "))
  export_temp <- append(export_temp, intragene_peak_sum)
}
export[nrow(export)+1,] <- export_temp
rownames(export)[nrow(export)] <- "intragenic regions"

# create data export file (2)
colnames(export) <- export[1,]
export <- export[c(2:nrow(export)),]
write.csv(export, paste0(path, "/export.csv"))
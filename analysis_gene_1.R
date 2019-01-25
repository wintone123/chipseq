# load separated bed files

# load library
library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(biomaRt)
library(chipseq)
library(BSgenome.Mmusculus.UCSC.mm10)  
library(Gviz)
library(dplyr)

print("================Let's go================")
# load necessary data
genome <- BSgenome.Mmusculus.UCSC.mm10
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "mmusculus_gene_ensembl",  # human (hsapiens_gene_ensembl)
                host = "jul2018.archive.ensembl.org")  # ensembl release 93 (july 2018)
ds <- useDataset("mmusculus_gene_ensembl", mart = mart)
prepareChipseq <- function(reads){
  frag_len = median(estimate.mean.fraglen(reads))
  reads_extended = resize(reads, width = frag_len)
  return(trim(reads_extended))
}
path <- "/mnt/c/chipseq/test2"  # input file path
load_list <- list.files(file.path(path, "split_files"))
dir.create(file.path(path,"output"), showWarnings = FALSE)  # create output folder

# gene list
gene_list <- c("Nanog","Sox2","Pou5f1")
gene_info_df <- getBM(attributes = c("ensembl_gene_id","external_gene_name", "chromosome_name", "strand",
                                     "start_position", "end_position", "exon_chrom_start", "exon_chrom_end",
                                     "5_utr_start", "5_utr_end", "3_utr_start", "3_utr_end"), 
                      mart = ds, filter = "external_gene_name", values = gene_list)
gene_df <- gene_info_df[,c(2,3)]
gene_df <- gene_df[!duplicated(gene_df),]
gene_df <- gene_df[order(gene_df$chromosome_name),]
for(i in c(1:nrow(gene_df))){
  assign(paste0(gene_df[i,]$external_gene_name, "_info"), 
         filter(gene_info_df, external_gene_name == gene_df[i,]$external_gene_name))
}
print("database loaded")
 
# analysis process
for(gene_num in c(1:nrow(gene_df))){
  print("+++++++++++++++Loop Start+++++++++++++++")
  # load info
  gene <- gene_df[gene_num,]$external_gene_name
  print(paste0("target gene: ", gene))
  gene_name <- paste0(gene, "_info")
  chrom <- gene_df[gene_num,]$chromosome_name
  
  # isolate promoter (+/- 200 bp from TSS)
  promoter <- eval(parse(text = gene_name))[c(1:6)]
  promoter <- promoter[!duplicated(promoter),]
  if(promoter$strand == 1){
    start <- promoter$start_position - 200
    end <- promoter$start_position + 200
  } else {
    start <- promoter$end_position - 200
    end <- promoter$end_position + 200
  }
  promoter_region <- GRanges(seqnames = Rle(paste0("chr", chrom)),
                             ranges = IRanges(start = start, end = end),
                             strand = Rle("*"))

  # isolate extend 5' UTR 
  utr_5 <- eval(parse(text = gene_name))[c(1:4,9,10)]
  utr_5 <- na.omit(utr_5[!duplicated(utr_5),])
  utr_5 <- utr_5[order(utr_5$"5_utr_start"),]
  for(i in c(1:nrow(utr_5))){
    if(i == 1){
      start <- utr_5[i,]$"5_utr_start"
      end <- utr_5[i,]$"5_utr_end"
      utr_5_list <- c(start, end)
    } else{
      if(utr_5[i,]$"5_utr_start" < start){
        start <- utr_5[i,]$"5_utr_start"
        utr_5_list[1] <- start
      } 
      if(utr_5[i,]$"5_utr_end" > end){
        end <- utr_5[i,]$"5_utr_end"
        utr_5_list[2] <- end
      }
    }
  }
  utr_5_regin <- GRanges(seqnames = Rle(paste0("chr", chrom)),
                         ranges = IRanges(start = utr_5_list[1], end = utr_5_list[2]),
                         strand = Rle("*"))
  
  # isolate extend 3' UTR
  utr_3 <- eval(parse(text = gene_name))[c(1:4,11,12)]
  utr_3 <- na.omit(utr_3[!duplicated(utr_3),])
  utr_3 <- utr_3[order(utr_3$"3_utr_start"),]
  for(i in c(1:nrow(utr_3))){
    if(i == 1){
      start <- utr_3[i,]$"3_utr_start"
      end <- utr_3[i,]$"3_utr_end"
      utr_3_list <- c(start, end)
    } else{
      if(utr_3[i,]$"3_utr_start" < start){
        start <- utr_3[i,]$"3_utr_start"
        utr_3_list[1] <- start
      } 
      if(utr_3[i,]$"3_utr_end" > end){
        end <- utr_3[i,]$"3_utr_end"
        utr_3_list[2] <- end
      }
    }
  }
  utr_3_regin <- GRanges(seqnames = Rle(paste0("chr", chrom)),
                         ranges = IRanges(start = utr_3_list[1], end = utr_3_list[2]),
                         strand = Rle("*"))

  # isolate extend exon and intron from extend exon
  exon <- eval(parse(text = gene_name))[c(1:4,7,8)]
  exon <- na.omit(exon[!duplicated(exon),])
  exon <- exon[order(exon$exon_chrom_start),]
  for(i in c(1:nrow(exon))){
    if(i == 1){
      start <- exon[i,]$exon_chrom_start
      end <- exon[i,]$exon_chrom_end
      exon_list <- c(start, end)
    } else{
      if(exon[i,]$exon_chrom_start > end){
        start <- exon[i,]$exon_chrom_start
        end <- exon[i,]$exon_chrom_end
        exon_list <- append(exon_list, c(start, end)) 
      } else{
        if(exon[i,]$exon_chrom_end >= end){
          end <- exon[i,]$exon_chrom_end
          exon_list[length(exon_list)] <- end
        }
      }
    }
  }
  for(i in c(1:(length(exon_list)/2))){
    if(i == 1){
      exon_extend <- data.frame("external_gene_name" = gene, "start" = exon_list[i*2-1], "end" = exon_list[i*2])
    } else{
      exon_extend[nrow(exon_extend)+1,] <- c(gene, exon_list[i*2-1], exon_list[i*2])
    }
  }
  exon_extend$start <- as.numeric(exon_extend$start)
  exon_extend$end <- as.numeric(exon_extend$end)
  exon_extend_regin <- GRanges(seqnames = Rle(paste0("chr", chrom)),
                               ranges = IRanges(start = exon_extend$start, end = exon_extend$end),
                               strand = Rle("*"))
  if(length(exon_list) >= 4){
    for(i in c(1:(length(exon_list)/2-1))){
      if(i == 1){
        intron_from_exon_extend <- data.frame("external_gene_name" = gene, "start" = exon_list[i*2]+1, "end" = exon_list[i*2+1]-1)
      } else{
        intron_from_exon_extend[nrow(intron_from_exon_extend)+1,] <- c(gene, exon_list[i*2]+1, exon_list[i*2+1]-1)
      }
    }
    intron_from_exon_extend$start <- as.numeric(intron_from_exon_extend$start)
    intron_from_exon_extend$end <- as.numeric(intron_from_exon_extend$end)
    intron_from_exon_extend_regin <- GRanges(seqnames = Rle(paste0("chr", chrom)),
                                             ranges = IRanges(start = intron_from_exon_extend$start, 
                                                              end = intron_from_exon_extend$end),
                                             strand = Rle("*"))
  } else{
    intron_from_exon_extend_regin <- GRanges(seqnames = Rle(paste0("chr", chrom)),
                                             ranges = IRanges(start = 0, end = 0),
                                             strand = Rle("*"))
  }
  
  # isolate extend intron and exon from extend intron
  for(i in c(1:nrow(exon))){
    if(i == 1){
      start <- exon[i,]$exon_chrom_start
      end <- exon[i,]$exon_chrom_end
      exon_list_2 <- c(start, end)
    } else{
      if(exon[i,]$exon_chrom_start > end){
        start <- exon[i,]$exon_chrom_start
        end <- exon[i,]$exon_chrom_end
        exon_list_2 <- append(exon_list_2, c(start, end))
      } else {
        if(exon[i,]$exon_chrom_start <= end){
          start <- exon[i,]$exon_chrom_start
          end <- exon[i,]$exon_chrom_end
          exon_list_2[length(exon_list_2)] <- start
          exon_list_2[length(exon_list_2)] <- end
        } else{
          start <- exon[i,]$exon_chrom_start
          exon_list_2[length(exon_list_2)-1] <- start
          }
        }
      }
    }

  for(i in c(1:(length(exon_list_2)/2))){
    if(i == 1){
      exon_from_intron_extend <- data.frame("external_gene_name" = gene, "start" = exon_list_2[i*2-1], "end" = exon_list_2[i*2])
    } else{
      exon_from_intron_extend[nrow(exon_from_intron_extend)+1,] <- c(gene, exon_list_2[i*2-1], exon_list_2[i*2])
    }
  }
  exon_from_intron_extend$start <- as.numeric(exon_from_intron_extend$start)
  exon_from_intron_extend$end <- as.numeric(exon_from_intron_extend$end)
  exon_from_intron_extend_regin <- GRanges(seqnames = Rle(paste0("chr", chrom)),
                                           ranges = IRanges(start = exon_from_intron_extend$start, 
                                           end = exon_from_intron_extend$end),
                                           strand = Rle("*"))
  if(length(exon_list_2) >= 4){
    for(i in c(1:(length(exon_list_2)/2-1))){
      if(i == 1){
        intron_extend <- data.frame("external_gene_name" = gene, "start" = exon_list_2[i*2]+1, "end" = exon_list_2[i*2+1]-1)
      } else{
        intron_extend[nrow(intron_extend)+1,] <- c(gene, exon_list_2[i*2]+1, exon_list_2[i*2+1]-1)
      }
    }
    intron_extend$start <- as.numeric(intron_extend$start)
    intron_extend$end <- as.numeric(intron_extend$end)
    intron_extend_regin <- GRanges(seqnames = Rle(paste0("chr", chrom)),
                                   ranges = IRanges(start = intron_extend$start, 
                                                              end = intron_extend$end),
                                   strand = Rle("*"))
  } else{
    intron_extend_regin <- GRanges(seqnames = Rle(paste0("chr", chrom)),
                                   ranges = IRanges(start = 0, end = 0),
                                   strand = Rle("*"))
  }

  #load files
  if(gene_num == 1){
    read_list <- vector()
    for(file_name in load_list){
      if(strsplit(file_name, split = "_", fixed = TRUE)[[1]][3] == paste0("chr", chrom, ".bed")){
        base_name <- strsplit(file_name, split = "_", fixed = TRUE)[[1]][1]
        assign(base_name, prepareChipseq(import.bed(file.path(path, "split_files", file_name))))
        read_list <- c(read_list, base_name)
      }
    }
    read_list <- sort(read_list)
  } else {
    if(chrom != gene_df[gene_num-1,]$chromosome_name){
      for(name in read_list){
        remove(name)
      }
      read_list <- vector()
      for(file_name in load_list){
        if(strsplit(file_name, split = "_", fixed = TRUE)[[1]][3] == paste0("chr", chrom, ".bed")){
          base_name <- strsplit(file_name, split = "_", fixed = TRUE)[[1]][1]
          assign(base_name, prepareChipseq(import.bed(file.path(path, "split_files", file_name))))
          read_list <- c(read_list, base_name)
        }
      }
      read_list <- sort(read_list)
    }
  }
  print("bed file loaded")

  # findoverlap 
  if(gene_num == 1){
    output <- data.frame(c("total peaks", "promoter peaks", "5' UTR peaks", "3' UTR peaks", "exon extend peaks", 
                           "exon from intron extend peaks", "intron extend peaks", "intron from exon extend peaks"))
    output[,ncol(output)+1] <- c("", rep(gene, nrow(output)-1))
    output[,ncol(output)+1] <- c(rep(paste0("chr", chrom), nrow(output)))
    for(name in read_list){
      print(paste0("sample: ", name))
      total_peak <- length(eval(parse(text = name)))
      print(paste0("total peaks: ", total_peak))
      promoter_peak <- sum(countOverlaps(eval(parse(text = name)), promoter_region))
      print(paste0("promoter peaks: ", promoter_peak))
      utr_5_peak <- sum(countOverlaps(eval(parse(text = name)), utr_5_regin))
      print(paste0("5' UTR peaks: ", utr_5_peak))
      utr_3_peak <- sum(countOverlaps(eval(parse(text = name)), utr_3_regin))
      print(paste0("3' UTR peaks: ", utr_3_peak))
      exon_extend_peak <- sum(countOverlaps(eval(parse(text = name)), exon_extend_regin))
      print(paste0("exon extend peaks: ", exon_extend_peak))
      exon_from_intron_extend_peak <- sum(countOverlaps(eval(parse(text = name)), exon_from_intron_extend_regin))
      print(paste0("exon from intron extend peaks: ", exon_from_intron_extend_peak))
      intron_extend_peak <- sum(countOverlaps(eval(parse(text = name)), intron_extend_regin))
      print(paste0("intron extend peaks: ", intron_extend_peak))
      intron_from_exon_extend_peak <- sum(countOverlaps(eval(parse(text = name)), intron_from_exon_extend_regin))
      print(paste0("intron from exon extend peaks: ", intron_from_exon_extend_peak))
      output[,ncol(output)+1] <- c(total_peak, promoter_peak, utr_5_peak, utr_3_peak, exon_extend_peak, 
                                  exon_from_intron_extend_peak, intron_extend_peak, intron_from_exon_extend_peak)
      print("---------------------------------------")
    }
    colnames(output) <- c("items", "gene", "chrom", read_list)
  } else{
    output_1 <- data.frame(c("total peaks", "promoter peaks", "5' UTR peaks", "3' UTR peaks", "exon extend peaks", 
                             "exon from intron extend peaks", "intron extend peaks", "intron from exon extend peaks"))
    output_1[,ncol(output_1)+1] <- c("", rep(gene, nrow(output_1)-1))
    output_1[,ncol(output_1)+1] <- c(rep(paste0("chr", chrom), nrow(output_1)))
    for(name in read_list){
      print(paste0("sample: ", name))
      total_peak <- length(eval(parse(text = name)))
      print(paste0("total peaks: ", total_peak))
      promoter_peak <- sum(countOverlaps(eval(parse(text = name)), promoter_region))
      print(paste0("promoter peaks: ", promoter_peak))
      utr_5_peak <- sum(countOverlaps(eval(parse(text = name)), utr_5_regin))
      print(paste0("5' UTR peaks: ", utr_5_peak))
      utr_3_peak <- sum(countOverlaps(eval(parse(text = name)), utr_3_regin))
      print(paste0("3' UTR peaks: ", utr_3_peak))
      exon_extend_peak <- sum(countOverlaps(eval(parse(text = name)), exon_extend_regin))
      print(paste0("exon extend peaks: ", exon_extend_peak))
      exon_from_intron_extend_peak <- sum(countOverlaps(eval(parse(text = name)), exon_from_intron_extend_regin))
      print(paste0("exon from intron extend peaks: ", exon_from_intron_extend_peak))
      intron_extend_peak <- sum(countOverlaps(eval(parse(text = name)), intron_extend_regin))
      print(paste0("intron extend peaks: ", intron_extend_peak))
      intron_from_exon_extend_peak <- sum(countOverlaps(eval(parse(text = name)), intron_from_exon_extend_regin))
      print(paste0("intron from exon extend peaks: ", intron_from_exon_extend_peak))
      output_1[,ncol(output_1)+1] <- c(total_peak, promoter_peak, utr_5_peak, utr_3_peak, exon_extend_peak, 
                                       exon_from_intron_extend_peak, intron_extend_peak, intron_from_exon_extend_peak)
      print("---------------------------------------")
    }
    colnames(output_1) <- c("items", "gene", "chrom", read_list)
    output <- rbind(output, output_1)
  }    
  print("+++++++++++++++Loop Done+++++++++++++++")
}
# export data
write.csv(output, paste0(path, "/output/output.csv"))
print(paste0("output: /output/output.csv"))
print("================Finish!================")

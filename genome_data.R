# create data export file (1)
colname <- vector()
export <- promoter_peaks_export <- data.frame(t(unlist(read_list)), stringsAsFactors = FALSE)
export_temp <- vector()
for(name in read_list){
  export_temp <- append(export_temp, length(eval(parse(text = name))))
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
egs_regions <- GRanges(seqnames = Rle(chromosome)),
                      ranges = IRanges(start = egs$start_position, end = egs$end_position),
                      strand = Rle("*"),
                      gene = egs$external_gene_name)
egs_total_length <- sum(width(reduce(egs_regions)))

# peaks in egs regions
export_temp <- vector()
for(name in read_list){
  egs_peak_sum <- sum(countOverlaps(egs_regions, eval(parse(text = name))))
  print(paste(name, "has", egs_peak_sum, "in egs regions", sep = " "))
  export_temp <- append(export_temp, egs_peak_sum)
}
export[nrow(export)+1,] <- export_temp
rownames(export)[nrow(export)] <- "egs regions"

# peaks in certain egs
peaks_in_egs <- function(gene_name){
  egs_index <- match(gene_name, egs_regions$gene)
  for(name in read_list){
    index <- match(name, read_list)
    peak_sum <- sum(countOverlaps(egs_regions[egs_index], eval(parse(text = name)))) / fold_list[index]
    print(paste(name, "has", peak_sum, "peaks in", gene_name, "egs", sep = " "))
  }
}
peaks_in_egs("")

# egs overlapping regions
egs <- egs[order(egs$start_position),]
for(i in c(1: nrow(egs))){
  if(i == 1){
    start <- egs[i,]$start_position
    end <- egs[i,]$end_position
    egs_overlap_list <- c(start, end)
  }
  else{
    if(egs[i,]$start_position > end){
      start <- egs[i,]$start_position
      end <- egs[i,]$end_position
      egs_overlap_list <- append(egs_overlap_list, c(start, end))
    }
    else{
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
egs_overlap_regions <- GRanges(seqnames = Rle(chromosome),
                               ranges = IRanges(start = egs_overlap$start, end = egs_overlap$end),
                               strand = Rle("*"))

# peaks in egs overlapp regions
export_temp <- vector()
for(name in read_list){
  egs_overlap_peak_sum <- sum(countOverlaps(egs_overlap_regions, eval(parse(text = name))))
  print(paste(name, "has", egs_overlap_peak_sum, "in egs overlap regions", sep = " "))
  export_temp <- append(export_temp, egs_overlap_peak_sum)
}
export[nrow(export)+1,] <- export_temp
rownames(export)[nrow(export)] <- "egs overlap regions"

# set promoter (+/- 200 bp around TSS)
egs$TSS <- ifelse(egs$strand == "1", egs$start_position, egs$end_position)
promoter_regions <- GRanges(seqnames = Rle(chromosome),
                            ranges = IRanges(start = egs$TSS - 200, end = egs$TSS + 200),
                            strand = Rle("*"),
                            gene = egs$external_gene_name)
promoter_total_length <- sum(width(reduce(promoter_regions)))

# peaks in promoter regions (+/- 200 bp around TSS)
export_temp <- vector()
for(name in read_list){
  promoter_peak_sum <- sum(countOverlaps(promoter_regions, eval(parse(text = name))))
  print(paste(name, "has", promoter_peak_sum, "in promoter regions", sep = " "))
  export_temp <- append(export_temp, promoter_peak_sum)
}
export[nrow(export)+1,] <- export_temp
rownames(export)[nrow(export)] <- "promoter regions"

# peaks in promoter regions (+/- 200 bp around TSS, df export)
promoter_peaks_export <- data.frame(t(unlist(c("gene_name", read_list))), stringsAsFactors = FALSE)
for(gene_index in c(1:length(promoter_regions))){
  peak_sum_list <- c(promoter_regions[gene_index]$gene)
  for(name_index in c(1:length(read_list))){
    peak_sum <- sum(countOverlaps(promoter_regions[gene_index], 
                                  eval(parse(text = read_list[name_index])))) / fold_list[name_index]
    peak_sum_list <- append(peak_sum_list, peak_sum)
  }
  promoter_peaks_export[nrow(promoter_peaks_export)+1,] <- peak_sum_list
}
colnames(promoter_peaks_export) <- promoter_peaks_export[1,]
promoter_peaks_export <- promoter_peaks_export[c(2:nrow(promoter_peaks_export)),]

# peaks in certain promoter with normalization (+/- 200 bp around TSS)
peaks_in_promoter <- function(gene_name){
  promoter_index <- match(gene_name, promoter_regions$gene)
  for(name in read_list){
    index <- match(name, read_list)
    peak_sum <- sum(countOverlaps(promoter_regions[promoter_index], eval(parse(text = name)))) / fold_list[index]
    print(paste(name, "has", peak_sum, "peaks in", gene_name, "promoter", sep = " "))
  }
}
peaks_in_promoter("")

# overlapping promoter with enriched regions
ovlp2 <- findOverlaps(enriched_regions, promoter_regions)
cat(sprintf("%d of %d promoters are overlapped by an enriched region.",
            length(unique(subjectHits(ovlp2))), length(promoter_regions)))
ovlp2b <- findOverlaps(promoter_regions, enriched_regions)
cat(sprintf("d% od d% enriched regions overlap a promoter.",
            length(unique(subjectHits(ovlp2))), length(enriched_regions)))
pos_TSS <- egs[unique(queryHits(findOverlaps(promoter_regions, enriched_regions))),]

# distribution of peaks around a subset of promoters
promoter_tiles <- sapply(1:nrow(pos_TSS), function(i)
    if(pos_TSS$strand[i] == "1")
      pos_TSS$TSS[i] + seq(-1000, 900, length.out = 20)
    else 
      pos_TSS$TSS[i] + seq(900, -1000, length.out = 20))
promoter_tiles <- GRanges(tilename = paste(rep(pos_TSS$ensembl_gene_id, each = 20), 1:20, sep = "-"),
                          seqnames = Rle(rep(paste0("chr", pos_TSS$chromosome_name), each = 20)),
                          ranges = IRanges(start = as.vector(promoter_tiles), width = 100),
                          strand = Rle(rep("*", length(as.vector(promoter_tiles)))),
                          seqinfo = si)

# peaks in certain gene with normalization
peaks_in_gene <- function(gene_name){
  index <- match(gene_name, egs$external_gene_name)
  assign(gene_name, GRanges(seqnames = Rle(chromosome),
                            ranges = IRanges(start = egs[index,4], end = egs[index,5]),
                            strand = Rle("*"), 
                            gene = egs[index,2]))
  for(name in read_list){
    index <- match(name, read_list)
    peak_sum <- sum(countOverlaps(eval(parse(text = gene_name)), eval(parse(text = name)))) / fold_list[index]
    print(paste(name, "has", peak_sum, "peaks in", gene_name, "regions", sep = " "))
  }
}
peaks_in_gene("")

# exon isolation
exon <- getBM(attributes = c("ensembl_gene_id","external_gene_name",
                            "chromosome_name", "exon_chrom_start",
                            "exon_chrom_end", "strand"), mart = ds,
                            filter = "chromosome_name", 
                            values = strsplit(chromosome, split = "r")[[1]][2])

# peaks in exon regions
exon_regions <- GRanges(seqnames = Rle(chromosome),
                        ranges = IRanges(start = exon$exon_chrom_start, end = exon$exon_chrom_end),
                        srtand = Rle("*"),
                        gene = exon$external_gene_name)
export_temp <- vector()
for(name in read_list){
  exon_peak_sum <- sum(countOverlaps(exon_regions, eval(parse(text = name))))
  print(paste(name, "has", exon_peak_sum, "in exon regions", sep = " "))
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
    }
    else{
      if(name_temp[i,4] > end){
        start <- name_temp[i,4]
        end <- name_temp[i,5]
        exon_extend_list <- append(exon_extend_list, c(start, end))
      }
      else{
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
                                strand = Rle("*"))
colnames(intron_from_exon_extend) <- intron_from_exon_extend[1,]
intron_from_exon_extend <- intron_from_exon_extend[c(2:nrow(intron_from_exon_extend)),]
intron_from_exon_extend$intron_start <- as.numeric(intron_from_exon_extend$intron_start)
intron_from_exon_extend$intron_end <- as.numeric(intron_from_exon_extend$intron_end)
intron_from_exon_extend_regions <- GRanges(seqnames = Rle(chromosome),
                                    ranges = IRanges(start = intron_from_exon_extend$intron_start, 
                                                     end = intron_from_exon_extend$intron_end),
                                    strand = Rle("*"))

# peaks in certain extend gene exon with normalization 
peaks_in_gene_exon <- function(gene_name){
  index_list <- which(exon_extend$external_gene_name == gene_name)
  assign(gene_name, GRanges(seqnames = Rle(chromosome),
                            ranges = IRanges(start = exon_extend[index_list,]$exon_start, end = exon_extend[index_list,]$exon_end),
                            strand = Rle("*"), 
                            gene = exon_extend[index_list,]$external_gene_name))
  for(name in read_list){
    index <- match(name, read_list)
    peak_sum <- sum(countOverlaps(eval(parse(text = gene_name)), eval(parse(text = name)))) / fold_list[index]
    print(paste(name, "has", peak_sum, "peaks in", gene_name, "exon regions", sep = " "))
  }
}
peaks_in_gene_exon("")

# peaks in certain gene intron from extend exon with normalization 
peaks_in_gene_intron <- function(gene_name){
  index_list <- which(intron_from_exon_extend$external_gene_name == gene_name)
  assign(gene_name, GRanges(seqnames = Rle(chromosome),
                            ranges = IRanges(start = intron_from_exon_extend[index_list,]$intron_start,
                                             end = intron_from_exon_extend[index_list,]$intron_end),
                            strand = Rle("*"), 
                            gene = intron_from_exon_extend[index_list,]$external_gene_name))
  for(name in read_list){
    index <- match(name, read_list)
    peak_sum <- sum(countOverlaps(eval(parse(text = gene_name)), eval(parse(text = name)))) / fold_list[index]
    print(paste(name, "has", peak_sum, "peaks in", gene_name, "exon regions", sep = " "))
  }
}
peaks_in_gene_exon("")

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
      }
      else{
        if(name_temp[i,5] <= end){
          start <- name_temp[i,4]
          end <- name_temp[i,5]
          intron_extend_list[length(intron_extend_list)-1] <- start
          intron_extend_list[length(intron_extend_list)] <- end
        }
        else{
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
utr_5_total_length <- sum(width(reduce(utr_5_regions)))
pos_utr_5 <- utr_5[unique(queryHits(findOverlaps(utr_5_regions, enriched_regions))),]

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
    }
    else{
      if(name_temp[i,]$"5_utr_start" <= start){
        if(name_temp[i,]$"5_utr_end" >= end){
          start <- name_temp[i,]$"5_utr_start"
          end <- name_temp[i,]$"5_utr_end"
          utr_5_list <- c(start, end)
        }
        else{
          start <- name_temp[i,]$"5_utr_start"
          utr_5_list[1] <- start
        }
      }
      else{
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
  }
  else{
    if(utr_5_extend[i,]$start > end){
      start <- utr_5_extend[i,]$start
      end <- utr_5_extend[i,]$end
      utr_5_overlap_list <- append(utr_5_overlap_list, c(start, end))
    }
    else{
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
utr_3_total_length <- sum(width(reduce(utr_3_regions)))
pos_utr_3 <- utr_5[unique(queryHits(findOverlaps(utr_3_regions, enriched_regions))),]

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
    }
    else{
      if(name_temp[i,]$"3_utr_start" <= start){
        if(name_temp[i,]$"3_utr_end" >= end){
          start <- name_temp[i,]$"3_utr_start"
          end <- name_temp[i,]$"3_utr_end"
          utr_3_list <- c(start, end)
        }
        else{
          start <- name_temp[i,]$"3_utr_start"
          utr_3_list[1] <- start
        }
      }
      else{
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
  }
  else{
    if(utr_3_extend[i,]$start > end){
      start <- utr_3_extend[i,]$start
      end <- utr_3_extend[i,]$end
      utr_3_overlap_list <- append(utr_3_overlap_list, c(start, end))
    }
    else{
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
    }
    else{
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
  print(paste(name, "has", intragene_peak_sum, "in intragenic regions", sep = " "))
  export_temp <- append(export_temp, intragene_peak_sum)
}
export[nrow(export)+1,] <- export_temp
rownames(export)[nrow(export)] <- "intragenic regions"

# create data export file (2)
colnames(export) <- export[1,]
export <- export[c(2:nrow(export)),]
write.csv(export, paste0(path, "/export.csv"))

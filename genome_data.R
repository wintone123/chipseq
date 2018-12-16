# egs isolation
listAttributes(mart)[c(1:5),]
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
egs_total_length <- sum(width(reduce(egs_regions)))

# remove egs overlapping regions


# peaks in egs regions
for(name in read_list){
  egs_peak_sum <- sum(countOverlaps(egs_regions, eval(parse(text = name))))
  print(paste(name, "has", egs_peak_sum, "in egs regions", sep = " "))
}

# set promoter (+/- 200 bp around TSS)
egs$TSS <- ifelse(egs$strand == "1", egs$start_position, egs$end_position)
promoter_regions <- GRanges(seqnames = Rle(paste0("chr", egs$chromosome_name)),
                            ranges = IRanges(start = egs$TSS - 200, end = egs$TSS + 200),
                            strand = Rle(rep("*", nrow(egs))),
                            gene = egs$external_gene_name)
promoter_total_length <- sum(width(reduce(promoter_regions)))

# peaks in promoter regions (+/- 200 bp around TSS)
for(name in read_list){
  promoter_peak_sum <- sum(countOverlaps(promoter_regions, eval(parse(text = name))))
  print(paste(name, "has", promoter_peak_sum, "in promoter regions", sep = " "))
}

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
  assign(gene_name, GRanges(seqnames = Rle(paste0("chr", egs[index,3])),
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
exon_regions <- GRanges(seqnames = Rle(paste0("chr", exon$chromosome_name)),
                        ranges = IRanges(start = exon$exon_chrom_start, end = exon$exon_chrom_end),
                        srtand = Rle(rep("*", nrow(exon))),
                        gene = exon$external_gene_name)
for(name in read_list){
  exon_peak_sum <- sum(countOverlaps(exon_regions, eval(parse(text = name))))
  print(paste(name, "has", exon_peak_sum, "in exon regions", sep = " "))
}

# peaks in certain gene exon with normalization
peaks_in_gene_exon <- function(gene_name){
  index_list <- which(exon$external_gene_name == gene_name)
  assign(gene_name, GRanges(seqnames = Rle(paste0("chr", exon[index_list,3])),
                            ranges = IRanges(start = exon[index_list,4], end = exon[index_list,5]),
                            strand = Rle("*",length(index_list)), 
                            gene = exon[index_list,2]))
  for(name in read_list){
    index <- match(name, read_list)
    peak_sum <- sum(countOverlaps(eval(parse(text = gene_name)), eval(parse(text = name)))) / fold_list[index]
    print(paste(name, "has", peak_sum, "peaks in", gene_name, "exon regions", sep = " "))
  }
}
peaks_in_gene_exon("")

# intron isolation from exon database
intron <- data.frame("ensembl_gene_id","external_gene_name","chromosome_name",
                 "strand", "intron_chrom_start","intron_chrom_end", stringsAsFactors=FALSE)
intron_list <- vector()
intron_temp_2 <- vector()
for(i in c(1:nrow(exon))){
  if(exon[i,6] == 1){
    if(i == 1){
      if(exon[i,2] == exon[i+1,2]){
        intron_list <- append(intron_list, exon[i,1])
        intron_list <- append(intron_list, exon[i,2])
        intron_list <- append(intron_list, exon[i,3])
        intron_list <- append(intron_list, exon[i,6])
        intron_list <- append(intron_list, exon[i,4])
        intron_list <- append(intron_list, exon[i,5])
      }
    }
    else if(i == nrow(exon)){
      if(exon[i,2] == exon[i-1,2]){
        intron_list <- append(intron_list, exon[i,4])
        intron_list <- append(intron_list, exon[i,5])
        intron_temp_1 <- intron_list[1:4]
        intron_list <- intron_list[5:length(intron_list)]
        for(j in c(1:(length(intron_list)/2-1))){
          intron_temp_2 <- append(intron_temp_1, intron_list[j*2])
          intron_temp_2 <- append(intron_temp_2, intron_list[j*2+1])
          intron[nrow(intron) + 1,] <- intron_temp_2 
          intron_temp_2 <- vector()
        }
        intron_list <- vector()
      }
    }
    else{
      if(exon[i,2] != exon[i-1,2] & exon[i,2] == exon[i+1,2]){
        intron_list <- append(intron_list, exon[i,1])
        intron_list <- append(intron_list, exon[i,2])
        intron_list <- append(intron_list, exon[i,3])
        intron_list <- append(intron_list, exon[i,6])
        intron_list <- append(intron_list, exon[i,4])
        intron_list <- append(intron_list, exon[i,5])
      }
      else if(exon[i,2] == exon[i-1,2] & exon[i,2] == exon[i+1,2]){
        intron_list <- append(intron_list, exon[i,4])
        intron_list <- append(intron_list, exon[i,5])
      }
      else if(exon[i,2] == exon[i-1,2] & exon[i,2] != exon[i+1,2]){
        intron_list <- append(intron_list, exon[i,4])
        intron_list <- append(intron_list, exon[i,5])
        intron_temp_1 <- intron_list[1:4]
        intron_list <- intron_list[5:length(intron_list)]
        for(j in c(1:(length(intron_list)/2-1))){
          intron_temp_2 <- append(intron_temp_1, intron_list[j*2])
          intron_temp_2 <- append(intron_temp_2, intron_list[j*2+1])
          intron[nrow(intron) + 1,] <- intron_temp_2 
          intron_temp_2 <- vector()
        }
        intron_list <- vector()
      }
    }
  }
  else if(exon[i,6] == -1){
    if(i == 1){
      if(exon[i,2] == exon[i+1,2]){
        intron_list <- append(intron_list, exon[i,1])
        intron_list <- append(intron_list, exon[i,2])
        intron_list <- append(intron_list, exon[i,3])
        intron_list <- append(intron_list, exon[i,6])
        intron_list <- append(intron_list, exon[i,4])
        intron_list <- append(intron_list, exon[i,5])
      }
    }
    else if(i == nrow(exon)){
      if(exon[i,2] == exon[i-1,2]){
        intron_list <- append(intron_list, exon[i,4], after = 4)
        intron_list <- append(intron_list, exon[i,5], after = 5)
        intron_temp_1 <- intron_list[1:4]
        intron_list <- intron_list[5:length(intron_list)]
        for(j in c(1:(length(intron_list)/2-1))){
          intron_temp_2 <- append(intron_temp_1, intron_list[j*2])
          intron_temp_2 <- append(intron_temp_2, intron_list[j*2+1])
          intron[nrow(intron) + 1,] <- intron_temp_2 
          intron_temp_2 <- vector()
        }
        intron_list <- vector()
      }
    }
    else{
      if(exon[i,2] != exon[i-1,2] & exon[i,2] == exon[i+1,2]){
        intron_list <- append(intron_list, exon[i,1])
        intron_list <- append(intron_list, exon[i,2])
        intron_list <- append(intron_list, exon[i,3])
        intron_list <- append(intron_list, exon[i,6])
        intron_list <- append(intron_list, exon[i,4])
        intron_list <- append(intron_list, exon[i,5])
      }
      else if(exon[i,2] == exon[i-1,2] & exon[i,2] == exon[i+1,2]){
        intron_list <- append(intron_list, exon[i,4], after = 4)
        intron_list <- append(intron_list, exon[i,5], after = 5)
      }
      else if(exon[i,2] == exon[i-1,2] & exon[i,2] != exon[i+1,2]){
        intron_list <- append(intron_list, exon[i,4], after = 4)
        intron_list <- append(intron_list, exon[i,5], after = 5)
        intron_temp_1 <- intron_list[1:4]
        intron_list <- intron_list[5:length(intron_list)]
        for(j in c(1:(length(intron_list)/2-1))){
          intron_temp_2 <- append(intron_temp_1, intron_list[j*2])
          intron_temp_2 <- append(intron_temp_2, intron_list[j*2+1])
          intron[nrow(intron) + 1,] <- intron_temp_2 
          intron_temp_2 <- vector()
        }
        intron_list <- vector()
      }
    }
  }
}
colnames(intron) <- intron[1,]
intron <- intron[c(2:nrow(intron)),]
intron$intron_chrom_start <- as.numeric(intron$intron_chrom_start)
intron$intron_chrom_end <- as.numeric(intron$intron_chrom_end)
intron$chromosome_name <- as.numeric(intron$chromosome_name)
intron$strand <- as.numeric(intron$strand)

# peaks in intron regions
intron_regions <- GRanges(seqnames = Rle(paste0("chr", intron$chromosome_name)),
                          ranges = IRanges(start = intron$intron_chrom_start, end = intron$intron_chrom_end),
                          srtand = Rle(rep("*", nrow(intron))),
                          gene = intron$external_gene_name)
for(name in read_list){
  intron_peak_sum <- sum(countOverlaps(intron_regions, eval(parse(text = name))))
  print(paste(name, "has", intron_peak_sum, "in intron regions", sep = " "))
}

# peaks in certain gene intron with normalization
peaks_in_gene_intron <- function(gene_name){
  index_list <- which(intron$external_gene_name == gene_name)
  assign(gene_name, GRanges(seqnames = Rle(paste0("chr", intron[index_list,3])),
                            ranges = IRanges(start = intron[index_list,5], end = intron[index_list,6]),
                            strand = Rle("*",length(index_list)), 
                            gene = intron[index_list,2]))
  for(name in read_list){
    index <- match(name, read_list)
    peak_sum <- sum(countOverlaps(eval(parse(text = gene_name)), eval(parse(text = name)))) / fold_list[index]
    print(paste(name, "has", peak_sum, "peaks in", gene_name, "intron regions", sep = " "))
  }
}
peaks_in_gene_intron("")

# 5' UTR isolation and setting
utr_5 <- getBM(attributes = c("ensembl_gene_id","external_gene_name",
                            "chromosome_name", "5_utr_start",
                            "5_utr_end", "strand"), mart = ds,
                            filter = "chromosome_name", 
                            values = strsplit(chromosome, split = "r")[[1]][2])
utr_5 <- na.omit(utr_5)
utr_5_regions <- GRanges(seqnames = Rle(paste0("chr", utr_5$chromosome_name)),
                         ranges = IRanges(start = utr_5$'5_utr_start', end = utr_5$'5_utr_end'),
                         strand = Rle(rep("*", nrow(utr_5))),
                         gene = utr_5$external_gene_name)
utr_5_total_length <- sum(width(reduce(utr_5_regions)))
pos_utr_5 <- utr_5[unique(queryHits(findOverlaps(utr_5_regions, enriched_regions))),]

# peaks in 5' UTR regions
for(name in read_list){
  utr_5_peak_sum <- sum(countOverlaps(utr_5_regions, eval(parse(text = name))))
  print(paste(name, "has", utr_5_peak_sum, "peaks in 5' UTR regions", sep = " "))
}

# 3' UTR isolation and setting
utr_3 <- getBM(attributes = c("ensembl_gene_id","external_gene_name",
                              "chromosome_name", "3_utr_start",
                              "3_utr_end", "strand"), mart = ds,
                              filter = "chromosome_name", 
                              values = strsplit(chromosome, split = "r")[[1]][2])
utr_3 <- na.omit(utr_3)
utr_3_regions <- GRanges(seqnames = Rle(paste0("chr", utr_3$chromosome_name)),
                         ranges = IRanges(start = utr_3$'3_utr_start', end = utr_3$'3_utr_end'),
                         strand = Rle(rep("*", nrow(utr_3))),
                         gene = utr_3$external_gene_name)
utr_3_total_length <- sum(width(reduce(utr_3_regions)))
pos_utr_3 <- utr_5[unique(queryHits(findOverlaps(utr_3_regions, enriched_regions))),]

# peaks in 3' UTR regions
for(name in read_list){
  utr_3_peak_sum <- sum(countOverlaps(utr_3_regions, eval(parse(text = name))))
  print(paste(name, "has", utr_3_peak_sum, "peaks in 3' UTR regions", sep = " "))
}

# intragenic regions isolation (gene distance > 200kb)
gd_list <- c(1)
for(i in c(1:nrow(egs))){
  if(i != nrow(egs)){
    if(egs[i+1,4] - egs[i,5] > 200000){
      gd_list <- append(gd_list, egs[i,4]-100000)
      gd_list <- append(gd_list, egs[i,5]+100000)
    }
  }
  else{
    if(egs[i,4] - egs[i-1,5] > 200000){
      gd_list <- append(gd_list, egs[i,4]-100000)
      gd_list <- append(gd_list, egs[i,5]+100000)
    }
  }
}
gd_list <- append(gd_list, seqlengths(genome)[6])
gd_temp <- data.frame("chromosome_name", "gd_start", "gd_end", stringsAsFactors=FALSE)
gd_temp_list <- c(chromosome)
for(i in c(1:(length(gd_list)/2))){
  gd_temp_list <- append(gd_temp_list, gd_list[i*2-1])
  gd_temp_list <- append(gd_temp_list, gd_list[i*2])
  gd_temp[nrow(gd_temp)+1,] <- gd_temp_list
  gd_temp_list <- c(strsplit(chromosome, split = "r")[[1]][2])
}
colnames(gd_temp) <- gd_temp[1,]
gd_temp <- gd_temp[c(2:nrow(gd_temp)),]
gd_temp$chromosome_name <- as.numeric(gd_temp$chromosome_name)
gd_temp$gd_start <- as.numeric(gd_temp$gd_start)
gd_temp$gd_end <- as.numeric(gd_temp$gd_end)
gd_regions <- GRanges(seqnames = Rle(gd_temp$chromosome_name), 
                      ranges = IRanges(start = gd_temp$gd_start, end = gd_temp$gd_end))

# peaks in intragenic regions (gene distance > 200kb)
for(name in read_list){
  intragene_peak_sum <- sum(countOverlaps(gd_regions, eval(parse(text = name))))
  print(paste(name, "has", intragene_peak_sum, "in intragenic regions", sep = " "))
}

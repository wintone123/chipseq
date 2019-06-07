chrom_track <- IdeogramTrack(chromosome = "chr1", genome = "mm10")
axis_track <- GenomeAxisTrack()
biomart_track <- BiomartGeneRegionTrack(chromosome = "chr1", genome = "mm10", name = "refseq",
                                        stacking = "squish", filter = list("with_refseq_mrna"=TRUE),
                                        collapseTranscripts = "longest", featureMap = fm, 
                                        biomart = mart)
cpg_track <- AnnotationTrack(cpg_chr1, chromosome = "chr1", genome = "mm10", name = "cpg",
                             stacking = "dense", col = "blue", size = 1)
rep1_track <- AnnotationTrack(rep1, chromosome = "chr1", genome = "mm10", name = "rep1",
                              stacking = "dense", col = "green", size = 2)
rep2_track <- AnnotationTrack(rep2, chromosome = "chr1", genome = "mm10", name = "rep1",
                              stacking = "dense", col = "green", size = 2)
overlap_track <- AnnotationTrack(overlap_chr1, chromosome = "chr1", genome = "mm10", name = "rep1",
                              stacking = "dense", col = "yellow", size = 2)
plotTracks(c(chrom_track, axis_track, biomart_track, cpg_track, rep1_track, rep2_track, overlap_track),
           from = 4240000, to = 5570000, transcriptAnnotation = "symbol")
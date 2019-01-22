# load peak file
rep1 <- import.bed("c:/chipseq/test3/rep1_best_chr6.bed")
rep2 <- import.bed("c:/chipseq/test3/rep2_best_chr6.bed")
control <- import.bed("c:/chipseq/test3/control_best_chr6.bed")
rep1_peak <- import_np("c:/chipseq/test3/rep1_best_peaks_chr6.narrowpeak")
rep2_peak <- import_np("c:/chipseq/test3/rep2_best_peaks_chr6.narrowpeak")


# pic
rep1_peak_track <- AnnotationTrack(rep1_peak, chromosome = "chr6", genome = "mm10", name = "rep1", col = "blue", fill = "blue", stacking = "dense")
rep2_peak_track <- AnnotationTrack(rep2_peak, chromosome = "chr6", genome = "mm10", name = "rep1", col = "blue", fill = "blue", stacking = "dense")
enriched_track <- AnnotationTrack(enriched, chromosome = "chr6", genome = "mm10", name = "rep1", col = "yellow", fill = "yellow", stacking = "dense")
axis_track <- GenomeAxisTrack()
genome_track <- IdeogramTrack(chromosome = "chr6", genome = "mm")
plotTracks(c(genome_track, axis_track, rep1_peak_track, rep2_peak_track, enriched_track),
           from = 122600000, to = 122700000)
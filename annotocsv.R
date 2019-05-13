library(BSgenome)
library(dplyr)

imput_path <- "/mnt/c/others/Mus_musculus.GRCm38.96.gtf"
out_path <- "/mnt/c/others/anno_mouse"
chr_list <- c(1:19, "X", "Y")

cat("loading file", "\n")
anno <- import.gff2(imput_path)
anno_df <- data.frame(chr = factor(anno@seqnames),
                      start = start(anno),
                      end = end(anno),
                      strand = anno@strand,
                      gene_name = anno$gene_name,
                      type = anno$type,
                      gene_id = anno$gene_id,
                      transcript_id = anno$transcript_id,
                      gene_biotype = anno$gene_biotype,
                      stringsAsFactors = FALSE)
anno_df_fil <- filter(anno_df, gene_biotype == "protein_coding" & type == "gene")

for (chrom in chr_list) {
    cat(paste0("processing chr", chrom), "\n")
    anno_chr <- filter(anno_df_fil, chr == chrom)
    write.csv(anno_chr, paste0(out_path, "/anno_chr" , chrom , ".csv"), row.names = FALSE)
}
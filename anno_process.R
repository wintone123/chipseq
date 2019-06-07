# library
library(GenomicRanges)

# functions
intron_maker <- function(input_fil) {  
  intron_df <- data.frame(start = NA, end = NA)
  cond <- TRUE
  s <- 1
  while (cond) {
    exon_list <- vector()
    for (i in s:length(input)) {
      if (i == length(input)) {
        cond <- FALSE
      }
      if (input$type[i] == "exon") {
        exon_list <- append(exon_list, c(start(input)[i], end(input)[i]))
      } else if (input$type[i] != "exon" & length(exon_list) != 0) {
        s <- i
        break()
      }
    }
    if (length(exon_list) > 2) {
      for (j in 1:(length(exon_list)/2-1)) {
        intron_df <- rbind(intron_df, data.frame(start = exon_list[j*2]+1, end = exon_list[j*2+1]-1))
      }
    }
    cat(nrow(intron_df),"\r")
  }
  return(intron_df)
}

# load file
path <- ""
gff2 <- ""
input <- import.gff3(file.path(path, gff2))

# filtering
input_fil <- input[input$type %in% c("")]
overlap_peaks <- function(GRanges_1, GRanges_2){
    df_1 <- data.frame(start(GRanges_1), end(GRanges_1))
    df_2 <- data.frame(start(GRanges_2), end(GRanges_2))
    start_out <- vector()
    end_out <- vector()
    loop_s <- 1
    for(i in c(1:nrow(df_1))){
        start_1 <- df_1[i,1]
        end_1 <- df_1[i,2]
        for(j in c(loop_s:nrow(df_2))){
            start_2 <- df_2[j,1]
            end_2 <- df_2[j,2]
            if (end_2 <= start_1){
                loop_s <- j
            } else if (start_2 >= end_1) {
                break
            } else{
                if(start_2 <= start_1){
                    start_out <- append(start_out, start_1)
                    if(end_2 <= end_1 & end_2 > start_1){
                        end_out <- append(end_out, end_2)
                    } else if(end_2 > end_1){
                        end_out <- append(end_out, end_1)
                    }
                } else if(start_2 > start_1 & start_2 < end_1){
                    start_out <- append(start_out, start_2)
                    if(end_2 <= end_1){
                        end_out <- append(end_out, end_2)
                    } else{
                        end_out <- append(end_out, end_1)
                    }
                }
            }
        }
    }
    new_GRanges <- GRanges(seqnames = Rle(GRanges_1@seqnames[1]),
                           ranges = IRanges(start = start_out, end = end_out),
                           strand = Rle("*"))
    return(new_GRanges)
}
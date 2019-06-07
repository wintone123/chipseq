# binning sequence as certain length
TileSequence <- function(seqname, start, end, tilewidth){
    start_list <- seq(start, end, by = tilewidth)
    end_list <- start_list + tilewidth -1
    if (start_list[length(start_list)] == end) {
        end_list <- end_list[c(1:length(end_list)-1)]
        end_list[length(end_list)] <- start_list[length(start_list)]
        start_list <- start_list[c(1:length(start_list)-1)]
    } else {
        end_list[length(end_list)] <- end
    }
    GRanges_temp <- GRanges(seqnames = Rle(seqname),
                          ranges = IRanges(start = start_list, end = end_list),
                          strand = Rle("*"))
    return(GRanges_temp)
}

# binning sequences in certain collection
TileCollection <- function(collection, tilewidth) {
    for (i in 1:nrow(collection)) {
        start <- collection$start_position[i]
        end <- collection$end_position[i]
        start_list <- seq(start, end, by = tilewidth)
        end_list <- start_list + tilewidth -1
        if (i == 1) {
            if (start_list[length(start_list)] == collection$end_position[i]) {
                end_list <- end_list[c(1:length(end_list)-1)]
                end_list[length(end_list)] <- start_list[length(start_list)]
                start_list <- start_list[c(1:length(start_list)-1)]
            } else {
                end_list[length(end_list)] <- end
            }
            data_df <- data.frame(seqname = rep(collection$chromosome_name[i], length(start_list)),
                                  start = start_list, 
                                  end = end_list,
                                  strand = rep(collection$strand[i], length(start_list)))
        } else {
            if (start_list[length(start_list)] == end) {
                end_list <- end_list[c(1:length(end_list)-1)]
                end_list[length(end_list)] <- start_list[length(start_list)]
                start_list <- start_list[c(1:length(start_list)-1)]
            } else {
                end_list[length(end_list)] <- end
            }
            data_df_temp <- data.frame(seqname = rep(collection$chromosome_name[i], length(start_list)),
                                       start = start_list, 
                                       end = end_list,
                                       strand = rep(collection$strand[i], length(start_list)))
            data_df <- rbind(data_df, data_df_temp)
        }
    }
    GRanges_temp <- GRanges(seqnames = Rle(data_df$seqname),
                            ranges = IRanges(start = data_df$start, end = data_df$end),
                            strand = Rle(data_df$strand))
    return(GRanges_temp)
}

# import narrowpeak file
ImportNp <- function(path) {
    temp <- read.delim2(path, quote = "/", head = FALSE)
    GRanges_temp <- GRanges(seqnames = Rle(temp[,1]),
                            ranges = IRanges(start = temp[,2], end = temp[,3]),
                            strand = Rle("*"),
                            name = as.character(temp[,4]),
                            score = as.numeric(temp[,5]),
                            -log10p = as.numeric(temp[,8]),
                            -log10q = as.numeric(temp[,9]))
    return(GRanges_temp)
}

OverlapPeaks <- function(GRanges_1, GRanges_2) {
    df_1 <- data.frame(start(GRanges_1), end(GRanges_1))
    df_2 <- data.frame(start(GRanges_2), end(GRanges_2))
    start_out <- vector()
    end_out <- vector()
    loop_s <- 1
    for (i in 1:nrow(df_1)) {
        start_1 <- df_1[i,1]
        end_1 <- df_1[i,2]
        for (j in loop_s:nrow(df_2)) {
            start_2 <- df_2[j,1]
            end_2 <- df_2[j,2]
            if (end_2 <= start_1) {
                loop_s <- j
            } else if (start_2 >= end_1) {
                break
            } else {
                if (start_2 <= start_1) {
                    start_out <- append(start_out, start_1)
                    if (end_2 <= end_1 & end_2 > start_1) {
                        end_out <- append(end_out, end_2)
                    } else if(end_2 > end_1) {
                        end_out <- append(end_out, end_1)
                    }
                } else if (start_2 > start_1 & start_2 < end_1) {
                    start_out <- append(start_out, start_2)
                    if (end_2 <= end_1) {
                        end_out <- append(end_out, end_2)
                    } else {
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

# binning sequences in certain collection
TileCollection <- function(collection, bin_num) {
    for (i in 1:nrow(collection)) {
        start <- collection$start_position[i]
        end <- collection$end_position[i]
        seq_list <- seq(start, end, length.out = bin_num + 1)
        start_list <- vector(length = bin_num)
        end_list <- vector(length =  bin_num)
        for (j in 1:(bin_num + 1)) {
            if (j == 1) {
                start_list[j] <- seq_list[j]
            } else if (j == bin_num + 1) {
                end_list[j-1] <- seq_list[j]
            } else {
                start_list[j] <- seq_list[j] + 1
                end_list[j-1] <- seq_list[j]
            }
        }
        if (i == 1) {
            data_df <- data.frame(seqname = rep(collection$chromosome_name[i], bin_num),
                                  start = start_list, 
                                  end = end_list,
                                  strand = rep(collection$strand[i], bin_num),
                                  name = rep(collection$external_gene_name[i], bin_num))
        } else {
            data_df_temp <- data.frame(seqname = rep(collection$chromosome_name[i], bin_num),
                                       start = start_list, 
                                       end = end_list,
                                       strand = rep(collection$strand[i], bin_num),
                                       name = rep(collection$external_gene_name[i], bin_num))
            data_df <- rbind(data_df, data_df_temp)
        }
    }
    GRanges_temp <- GRanges(seqnames = Rle(data_df$seqname),
                            ranges = IRanges(start = data_df$start, end = data_df$end),
                            strand = Rle(data_df$strand),
                            name = Rle(data_df$name))
    return(GRanges_temp)
}
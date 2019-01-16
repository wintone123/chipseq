# binning sequence as certain length
tilesequence <- function(seqname, start, end, tilewidth){
  start_list <- seq(start, end, by = tilewidth)
  end_list <- start_list + tilewidth -1
  if(start_list[length(start_list)] == end){
    end_list <- end_list[c(1:length(end_list)-1)]
    end_list[length(end_list)] <- start_list[length(start_list)]
    start_list <- start_list[c(1:length(start_list)-1)]
  } else{
    end_list[length(end_list)] <- end
  }
  GRanges_temp <- GRanges(seqnames = Rle(seqname),
                          ranges = IRanges(start = start_list, end = end_list),
                          strand = Rle("*"))
  return(GRanges_temp)
}

# import narrowpeak file
import_np <- function(path){
  temp <- read.delim2(path, quote = "/", head = FALSE)
  GRanges_temp <- GRanges(seqnames = Rle(temp[,1]),
                          ranges = IRanges(start = temp[,2], end = temp[,3]),
                          strand = Rle("*"),
                          name = as.character(temp[,4]),
                          score = as.numeric(temp[,5]))
  return(GRanges_temp)
}
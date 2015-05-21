# Clean workspace
rm(list = ls())

#Arguments passed by shell (arg1: path)
args <- commandArgs(TRUE)

split.Rd.line <- function(x, n.char.max){
  # x <- "This is a line that describes a, function, that might be too long. Or not, I really do not know if this is working"
  # x <- "%This is a line that describes a, function, that might be too long. Or not, I really do not know if this is working"
  # x <- "This is a line that describes a, function, that might be too long.# Or not, I really do not know if this is working"
  # x <- "%This is a line that describes a, function, that might be too long.# Or not, I really do not know if this is working"
  # x <- "This is a line that \"de- scribes\" a, function, that might be too long. Or \'not, I\' really do not know if this is working"
  # n.char.max <- 25
  
  # Simple cases:
  # Empty string
  if(nchar(x) <= n.char.max) return(x)
  # Only hash
  if(gsub(" ", "", gsub("#", "", x)) == ""){
    return(substring(x, 1, n.char.max))
  }
    
  res <- NULL
  
  comment.line <- grepl("^%", x)
  url.line     <- grepl("\\url", x)
  inline.comment <- FALSE
  
  while(nchar(x) > n.char.max){
    
    whitespace.loc  <- gregexpr(" ", x)[[1]]
    quote.loc       <- gregexpr("\"|\'", x)[[1]]
    
    if(whitespace.loc[1] == -1){
      split.idx <- n.char.max
    }else{
      split.idx <- whitespace.loc[max(which((n.char.max - whitespace.loc < 0) == FALSE))]
      if(split.idx == 1){split.idx = n.char.max}
    }
    
    if(length(quote.loc) > 1){
      if(split.idx >= quote.loc[1] & split.idx <= quote.loc[2]){
        split.idx <- quote.loc[1] - 1
        if(split.idx <= 1){
          split.idx <- quote.loc[2] + 1
        }
      }
    }
    
    if(inline.comment){
        inline.comment.loc <- -1
      }else{
        inline.comment.loc <- min(gregexpr("#", x)[[1]])
      }
    if(inline.comment.loc > 0 & inline.comment.loc <= split.idx & !comment.line & !url.line){
      split.idx <- inline.comment.loc
      x         <- gsub("#", "", x)
      inline.comment <- TRUE
    }
    res <- c(res, substring(x, 1, split.idx-1))
    leading.character <- if (comment.line) "%" else {if (inline.comment) "#" else ""}
    x <- paste(leading.character, substring(x, split.idx, nchar(x)), sep = "")
    if(x == "") break
  }
  return(c(res, x))
}

rd.files <- list.files(args[1], pattern = "*.Rd")

for( file in rd.files){
  # file <- rd.files[1]
  # file <- rd.files[95]
  in.file <- scan(paste0(args[1], "/", file), what = character(), sep = "\n")
  out.file <- NULL
  for( i in 1:length(in.file)){
    # i <- 9
    out.file <- c(paste0(args[1], "/", out.file), split.Rd.line(in.file[i], 79))
    cat(".")
  }
  write(out.file, file, sep = "\n")
  cat("\n", match(file, rd.files), "\n")
}

cat("\n")




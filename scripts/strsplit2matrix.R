strsplit2matrix <- function (x, split, fixed = FALSE, perl = FALSE) 
{
  if (is.factor(x)) {
    x <- as.character(labels(x))
  }
  if (is.list(x)) {
    x <- x[[1]]
  }
  x <- strsplit(x, split = split, fixed = fixed, perl = perl)
  n.splits <- sapply(x, length)
  if (length(unique(n.splits)) != 1) {
    max.length <- max(n.splits)
    min.length <- min(n.splits)
    x <- lapply(x, function(x, n) {
      c(x, rep(NA, n - length(x)))
    }, max.length)
  }
  x <- as.data.frame(do.call("rbind", x), stringsAsFactors = FALSE)
  return(x)
}

fishersMethod <- function(x, df=2*length(x[!is.na(x)])) {
	x <- x[!is.na(x)]
	pchisq(-2*sum(log(x)), df, lower.tail=FALSE)
}
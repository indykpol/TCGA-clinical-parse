fishersMethod <- function(x) {
	x <- x[!is.na(x)]
	pchisq(-2*sum(log(x)), df=length(x), lower.tail=FALSE)
}
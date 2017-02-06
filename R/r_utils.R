# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:


# create date code of the form YYMMDD
# inputs:
#    d (optional): date code in the format of Sys.Date() for which to generate the date code, defaulting to 'today'
date_code <- function(d=NA){
	# reference http://www.r-cookbook.com/node/17
	if(is.na(d)) d <- Sys.Date()
	pattern <- '20([[:digit:]]{2})-([[:digit:]]{2})-([[:digit:]]{2})'
	paste(
		sub(pattern, '\\1', d), sub(pattern, '\\2', d), sub(pattern, '\\3', d),
		sep="")
}


#' @importFrom magrittr %>%
NULL


#' Returns TRUE if no values are duplicated
#'
#' This function tests if values in the sepcified columns have no
#' duplicate values. This the inverse of
#' \code{\link[base]{duplicated}}. This is a convenience function
#' meant to be used as a predicate in an \code{\link{assertr}}
#' verify statment.
#'
#' Warning: Since this uses \code{\link[base]{duplicated}}, the columns
#' are pasted together with '\\t', this will fail if any of the
#' specified columns contain '\\t' characters.
#'
#' @param ... columns from the data.frame
#' @return A vector of the same length that is TRUE when the values in
#' the specified columns are distinct
#' @seealso \code{\link{duplicated}}
#' @examples
#' z <- data.frame(a=c(1,2,3,4,5,6), b=c(1,1,1,3,3,3), c=c(1,2,3,1,2,3))
#' zz <- z %>% verify(no_duplicates(a))
#' zz <- z %>% verify(no_duplicates(b))    # verification failed! (4 failures)
#' zz <- z %>% verify(no_duplicates(c))    # verification failed! (3 failures)
#' zz <- z %>% verify(no_duplicates(a,b))
#' zz <- z %>% verify(no_duplicates(a,c))
#' zz <- z %>% verify(no_duplicates(b,c))
#' zz <- z %>% verify(no_duplicates(a,b,c))
#'
#' @export
no_duplicates <- function(...){
    args <- list(...)
    args$incomparables=F
    !duplicated.data.frame(args)
}

get_random_intervs <- function(len, n=1){
    as.data.frame(tibble(chrom=gintervals.all()$chrom[1], start=sample(1:(gintervals.all()$end[1] - len), n), end=start+len))
}

get_random_seq <- function(len, n=1, upper_case=TRUE){	
	# For some reason misha doesn't recognize GROOT and ALLGENOME
	.gcheckroot()
	GROOT <- .GlobalEnv$GROOT	
	ALLGENOME <- .GlobalEnv$ALLGENOME
	
    s <- gseq.extract(get_random_intervs(len=len, n=n))
    if (upper_case){
    	s <- toupper(s)
    }
    return(s)
}

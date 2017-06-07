get_random_intervs <- function(len, n=1){
    tibble(chrom=gintervals.all()$chrom[1], start=sample(1:(gintervals.all()$end[1] - len), n), end=start+len)
}

get_random_seq <- function(len, n=1){
    gseq.extract(get_random_intervs(len, n=n))
}

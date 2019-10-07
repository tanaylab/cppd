#' Filter for putative enhancers that are within gene tads
#' 
#' @param genes character vector of gene symbols
#' @param intervals putative enhancer intervals
#' @param tads tad intervals: intervals set with an extra id column
#' 
#' @return intervals from \code{intervals} with an extra 'id' column for TAD id and 'gene' with the genes from \code{genes} within the TAD
#' 
#' @export
cppd.gene_enhancers <- function(genes, intervals, tads){
	gene_intervals <- gintervals.load('intervs.global.tss') %>% 
		filter(geneSymbol %in% genes) %>% 
		select(chrom, start, end, strand, gene=geneSymbol) %>%
		as_tibble()
	gene_intervals <- gene_intervals %>% 
		gintervals.neighbors1(tads) %>% 
		filter(dist == 0) %>% 
		select(chrom, start, end, strand, gene, id)
	intervals <- intervals %>% 
		gintervals.neighbors1(tads) %>% 
		filter(dist == 0) %>% 
		select(-(chrom1:end1), -dist)	    
	gene_intervals <- gene_intervals %>% distinct(id, gene) %>% group_by(id) %>% nest(gene = gene)
	intervals <- intervals %>%
		filter(id %in% gene_intervals$id) %>% 
		left_join(gene_intervals, by='id') %>% 
		filter(!is.null(gene)) %>% 
		distinct(chrom, start, end, .keep_all=TRUE)
	
	return(intervals)
}

#' Get promoter region of genes
#' 
#' @param genes character vector of gene symbols
#' @param upstream number of base pairs upstream to the gene
#' @param downstream number of base pairs downstream to the gene
#' 
#' @return intervals set with the gene promoters
#' 
#' @export
cppd.gene_promoters <- function(genes, upstream=500, downstream=50){
    gintervals.load('intervs.global.tss') %>%
    	filter(geneSymbol %in% genes) %>% 
    	mutate(start = ifelse(strand == 1, start - upstream, start - downstream), 
    		   end = ifelse(strand == 1, end + downstream, end + upstream)) %>% 
    	as_tibble()
}


#' Get CpGs within a region
#' 
#' @param intervals intervals set
#' 
#' @return intervals set of CpGs within the given intervals
cppd.intervals_cpgs <- function(intervals){
    gintervals.neighbors1(giterator.intervals(iterator="intervs.global.seq_CG", intervals=intervals), intervals) 
}
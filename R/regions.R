#' @export
cppd.gene_enhancers <- function(genes, intervals, tads){
	gene_intervals <- gintervals.load('intervs.global.tss') %>% 
		filter(geneSymbol %in% genes) %>% 
		select(chrom, start, end, strand, gene=geneSymbol) %>%
		as.tibble()
	gene_intervals <- gene_intervals %>% 
		gintervals.neighbors1(tads) %>% 
		filter(dist == 0) %>% 
		select(chrom, start, end, strand, gene, id)
	intervals <- intervals %>% 
		gintervals.neighbors1(tads) %>% 
		filter(dist == 0) %>% 
		select(-(chrom1:end1), -dist)	
	gene_intervals <- gene_intervals %>% distinct(id, gene) %>% group_by(id) %>% nest(gene, .key='gene')
	intervals <- intervals %>%
		filter(id %in% gene_intervals$id) %>% 
		left_join(gene_intervals, by='id') %>% 
		filter(!is.null(gene)) %>% 
		distinct(chrom, start, end, .keep_all=TRUE)
	
	return(intervals)
}

#' @export
cppd.gene_promoters <- function(genes, upstream=500, downstream=50){
    gintervals.load('intervs.global.tss') %>%
    	filter(geneSymbol %in% genes) %>% 
    	mutate(start = ifelse(strand == 1, start - upstream, start - downstream), 
    		   end = ifelse(strand == 1, end + downstream, end + upstream)) %>% 
    	tbl_df()
}

choose_probes <- function(cands, regions, exp_regions, TM_range, probes_ofn=NULL, all_probes_ofn=NULL, regs_annots=NULL, n_probes=NULL, kmer_len=15, downsample=FALSE, rm_revcomp=TRUE){
    
    cands <- annotate_cands(cands, exp_regions)

    loginfo('counting CpGs per strand')
    cands <- get_min_cgs_per_strand(cands)  

    loginfo('sorting candidates')
    cands <- sort_cands(cands, TM_range)
    
    loginfo('Analyzing single CpG candidates')
    # single CpG
    probes00 <- get_probes(cands, 0, 0, TM_range)
    probes11 <- get_probes(cands, 1, 1, TM_range)
    probes10 <- get_probes(cands, 1, 0, TM_range)
    probes01 <- get_probes(cands, 0, 1, TM_range)

    # two CpGs
    loginfo('Analyzing 2 CpGs candidates')
    # Try to get 2 probes from the 0 strand
    probes20 <- cands %>% filter(strand == -1) %>% get_probes(2, 0, TM_range, 2)
    probes02 <- cands %>% filter(strand == 1) %>% get_probes(0, 2, TM_range, 2)

    # take the important 2-0/0-2 regions and discard the rest
    left_regs <- probes20 %>% group_by(chrom, start_reg, end_reg) %>% filter(keep, n() == 1) %>% ungroup %>% distinct(chrom, start_reg, end_reg)
    probes20 <- probes20 %>% group_by(chrom, start_reg, end_reg) %>% filter(n() == 2) %>% ungroup
    probes20_cgs <- cands %>% inner_join(left_regs) %>% ungroup %>% get_probes(2, 0, TM_range, 1)


    left_regs <- probes02 %>% group_by(chrom, start_reg, end_reg) %>% filter(keep, n() == 1) %>% ungroup %>% distinct(chrom, start_reg, end_reg)
    probes02 <- probes02 %>% group_by(chrom, start_reg, end_reg) %>% filter(n() == 2) %>% ungroup
    probes02_cgs <- cands %>% inner_join(left_regs) %>% ungroup %>% get_probes(0, 2, TM_range, 1)

    # high number of CpGs in important regions
    high_cg_cands <- cands %>% group_by(chrom, start_reg, end_reg) %>% filter(min(cg_num) >= 2, keep) %>% ungroup

    probes <- bind_rows(probes00, probes11, probes10, probes01, probes20, probes02, probes20_cgs, probes02_cgs) 

    loginfo('adding sequence')    
    probes <- add_probe_seq(probes)

    if (rm_revcomp){
        loginfo('checking reverse complementarity')   
        p1 <- probes %>% arrange(-keep, chrom, start_reg, end_reg) %>% add_kmer_revcomp(k=kmer_len)  
        probes <- p1 %>% group_by(chrom, start_reg, end_reg) %>% filter(!any(revcomp_dup)) %>% ungroup %>% select(-revcomp_dup)
    } else {
        logging::logwarn('rm_revcomp is FALSE: not checking for reverse complementarity')
        p1 <- probes        
    }
    

    loginfo('number of probes: %d', nrow(probes))

    if (downsample && !is.null(n_probes) && n_probes > nrow(probes)){
        probes <- downsample_probes(probes, n_probes, all_probes_ofn=all_probes_ofn)
    }
    
    if (!is.null(regs_annots)){
        loginfo('adding annotations')
        regs_annot <- read_csv(regs_annots) %>% select(chrom, start, end, type)
        probes_annot <- probes %>% select(chrom, start, end, strand) %>% gintervals.neighbors1(regs_annot %>% select(chrom, start, end, type), maxneighbors=10, maxdist=300) %>% group_by(chrom, start, end, strand) %>% summarise(type = unique(type) %>% paste(collapse=','))
        probes <- probes %>% left_join(probes_annot %>% rename(annot_type=type))
    }

    # add original region ids
    probes <- regions %>% 
        gintervals.neighbors1(select(probes, chrom, start=start_reg, end=end_reg)) %>%
        filter(dist == 0) %>% 
        group_by(chrom1, start1, end1) %>% 
        summarise(reg_id = paste(id, collapse='_')) %>% 
        ungroup %>% 
        select(chrom=chrom1, start_reg=start1, end_reg=end1, reg_id) %>% 
        right_join(probes) %>% 
        select(chrom, start, end, strand, keep, start_reg, end_reg, reg_id, everything())
    

    if (!is.null(probes_ofn)){
        write_csv(probes, probes_ofn)
    }

    return(probes)   
}

##########################################################################

annotate_cands <- function(cands, exp_regs){        
    loginfo('counting CpGs')
    cands <- add_cg_num(cands)
    loginfo('adding CG content')
    cands <- add_GC_cont(cands) 

    # add keep annotation     
    cands <- cands %>% left_join(exp_regs %>% select(chrom, start_reg, end_reg, keep))
    return(cands)
}

add_cg_num <- function(cands){ 
    cands <- cands %>% mutate(cg_num = str_count(seq, 'CG'))
    return(cands)
}

add_GC_cont <- function(cands){
    cands <- cands %>% mutate(GC = str_count(seq, 'G') + str_count(seq, 'C'))
    return(cands)

}

get_min_cgs_per_strand <- function(cands){
    strands_min_cgs <- cands %>% group_by(chrom, start_reg, end_reg, strand, keep) %>% summarise(min_cgs_strand = min(cg_num)) %>% ungroup %>% mutate(strand = ifelse(strand == 1, 'min_cgs_plus', 'min_cgs_minus')) %>% ungroup %>% spread(strand, min_cgs_strand)
    
    cands <- cands %>% left_join(strands_min_cgs %>% select(-keep), by=c('chrom', 'start_reg', 'end_reg'))
    return(cands)
}

sort_cands <- function(cands, TM_range){
    cands %>% mutate( 
            max_TM_diff = ifelse(between(max_TM, TM_range[1], TM_range[2]), 0, pmin(abs(max_TM - TM_range[1]), abs(max_TM - TM_range[2]))),  
            min_TM_diff = ifelse(between(min_TM, TM_range[1], TM_range[2]), 0, pmin(abs(min_TM - TM_range[1]), abs(min_TM - TM_range[2]))), 
            best_TM_diff = pmin(max_TM_diff, min_TM_diff)) %>% 
        select(-max_TM_diff, -min_TM_diff) %>%
        arrange(chrom, start_reg, end_reg, strand, cg_num, best_TM_diff, kmers, homopol, revcomp, strand) 
}

get_probes <- function(sorted_cands, min_p_cgs, min_m_cgs, TM_range, num_probs=1){        
    sorted_cands <- sorted_cands %>% filter(cg_num <= max(min_p_cgs, min_m_cgs), min_cgs_plus == min_p_cgs, min_cgs_minus == min_m_cgs) %>%
        get_best_probe(num_probs) %>%   
        mutate(type = paste0(min_p_cgs, '_', min_m_cgs))
}

get_best_probe <- function(sorted_cands, num_probs=1){  
    if (num_probs == 1){        
        return(sorted_cands %>%
            group_by(chrom, start_reg, end_reg, strand) %>% 
            slice(1) %>%
            ungroup)
    }

    if (num_probs == 2){        
        return(sorted_cands %>% 
            group_by(chrom, start_reg, end_reg, strand) %>% 
            mutate(s = start[1], e = end[1]) %>% 
            filter((start == s & end == e) | !((start > s & start < e) |  (end > s & end < e))) %>% 
            slice(1:2) %>% 
            select(-s, -e) %>%
             ungroup)
    }    
}

add_probe_seq <- function(probes){    
    probes %>% select(-id, -seq) %>% mutate(genomic_seq = toupper(gseq.extract(.)), seq = gsub('CG', 'AG', genomic_seq), seq = gsub('C', 'T', seq))    
}



downsample_probes <- function(probes, n_probes, all_probes_ofn=NULL){
    loginfo(sprintf('downsampling %s probes to %s', nrow(probes), n_probes))
    if (!is.null(all_probes_ofn)){
        write_csv(probes, all_probes_ofn)
    }
    probes <- probes %>% group_by(chrom, start_reg, end_reg) %>% mutate(n_reg = n()) %>% ungroup
    probes_ds <- probes %>% arrange(-keep, chrom, start_reg, end_reg) %>% slice(1:n_probes)

    # make sure that we have all the probes from each region
    probes_ds <- probes_ds %>% left_join(probes %>% select(chrom, start_reg, end_reg, n_reg_org=n_reg)) %>% filter(n_reg == n_reg_org) %>% select(-n_reg_org, -n_reg)

    return(probes_ds)
}

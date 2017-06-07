
get_candidates <- function(regs, probe_len, step, expand, threads, bowtie_bin, bissli2_idx, max_mapping_k, jellyfish_unmethylated_db, jellyfish_methylated_db, jellyfish_db_type, kmer_k, jellyfish_bin, TM_range, filter_cands=TRUE, max_map, max_homopol, max_revcomp, max_kmers, ...){

    if (threads > 1){
        doMC::registerDoMC(threads)
    }
    
    regs_exp <- do.call_ellipsis(regions2seq, list(sprobes_all=regs, expand=expand), ...) %>% mutate(start_reg=start, end_reg=end)
    # regions2seq(regs, expand=expand, ...) %>% mutate(start_reg=start, end_reg=end)

    max_str_len <- max(str_length(regs_exp$seq) )
    cands <- map_df(seq(1, max_str_len-probe_len, step), ~ get_cand_seq(regs_exp, ., probe_len)) 

    init_cands_num <- nrow(cands)
    loginfo('initial number of candidates: %d', init_cands_num)

    # add the converted sequences
    cands <- cands %>% mutate(seqc = convert_seq(seq, methylated=FALSE), seqcm = convert_seq(seq, methylated=TRUE))    

    loginfo('calculating TM')
    cands <- cands %>% mutate(TM_c = calc_TM(seqc), TM_cm = calc_TM(seqcm)) %>% 
            mutate(min_TM = pmin(TM_c, TM_cm), max_TM = pmax(TM_c, TM_cm)) %>% 
            select(-TM_c, -TM_cm)
    if (filter_cands){
        cands <- cands %>% filter(min_TM >= TM_range[1], max_TM <= TM_range[2])
        loginfo('candidates left: %s (%f)', nrow(cands), nrow(cands) / init_cands_num)
    }

    loginfo('calculating homopolymer')
    cands <- cands %>% mutate(homopol_c = count_max_homopolymer(seqc), homopol_cm = count_max_homopolymer(seqcm))
    if (filter_cands){
        cands <- cands %>% merge_m_um_stats('homopol', pmax) %>% filter(homopol <= max_homopol)
        loginfo('candidates left: %d (%f)', nrow(cands), nrow(cands) / init_cands_num)
    }

    loginfo('calculating count_max_revcomp')
    cands <- cands %>% mutate(revcomp_c = count_max_revcomp(seqc), revcomp_cm = count_max_revcomp(seqcm))
    if (filter_cands){
        cands <- cands %>% merge_m_um_stats('revcomp', pmax) %>% filter(revcomp <= max_revcomp)
        loginfo('candidates left: %d (%f)', nrow(cands), nrow(cands) / init_cands_num)
    }

    loginfo('calculating mapping')
    cands <- cands %>% mutate(
        mapping_c = count_mappings(seqc, threads=threads, bowtie_bin=bowtie_bin, bissli2_idx=bissli2_idx, max_k=max_mapping_k), 
        mapping_cm = count_mappings(seqcm, threads=threads, bowtie_bin=bowtie_bin, bissli2_idx=bissli2_idx, max_k=max_mapping_k))
    if (filter_cands){
        cands <- cands %>% merge_m_um_stats('mapping', pmax) %>% filter(mapping <= max_map)
        loginfo('candidates left: %d (%f)', nrow(cands), nrow(cands) / init_cands_num)
    }

    loginfo('calculating kmers')
    cands <- cands %>% mutate(
        kmers_c = count_genome_kmers(seqc, jellyfish_db=jellyfish_unmethylated_db, k=kmer_k, jellyfish_bin=jellyfish_bin, jellyfish_db_type=jellyfish_db_type, threads=threads), 
        kmers_cm = count_genome_kmers(seqcm, jellyfish_db=jellyfish_methylated_db, k=kmer_k, jellyfish_bin=jellyfish_bin, jellyfish_db_type=jellyfish_db_type, threads=threads))

    if (filter_cands){
        cands <- cands %>% merge_m_um_stats('kmers', pmax) %>% filter(kmers <= max_kmers)
        loginfo('candidates left: %d (%f)', nrow(cands), nrow(cands) / init_cands_num)
    }

    return(list(cands=cands, regs_exp=regs_exp))   

}


############################################################################
get_cand_seq <- function(regs_exp, offset, probe_len){
    regs_exp %>% mutate(seq = substring(seq, offset, offset + probe_len - 1), start = start + offset - 1, end = start + probe_len) %>% filter(str_length(seq) == probe_len) %>% select(-cgs)
}

convert_seq <- function(seqs, methylated=FALSE){
    seqs <- toupper(seqs)
    if (methylated) {
        return(gsub('C(?!G)', 'T', seqs, perl=T))
    } else {
        return(gsub('C', 'T', seqs))
    }
}

calc_TM <- function(seqs){
    data_frame(seq=toupper(seqs)) %>% 
        mutate(
            'A' = str_count(seq, 'A'), 
            'T' = str_count(seq, 'T'), 
            'G' = str_count(seq, 'G'), 
            'C' = str_count(seq, 'C'), 
            len=str_length(seq), 
            TM = ifelse(
                len < 14, 
                (A + T) * 2 + (G + C) * 4, 
                64.9 + 41 * (G + C- 16.4) / (A + T + G + C))
            ) %>%
        .$TM      
}

count_max_homopolymer <- function(seqs){
    str_locate_all(seqs, '(C+|A+|T+|G+)') %>%  map_dbl(function(x) max(x[, 2] - x[, 1]) + 1)
}

count_max_revcomp <- function(seqs) {
    seqs1 <- strsplit(seqs, '')
    seqs2 <- strsplit(gseq.rev_comp(seqs), '')
    map2_dbl(seqs1, seqs2, function(s1, s2) {min(which(s1 != s2))})    
}

write_fasta <- function(seqs, fn=tempfile(fileext = '.fa')){
    lines <- map2(seqs, paste0('>read', 1:length(seqs)), function(x, y) c(y, x)) %>% do.call('c', .) %>% tibble(seq=.)
    fwrite(lines, fn, col.names=F)
    return(fn)
}

count_mappings <- function(seqs, threads=10, bowtie_bin='/net/mraid14/export/data/users/eladch/tools/CO6/bowtie2/2.2.6/bin/bowtie2-align-l', bissli2_idx='/net/mraid14/export/data/tools/bissli2/hg19/hg19', max_k=5){
    fasta_fn <- write_fasta(seqs)
    on.exit(qq('rm -f @{fasta_fn}'))

    cmd <- qq('@{bowtie_bin} --wrapper basic-0 -f -U @{fasta_fn} -x @{bissli2_idx}  --quiet -k @{max_k} -p @{threads} --reorder | samtools view -F 4 | awk \'{print $1}\'')
    mappings <- fread(cmd, col.names='read') %>% count(read)

    n_mapping <- tibble(seq=seqs, read=paste0('read', 1:length(seqs))) %>% left_join(mappings, by='read') %>% mutate(n = ifelse(is.na(n), 0, n)) %>% .$n

    return(n_mapping)   
    
}

merge_m_um_stats <- function(cands, fields, func){
    for (field in fields){
        c_f <- qq('@{field}_c')
        cm_f <- qq('@{field}_cm')
        cands[field] <- func(cands[[c_f]], cands[[cm_f]]) 
        cands <- cands %>% select(-one_of(c(c_f, cm_f)))   
    }    
    return(cands)    
}

############################################################################
split_intervals <- function(intervs){
  intervs <- intervs %>% mutate(l = end - start)
  intervs %>% mutate(end = end - ceiling(l/2)) %>% bind_rows(intervs %>% mutate(start = start + floor(l/2))) %>% select(-l)
}

regions2seq <- function(sprobes_all, expand, max_len=651, exp_tab_fn=NULL, min_cgs=1){        
    sprobes <- sprobes_all %>% mutate(start = start - expand, end = end + expand) %>% as.data.frame %>% gintervals.force_range() %>% gintervals.canonic() %>% mutate(l = end - start, split = l > as.numeric(max_len)) %>% tbl_df
    
    loginfo("max_len: %s", max_len)
    loginfo("expand: %s", expand)
    sprobes <- sprobes %>% filter(!split) %>% bind_rows(split_intervals(sprobes %>% filter(split))) %>% select(-split, -l) %>% mutate(new_id = 1:n())

    exp_tab <- sprobes_all %>% as.data.frame %>% gintervals.neighbors1(sprobes) %>% select(-dist)
    if (!is.null(exp_tab_fn)){        
        exp_tab %>% mutate_if(is.double, as.integer) %>% write_tsv(exp_tab_fn)    
    }    

    sprobes <- exp_tab %>% group_by(chrom=chrom1, start=start1, end=end1, new_id) %>% summarise(keep=any(keep)) %>% ungroup %>% select(chrom, start, end, id=new_id, keep)

    loginfo("extracting sequences\n")    
    sprobes <- sprobes %>% mutate(strand = 1) %>% select(chrom, start, end, strand, id, keep) 
    sprobes$seq <- gseq.extract(sprobes) %>% toupper 

    loginfo("counting CpGs\n")
    sprobes <- sprobes %>% mutate(cgs = str_count(seq, 'CG'))   

    loginfo("adding reverse strand\n")
    sprobes_minus <- sprobes %>% mutate(strand = -1) %>% select(chrom, start, end, strand, id, keep, cgs)
    sprobes_minus$seq <- gseq.extract(sprobes_minus) %>% toupper    

    sprobes <- bind_rows(sprobes, sprobes_minus)    

    return(sprobes)
}
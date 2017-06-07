create_genome_kmer_db <- function(fasta_files, db_ofn, jellyfish_bin, jellyfish_params='', k=14){
    system(qq('@{jellyfish_bin} count -o @{db_ofn} @{jellyfish_params} @{paste(fasta_files, collapse=" ")}'))
}

seqs2kmers <- function(seqs_df, k, group_columns='id'){    
    seq_len <- str_length(seqs_df$seq[1])
    seq_mat <- seqs_df %>% extract(seq, into=paste0('c', 1:seq_len), paste(rep('(.)', seq_len), collapse=''))
    kmers <- tibble(s = 1:(seq_len - k + 1), e = s + k - 1)  %>%
        plyr::adply(1, function(x) 
            seq_mat %>% 
                select(one_of(c(group_columns, paste0('c', x$s:x$e) ))) %>% 
                unite('kmer', one_of( paste0('c', x$s:x$e) ), sep='')) %>% 
        select(-s, -e) %>% tbl_df
    return(kmers)
}

count_genome_kmers <- function(seqs, jellyfish_db, k, jellyfish_db_type, jellyfish_bin=NULL, threads=1){
    if (jellyfish_db_type == 'tsv'){
        counts <- count_genome_kmers_tsv(seqs, jellyfish_db, k, threads)
    } else {
        counts <- count_genome_kmers_jf(seqs, jellyfish_db, k, jellyfish_bin, threads)
    }
    return(counts)
}

##############################################################################
get_n_kmers_jf <- function(df, k, jellyfish_db, jellyfish_bin){
    kmers <- seqs2kmers(df, k)
    kmer_fn <- write_fasta(kmers$kmer)
    on.exit(qq('rm -f @{kmer_fn}'))
    kmer_counts <- fread(qq('@{jellyfish_bin} query -s @{kmer_fn} @{jellyfish_db}'), col.names=c('kmer', 'n_kmer'))
    kmers %>% left_join(kmer_counts, by='kmer') %>% group_by(id) %>% summarise(n_kmer = max(n_kmer)) %>% right_join(df, by='id')    
}

count_genome_kmers_jf <- function(seqs, jellyfish_db, k, jellyfish_bin, threads){
    seqs_df <- tibble(id = 1:length(seqs), seq=seqs)    
    kmers <- seqs_df %>% mutate(chunk = ntile(id, threads)) %>% ddply(.(chunk), function(x) get_n_kmers_jf(x, k, jellyfish_db, jellyfish_bin), .parallel=TRUE)    
    n_kmers <- kmers %>% arrange(id) %>% select(n_kmer) %>% .$n_kmer
    return(n_kmers)   
}

##############################################################################
get_n_kmers_tsv <- function(df, k, db){    
    kmers <- seqs2kmers(df, k)     
    kmer_counts <- kmers %>% left_join(db, by='kmer')
    kmer_counts %>% group_by(id) %>% summarise(n_kmer = max(n_kmer, na.rm=TRUE))      
}

count_genome_kmers_tsv <- function(seqs, jellyfish_db, k, threads){
    loginfo('reading db %s', jellyfish_db)
    db <- fread(jellyfish_db, col.names=c('kmer', 'n_kmer')) %>% tbl_df    
    seqs_df <- tibble(id = 1:length(seqs), seq=seqs)

    loginfo('getting kmers')      
    kmers <- seqs_df %>% mutate(chunk = ntile(id, threads)) %>% ddply(.(chunk), function(x) get_n_kmers_tsv(x, k, db), .parallel=TRUE)
    
    loginfo('counting')  
    n_kmers <- kmers %>% arrange(id) %>% select(n_kmer) %>% .$n_kmer    
    return(n_kmers)   
}

##############################################################################
add_kmer_revcomp <- function(df, k=15){
    # record the order of probes
    ids <- df %>% distinct(chrom, start, end) %>% mutate(id = 1:n())

    loginfo('getting kmers')        
    df_kmers <- df %>% select(chrom, start, end, strand, start_reg, end_reg, seq) %>% seqs2kmers(k, group_columns=c('chrom', 'start', 'end', 'strand', 'start_reg', 'end_reg'))
    
    loginfo('reverse complementing')
    revcomp_kmers <-  df_kmers %>% mutate(kmer = gseq.rev_comp(kmer)) %>% left_join(ids) %>% select(kmer, kmer_id=id)
    

    loginfo('looking for matches')
    probes2rm <- df_kmers %>% 
        left_join(ids) %>% 
        # find kmers that are the same as any reverse complement kmer
        left_join(revcomp_kmers, by='kmer') %>%         
        # remove self matches and probes with no match
        filter(id != kmer_id) %>% 
        filter(!is.na(kmer_id)) %>% 
        # remove only probes that were after the matching probe in the initial data frame
        filter(id > kmer_id) %>% 
        select(id = kmer_id) %>% 
        left_join(ids) %>% 
        distinct(chrom, start, end)
    
    # mark the duplicate prboes with revcomp_dup == TRUE
    df <- df %>% left_join(probes2rm %>% mutate(revcomp_dup = TRUE)) %>% mutate(revcomp_dup = !is.na(revcomp_dup))
    return(df)
}

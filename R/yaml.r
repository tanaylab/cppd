
read_yaml <- function(yaml_fn, eval_expr=TRUE, handlers=NULL){
    suppressPackageStartupMessages(library('stringr'))
    yaml <- yaml::yaml.load_file(yaml_fn, handlers = list(expr = function(x) paste0("!expr ", as.character(x))))    
    yaml <- .apply_handlers(yaml, handlers)
    yaml <- expand_yaml(yaml, eval_expr=eval_expr)
    return(yaml)
}

.apply_handlers <- function(yaml, handlers){
    if (!is.null(handlers)){
        for (name in names(handlers)){
            if (!is.null(yaml[[name]])){
                yaml[[name]] <- handlers[[name]](yaml[[name]])
            }            
        }
    }
    yaml
}

expand_yaml <- function(yaml, eval_expr=TRUE){
    suppressPackageStartupMessages(library('stringr'))
    incl_yaml <- list()
    if ('include' %in% names(yaml)){        
        for (incl_fn in yaml$include){
            incl_yaml <- plyr::defaults(incl_yaml, read_yaml(incl_fn, eval_expr=eval_expr))    
        }        
    }    
    e <- new.env()    
    e$yaml <- plyr::defaults(incl_yaml, yaml)
    yaml <- rapply(e$yaml, function(x) .interp(x, 'yaml', e), how='replace')        
    return(yaml)
}

.interp <- function(s, obj, envir, collapse=FALSE, eval_expr=TRUE){
    if (is.character(s)){
        res <- map(s, function(.x){       
                r <- GetoptLong::qq(gsub('\\$\\{', paste0('@\\{', obj, '$'), .x), envir=envir, collapse=collapse)    
                if (str_detect(s, '!expr ') && eval_expr){            
                    r <- gsub('!expr ', '', r) %>% parse(text=.) %>% eval
                }
                return(r)    
            })                
        
        if (length(res) == 1){
            res <- res[[1]]
        } else if (all(map_lgl(res, is.character))){
            res <- map_chr(res, ~ .)
        }

    } else {
        res <- s
    }
    return(res)    
}





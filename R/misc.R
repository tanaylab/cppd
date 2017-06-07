do.call_ellipsis <- function(f, additional_params=list(), ...){
    f_args <- names(as.list(args(f)))
    elipsis <- list(...)        
    if (!is.null(names(elipsis))){
        new_elipsis <- list()        
        for (x in names(elipsis)[names(elipsis) %in% f_args]) { 
            new_elipsis[[x]] <- elipsis[[x]] 
        }                
        do.call(f, c(additional_params, new_elipsis))
    } else {
        do.call(f, additional_params)
    }
}
# library imports ------------------------------------------------
#' @import tidyr
#' @import dplyr
#' @import purrr
#' @import ggplot2
#' @import misha
#' @importFrom data.table fwrite
#' @importFrom GetoptLong qq
#' @importFrom logging loginfo


.onLoad <- function(libname, pkgname) {        
    logging::basicConfig()
}



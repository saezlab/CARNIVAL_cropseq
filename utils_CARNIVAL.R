library(dplyr)
library(purrr)

#!! NB THE CODE BELOW WAS COPIED FROM VIGNETTE FOR DOROTHEA (written by Alberto Valdeolivas)
#see" https://github.com/alberto-valdeolivas/DoRothEA_Vignette

#' df2regulon: Function to group DoRothEA regulons
#'
#' This function takes a data frame containing the TF-target interactions
#' from DoRothEA and returns its associated regulons. 
#'
#' @import viper
#' @import dplyr
#' @import purrr
#' @param df A data frame containing the TF-target interactions from DoRothEA,
#' as stored in https://github.com/saezlab/ConservedFootprints/tree/master/data
#'
#' @return Object of class regulon. Check  
#' \code{\link[=viper]{viper::viper()}} for further information
df2regulon = function(df) {
  regulon = df %>%
    split(.$tf) %>%
    map(function(dat) {
      tf = dat %>% distinct(tf) %>% pull()
      targets = setNames(dat$mor, dat$target)
      likelihood = dat$likelihood
      list(tfmode = targets, likelihood = likelihood)
    })
  return(regulon)
}

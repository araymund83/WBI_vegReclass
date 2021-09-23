convertToVegTypeBCR4 <- function(DT, pureCutoff = 0.8,
                                     deciSp = c("Popu_tre", "Popu_bal","Betu_pap"),
                                     coniSp = c("Pinu_ban", "Pinu_con", "Abie_bal", "Abie_las"),
                                     wetland  = c("Pinu_mar", "Lari_lar")) {

  if (.sumRelBs(c("Pice_mar", "Pice_gla", "Abie_bal"), DT) >= pureCutoff){
    "Spruce/Fir"
    ## Spruce/Fir dominant : spruce and fir stands where the combine spruce
    ##and fir component is 80% or more
  } else if(.sumRelBs(c("Pinu_ban", "Pinu_con"), DT) >= pureCutoff){
    "Pine"
    ##Lodgepole (Pinu_con) & jackPine (Pinu_ban) 80% or more
  }
  if (.sumRelBs(deciSp, DT) >= pureCutoff){
    "Decid"
    ##Deciduous : Stands where deciduous tree component is 80% or more
  } else if(.sumRelBs(coniSp, DT) >= pureCutoff){
    "Coni"
  } else if (.sumRelBs(deciSp, DT) < pureCutoff ||
             .sumRelBs(coniSp, DT) < pureCutoff){
    "Mixed"
  }else if (.sumRelBs("Pice_mar", DT) &&
            .sumRelBs("Lari_lar", DT)  >= pureCutoff){
    "Forested Wetland"
  } else
    ## just in case there are any not covered
    NA_character_
}
#' internal function that sums relative biomasses for species matching a character string,
#' but that can be appear duplicated in another species coding column.
#' @param sppToMatch character string of species to match against for summing B.
#' @param DT data.table with columns 'speciesCode', 'relB'.
.sumRelBs <- function(sppToMatch, DT) {
  DT[speciesCode %in% sppToMatch, relB] %>% sum()
}

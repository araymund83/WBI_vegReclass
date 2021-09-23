convertToVegTypeBCR4 <- function(DT, pureCutoff = 0.8,
                                     deciSp = c("Popu_tre", "Popu_bal","Betu_pap"),
                                     coniSp = c("Pinu_ban", "Pinu_con", "Abie_bal", "Abie_las"),
                                     wetland  = c("Pinu_mar", "Lari_lar")) {

  if (.sumRelBs(c("Pice_mar", "Lari_lar", "Pinu_ban"), DT) >= pureCutoff){
    "Stunted Coniferous"
    #	Tree cover sparse at 40% or less with stunted spruce, larch and jack
    #pine species making up 80% or more.
  } else if(.sumRelBs(coniSp, DT) >= pureCutoff){
    "Coniferous mix"
    #Coniferous trees component is 80% or more.
  }
  if (.sumRelBs(coniSp, DT) <= pureCutoff){
    "Coni/Decid mix"
    ##Stands where coniferous component is not 80% or more.

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

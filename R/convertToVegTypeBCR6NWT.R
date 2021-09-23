convertToVegTypeBCR6NWT <- function(DT, pureCutoff = 0.8,
                                     deciSp = c("Popu_tre", "Popu_bal","Betu_pap"),
                                     coniSp = c("Pinu_ban", "Pinu_con"),
                                     wetland  = c("Pinu_ban", "Pinu_con")) {

  ## TODO: use factors or integers for vegClass ???
  if (.sumRelBs(C ("Pice_gla", "Pice_mar"), DT) >= pureCutoff){
    "Spruce"
  } else if(.sumRelBs(deciSp, DT) >= pureCutoff){
    "Deci"
  } else if (.sumRelBs(coniSp, DT) >= pureCutoff){
    "Coni Mix"
  } else if (.sumRelBs(deciSp, DT) < pureCutoff ||
             .sumRelBs(coniSp, DT) < pureCutoff){
    "Mixed"
  } else if (.sumRelBs("Pice_mar", DT) &&
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

convertToVegType <- function(DT, pureCutoff = 0.8,
                             pine  = c("Pinu_ban", "Pinu_con"),
                             deciSp = c("Popu_tre", "Popu_bal","Betu_pap"),
                             coniSp = c("Pinu_ban", "Pinu_con", "Abie_bal", "Abie_las")) {

  ## TODO: use factors or integers for vegClass ???
  zeroCutoff = 0.05
  if (DT$leading %in% pine &&
      .sumRelBs(deciSp, DT) < 0.20) {
    ## Pine dominant: leading species are Pinu_ban(jackpine) &/or Pinu_con (Lodgepole pine)
    ## deciduous species less than 20%
    "Pine"
  } else if (DT$leading %in% c("Pice_mar") &&
             (.sumRelBs("Lari_lar", DT) < zeroCutoff ||
              .sumRelBs(c("Popu_tre", "Popu_bal", "Abie_bal", "Pinu_ban", "Pinu_con"), DT) >= zeroCutoff)) {
    #   ## Black Spruce (BkSp): leading species is Pice_mar (black spruce)
    #   ### deciduous species is more than 0
    "BkSp"
  } else if (.sumRelBs(c("Pice_gla", "Abie_bal"), DT) >= 0.8) {
    ## White Spruce (WhSp): dominated by Pice_gla (white spruce) & Abie_bal (balsam fir) > 80%
    "WhSp"
  } else if (.sumRelBs(deciSp, DT) >= 0.8) {
    ## Deciduous: Popu_tre(trembling aspen) &/or Popu_bal (balsam poplar) &/or
    ## Betu_pap(white birch) > 80%
    "Decid"
  } else if (.sumRelBs(deciSp, DT) >= 0.2 &&
             .sumRelBs(c(coniSp, "Pice_gla", "Pice_mar"), DT) >= 0.2) {
    ## Mixedwood : Deciduous species more than 20%
    "Mixed"
  } else if ((DT$leading %in% c("Pice_mar") &&
              (.sumRelBs("Lari_lar", DT) >= zeroCutoff ||
               .sumRelBs(c("Popu_tre", "Popu_bal", "Abie_bal", "Pinu_ban", "Pinu_con"), DT) < zeroCutoff))
             ||
             DT$leading == "Lari_lar")
  { ## Wetland :Black spruce dominated and  Deciduous species more than 20%
    "BkSpWet"
  } else if (.sumRelBs(deciSp, DT) < 0.2) {
    "ConMix"
  } else if ((DT$leading %in% c("Popu_tre", "Popu_bal", "Betu_pap")) &&
             .sumRelBs(c("Pinu_ban", "Lari_lar", "Pice_mar", "Pice_gla"), DT) >= zeroCutoff) {
    "DecWet"
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



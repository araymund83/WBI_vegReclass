#' Reclassification of pixel ages into age-classes for the WB_project.
#' 
#' @param DT A `cohortData` with `pixelGroup` from LandR
#' @param ageBreaks TODO
#' @param ageLabels TODO
#' 
#' @return a `data.table` with an extra column \code{"ageGroup"}
ageReclass <- function(DT,
                       ageBreaks = c(seq(0, 140, 10), 400),
                       ageLabels = c("0-10", "11-20", "21-30", "31-40", "41-50", "51-60",
                                     "61-70","71-80","81-90", "91-100", "101-110", "111-120",
                                     "121-130", "131-140", "141+")) {
  stopifnot(length(ageLabels) == length(ageBreaks) - 1)
  DT[, ageGroup := cut(ageMax, breaks = ageBreaks, right =  FALSE, labels = ageLabels)]
  return(DT)
}

library(SpaDES.core)
library(SpaDES.tools)
outSimSppLayers <- loadSimList("/home/araymundo/GITHUB/WB_LandCoverClass/inputs/simOutSpeciesLayers_ABBC.qs")

# Ana's work below --------------------------------------------------------------------------------

## Add the pixel number to the cohortData
pixelCohortData <- LandR::addNoPixel2CohortData(outSimSppLayers$cohortData,
                                                outSimSppLayers$pixelGroupMap,
                                                doAssertion = getOption("LandR.assertions", TRUE))

## Add vegetation type column to the cohortData table
vegTypeTable <- LandR::vegTypeGenerator(pixelCohortData, vegLeadingProportion = 0.7,
                                        mixedType = 2, sppEquiv = outSimSppLayers$sppEquiv,
                                        sppEquivCol = "WB", pixelGroupColName = "pixelGroup")
vegTypeTable2 <- data.table::copy(vegTypeTable)


vegTypeMap <- LandR::vegTypeMapGenerator(pixelCohortData, vegLeadingProportion = 0.8,
                                         mixedType = 2, sppEquiv = outSimSppLayers$sppEquiv,
                                         pixelGroupMap = outSimSppLayers$pixelGroupMap,
                                         sppEquivCol = "WB", pixelGroupColName = "pixelGroup")

##subset the pixelCohortData and create a new column with the max age per pixelGroup
newAgeCD <- vegTypeTable[, list(ageMax = max(age)), by = "pixelGroup"]  ## TODO: BIOMASS WEIGHTED MEAN ?

##apply reclass function to the subset data.table
newAgeCD <- ageReclass(newAgeCD)

## create a new  RasterLayer of Age reclassified
ageGroupRas <- rasterizeReduced(reduced = newAgeCD,
                                fullRaster = outSimSppLayers[["pixelGroupMap"]],
                                mapcode = "pixelGroup", newRasterCols =  "ageGroup")
writeRaster(ageRas, "ageRas.tif")

vegTypeRas <- rasterizedReduced(reduced =)

## subset the pixelCohortData and create a new column with the sum of Biomass(B)
## per pixelGroup and RelB per pixelGroup & speciesCode
sumBCD <- vegTypeTable2[, list(sumB = sum(B)), by =  "pixelGroup"]

## create the sumB RasterLayer
sumBRas <- rasterizeReduced(reduced = sumBCD,
                            fullRaster = outSimSppLayers[["pixelGroupMap"]],
                            mapcode = "pixelGroup", newRasterCols =  "sumB")
## subset the pixelCohortData and create a new column with the sum of Biomass(B) & relB
## per pixelGroup and RelB per pixelGroup & speciesCode
vegTypeTable2[, sumB := sum(B), by = .(pixelGroup)]
vegTypeTable2[, relB := sum(B)/sumB, by = .(pixelGroup, speciesCode)]
vegTypeTable2[is.na(relB) & sumB == 0, relB := 0]

if (any(is.na(vegTypeTable2$relB))) {
  stop("Missing values in relative Biomass")
}

## subset to a smaller DT
vegTypes <- unique(vegTypeTable2[B > 0, .(pixelGroup,leading, speciesCode, B, relB, sumB)])


# a <- DT[ , list(vegClass), by = "pixelGroup"]
# rasterTemplate <- outSimSppLayers$pixelGroupMap
# rasterTemplate[!is.na(rasterTemplate)] <- 0
#
# r <- SpaDES.tools::rasterizeReduced(reduced = a,
#                                     fullRaster = outSimSppLayers[["pixelGroupMap"]],
#                                     mapcode = "pixelGroup", newRasterCols = "vegClass")

# Alex testing --------------------------------------------------------------------------------
DT <- data.table::copy(vegTypes)
data.table::setkeyv(DT, cols = "pixelGroup")

## TODO: use DT instead of DT2
DT2 <- DT[pixelGroup >= 3200000,]
#DT2 <- DT[pixelGroup == 500,]
DT2[, vegClass := convertToVegType(.SD, pureCutoff = 0.5,
                                   deciSp = c("Popu_tre", "Popu_bal","Betu_pap"),
                                   coniSp = c("Pinu_ban", "Pinu_con", "Abie_bal", "Abie_las")),
   by = pixelGroup, .SDcols = c("leading", "speciesCode", "relB")]


DTNA <- DT[is.na(vegClass)]

# end Alex testing ----------------------------------------------------------------------------

################################################################################
#                        non-forested pixels
################################################################################
uniqueLCCClasses <- na.omit(unique(simOutPreamble$rstLCC))
nontreeClasses <- sort(uniqueLCCClasses[!uniqueLCCClasses %in% forestedLCCClasses])
treepixels <- simOutPreamble$rstLCC[] %in% forestedLCCClasses
nontreePixels <- which(simOutPreamble$rstLCC[] %in% nontreeClasses)
vals <- raster::getValues(simOutPreamble$rstLCC)
### get all NA's
lcc05NA <- is.na(simOutPreamble$rstLCC[])




nonForestRas <- raster::setValues(x = simOutPreamble$rstLCC,)

nonForesReclassTB <- Cache(prepInputs, url = paste0("https://drive.google.com/file/",
                                                    "d/17IGN5vphimjWjIfyF7XLkUeD-ze",
                                                    "Kruc1/view?usp=sharing"),
                           destinationPath = Paths$inputPath,
                           purge = 7,
                           fun = "data.table::fread",
                           userTags = "WBnonForest_LCC05")

reclassMatrix <- usefulFuns::makeReclassifyMatrix(table = nonForesReclassTB,
                                                  originalCol = "LCC05_Class",
                                                  reclassifiedTo = "nonForest_Class")
nonForestRas <- raster::reclassify(x = simOutPreamble$rstLCC, rcl = reclassMatrix[, -1])


# # get xy of all pixels that are not forested classes
nontreeIndex <- which(simOutPreamble$rstLCC[] %in% nontreeClasses)

# #extract pixel numbers of all xy from LCC05
 nontreeLocations <- xyFromCell(simOutPreamble$rstLCC, nontreeIndex)
 lcc05nonForest <-  as.data.table(raster::extract(simOutPreamble$rstLCC, nontreeLocations, cellnumbers = TRUE))

 lccnonForestLayer <- nontreePixels
 lccnonForestLayer[!is.na(lccnonforest)]

 # countnonForest <- lcc05nonForest[, .N, by= cells]
#
# lcc05Vals <- data.table(LCC = getValues(simOutPreamble$rstLCC), pixelID = 1:ncell(simOutPreamble$rstLCC))
# nonForestVals <- lcc05Vals[, LCC > 15]
# nonForestRas <- reclassify()
#
# nontreeVals <- lcc05Vals[, lcc > 15]
#
#
# lnontreeLocations<- xyFromCell(simOutPreamble$rstLCC, nontreeClasses)
lcc05nonForest <- data.table(lcc = getValues(simOutPreamble$rstLCC), pixelID = 1:ncell(simOutPreamble$rstLCC))

remapDT <- data.table::as.data.table(expand.grid(LCC2005 = c(NA_integer_, sort(uniqueLCCClasses)),
                                                 nonForest))
remapDT[LCC2005 == "NA", nonForest := 99]
remapDT[LCC2005 == 0, nonForest := NA_integer_]
remapDT[LCC2005 %in% c(1:15), nonForest := 0]
remapDT[LCC2005 %in% c(37,38,39), nonForest := 1]
remapDT[LCC2005 == 19, nonForest := 2]
remapDT[LCC2005 %in% c(25,33,36), nonForest := 3]
remapDT[LCC2005 %in% c(17,18,21,23,26:29), nonForest := 4]
remapDT[LCC2005 %in% c(16,20,22), nonForest := 5]
remapDT[LCC2005 %in% c(24,30:32), nonForest := 6]
remapDT[LCC2005 %in% treeClassesToReplace, nonForest := 0]

lcc05nonForest <- lcc05nonForest[remapDT, on = c("lcc" = "LCC2005")]

lcc05nonForestRas <- raster::setValues(x =simOutPreamble$rstLCC, values = lcc05nonForest)

message("Overlaying land cover maps...")
LCClarge <- Cache(overlayLCCs,
                      LCCs = list( LCC2005 = simOutPreamble$rstLCC),
                      forestedList = list(nonForest = 0, LCC2005 = forestedLCCClasses),
                      outputLayer = "LCC2005",
                      remapTable = remapDT,
                      classesToReplace = c(treeClassesToReplace, 99),
                      availableERC_by_Sp = NULL)
message("...done.")
treePixelsLCC <- which(sim$LCClarge[] %in% P(sim)$treeClassesLCC)
nonTreePixels <- which(sim$LCClarge[] %in% nontreeClassesLCC)

sim$nonTreePixels <- nonTreePixels
## Update rasterToMatch layer with all trees
ml[[ml@metadata[ml@metadata$rasterToMatch == 1, ]$layerName]][sim$nonTreePixels] <- NA
sim$rasterToMatch <- postProcess(rasterToMatch(ml), studyArea = sim$studyArea, filename2 = NULL)





## Everything in this file and any files in the R directory are sourced during `simInit()`;
## all functions and objects are put into the `simList`.
## To use objects, use `sim$xxx` (they are globally available to all modules).
## Functions can be used inside any function that was sourced in this module;
## they are namespaced to the module, just like functions in R packages.
## If exact location is required, functions will be: `sim$.mods$<moduleName>$FunctionName`.
defineModule(sim, list(
  name = "WBI_vegReclass",
  description = paste("the aim of this module is to make a reclassification of cohortData into predefined vegetation classes for
  six subpregions whitin the WBI project"),
  keywords = c("WB", "habitat types"),
  authors =  c(
    person("Alex M", "Chubaty", email= "achubaty@for-cast.ca", role = "aut"),
    person("Ana", "Raymundo", email= "angeles-ana-paula.raymundo-sanchez.1@ulaval.ca", role = "aut")
  ),
  childModules = character(0),
  version = list(SpaDES.core = "1.0.6.9020", WBI_vegReclass = "0.0.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = deparse(list("README.txt", "WBI_vegReclass.Rmd")),
  reqdPkgs = list("fasterize", "SpaDES.tools", "googledrive", "data.table", "raster", "reproducible", "LandR"),
  parameters = rbind(
    defineParameter(".plots", "character", "screen", NA, NA,
                    "Used by Plots function, which can be optionally used here"),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA,
                    "Describes the simulation time at which the first plot event should occur."),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "Describes the simulation time interval between plot events."),
    defineParameter("bcr", "numeric", 6, NA, NA,
                    paste("BCR region used for subdivision of the WB study area when",
                          "using reclassifying land cover classes.Options are:",
                          "4", "6", "7", "8",
                          "bcr4BC", "bcr6BC", "bcr6AB", "bcr6SK", "bcr7SK", "bcr8SK",
                          "bcr6MB", "bcr7MB", "bcr8MB", "bcr4NT", "bcr6NT", "bcr7NT",
                          "bcr7NU", "bcr4YU")),
    defineParameter("reclassTime", "numeric", 1, NA, NA,
                    "Describes the simulation time at which the reclassify event should occur."),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "Describes the simulation time at which the first save event should occur."),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events."),
    defineParameter("studyAreaName", "character", "BC", NA, NA,
                    paste("study area name for each of the provinces in WB. Options are: BC, AB",
                          "SK", "MB", "NT", "NU", "YU")),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    paste("Should caching of events or module be activated?",
                          "This is generally intended for data-type modules, where stochasticity",
                          "and time are not relevant"))
  ),
  inputObjects = bindrows(
    expectsInput("cohortData",  "data.table",
                 desc = paste("Initial community table, created from available biomass (g/m2)",
                              "age and species cover data, as well as ecozonation information",
                              "Columns: B, pixelGroup, speciesCode")),
    expectsInput("sppEquivCol", "data.table",
                 desc = "The column in sim$speciesEquivalency data.table to use as a naming convention" ),
    expectsInput("sppEquiv", "data.table",
                 desc = "The column in sim$speciesEquivalency data.table to use as a naming convention" ),
    expectsInput("studyArea", objectClass = "SpatialPolygonsDataFrame",
                  desc = "study area used for REPORTING. This shapefile is created in the WBI preamble module"),
    expectsInput("pixelGroupMap", "RasterLayer",
                 desc = "Initial community map that has mapcodes match initial community table"),
    expectsInput("rstLCC", "RasterLayer",
                 desc = "Initial Land Cover 2005 classes")

    ),

  outputObjects = bindrows(
    createsOutput("vegTypesRas", "RasterLayer",
                  desc = paste("reclassification of cohort data into pre-defined",
                               "vegetation classes for the WBI project")),
    createsOutput("nonForestRas", "RasterLayer",
                  desc = "reclassification of non forested pixels"),

    createsOutput("ageRas", "RasterLayer",
                  desc = "reclassification of age from cohorData"),


  )
))

## event types
#   - type `init` is required for initialization

doEvent.WBI_vegReclass = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      # schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$reclassTime, "WBI_vegReclass", "reclass")
      #sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "WBI_vegReclass", "plot")
      #sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "WBI_vegReclass", "save")
    },

   reclass = {

      bcrSApixels<- raster::extract(sim$pixelGroupMap,
                                    sim$studyArea,
                                    cellnumbers = TRUE, df = TRUE)

      ##assign the BCR number to the ID
      bcrSApixels$ID <- as.factor(bcrSApixels$ID)
      sim$studyArea$BCR <- as.factor(sim$studyArea$BCR)
      levels(bcrSApixels$ID) <- levels(sim$studyArea$BCR)

      ## rename columns for easier reference
      names(bcrSApixels) <- c("bcrID", "pixelID", "pixelGroup")
      bcrSApixels <- na.omit(bcrSApixels)

      ##merge both DT to have the BCR id in the cohortData
      pixelCohortData <- LandR::addNoPixel2CohortData(sim$cohortData,
                                                      sim$pixelGroupMap,
                                                      doAssertion = getOption('LandR.assertions', TRUE))

      ##merge both DT to have the BCR id in the cohortData
      bcrCohortData<- merge(bcrSApixels, pixelCohortData, by = 'pixelGroup', all = TRUE)
      bcrCohortData <- as.data.table(bcrCohortData)

      ## Add vegetation type column to the bcr cohortData table
      vegTypeTable <- LandR::vegTypeGenerator(bcrCohortData, vegLeadingProportion = 0.8,
                                              mixedType = 2, sppEquiv = sim$sppEquiv,
                                              sppEquivCol = sim$sppEquivCol,
                                              pixelGroupColName = 'pixelGroup')

      ## subset the pixelCohortData and create a new column with the sum of Biomass(B) & relB
      ## per pixelGroup and RelB per pixelGroup & speciesCode
      vegTypeTable[, sumB := sum(B), by = .(pixelGroup)]
      vegTypeTable[, relB := sum(B)/sumB, by = .(pixelGroup, speciesCode)]
      vegTypeTable[is.na(relB) & sumB == 0, relB := 0]

      ## check for missing values in B
      if (any(is.na(vegTypeTable$relB))) {
        stop("Missing values in relative Biomass")
      }

      ## subset to a smaller DT
      vegTypes <- unique(vegTypeTable[B > 0, .(pixelGroup,leading, speciesCode, bcrID, relB, sumB)])

      #DT <- data.table::copy(vegTypes)
      data.table::setkeyv(vegTypes, cols = "pixelGroup")

      ##subset cohortData for an specific bcr
      bcrSelect <- P(sim)$bcr

     # vegTypes <- vegTypes[vegTypes$bcrID == bcrSelect]
      vegTypes <- vegTypes[bcrID %in% bcrSelect ]

      ##subset pixelGroupMap for an specific bcr
      ## get pixelGroup per bcr polygon
      sim$studyArea <-st_as_sf(sim$studyArea)
      bcrPixelGroupRas <- fasterize::fasterize(sim$studyArea, sim$pixelGroupMap, field = "BCR")

      bcrSelectSA <- sim$studyArea[sim$studyArea$BCR %in% bcrSelect, ]

      ## reclass rules per studyArea and
      if (P(sim)$studyAreaName == "AB" || P(sim)$studyAreaName == "BC"
          && P(sim)$bcr == 6){
        vegTypes[, vegClass:= convertToVegType(.SD, pureCutoff = 0.8,
                                               deciSp = c("Popu_tre", "Popu_bal", "Betu_pap"),
                                               coniSp = c("Pinu_ban", "Pinu_con", "Abie_bal", "Abie_las")),
                 by = pixelGroup, .SDcols = c("speciesCode", "relB", "leading")]

      } else if (P(sim)$studyAreaName == "BC" && P(sim)$bcr == 4){
        vegTypes[, vegClass:= convertToVegTypeBCR4(.SD, pureCutoff = 0.8,
                                                   deciSp = c("Popu_tre", "Popu_bal", "Betu_pap"),
                                                   coniSp = c("Pinu_ban", "Pinu_con", "Abie_bal", "Abie_las"),
                                                   wetland = c("Pinu_ban", "Pinu_con")),
                 by = pixelGroup, .SDcols = c("speciesCode", "relB", "leading", "bcrID")]

      } else if (P(sim)$studyAreaName == "MB" || P(sim)$studyAreaName == "SK"
                 && P(sim)$bcr == 6){
        vegTypes[, vegClass:= convertToVegTypeBCR6MBSK(.SD, pureCutoff = 0.8,
                                                       deciSp = c("Popu_tre", "Popu_bal","Betu_pap"),
                                                       coniSp = c("Pinu_ban", "Pinu_con", "Abie_bal", "Abie_las"),
                                                       wetland = c("Pinu_ban", "Pinu_con")),
                 by = pixelGroup, .SDcols = c("speciesCode", "relB", "leading", "bcrID")]
      } else if (P(sim)$studyAreaName == "MB" && P(sim)$bcr == 7){
        vegTypes[, vegClass:= convertToVegTypeBCR7(.SD, pureCutoff = 0.8,
                                                   deciSp = c("Popu_tre", "Popu_bal", "Betu_pap"),
                                                   coniSp = c("Pinu_ban", "Pinu_con", "Abie_bal", "Abie_las"),
                                                   wetland  = c("Pinu_mar", "Lari_lar")),
                 by = pixelGroup, .SDcols = c("speciesCode", "relB", "leading", "bcrID")]
      } else if (P(sim)$studyAreaName == "MB" && P(sim)$bcr == 8){
        vegTypes[, vegClass:= convertToVegTypeBCR8(.SD, pureCutoff = 0.8,
                                                   deciSp = c("Popu_tre", "Popu_bal", "Betu_pap"),
                                                   coniSp = c("Pinu_ban", "Pinu_con", "Abie_bal", "Abie_las"),
                                                   wetland  = c("Pinu_mar", "Lari_lar")),
                 by = pixelGroup, .SDcols = c("speciesCode", "relB", "leading", "bcrID")]


      }else if (P(sim)$studyAreaName == "NT" && P(sim)$bcr == 4){
        vegTypes[, vegClass:= convertToVegTypeBCR4(.SD, pureCutoff = 0.8,
                                                   deciSp = c("Popu_tre", "Popu_bal", "Betu_pap"),
                                                   coniSp = c("Pinu_ban", "Pinu_con", "Abie_bal", "Abie_las"),
                                                   wetland = c("Pinu_ban", "Pinu_con")),
                 by = pixelGroup, .SDcols = c("speciesCode", "relB", "leading", "bcrID")]

      }

      vegTypes$vegClass <- as.factor(vegTypes$vegClass)

      vegTypesCD <- vegTypes[, list(vegType = unique(vegClass)), by = "pixelGroup"]


      vegTypesRas <- SpaDES.tools::rasterizeReduced(reduced = vegTypesCD,
                                                    fullRaster = sim$pixelGroupMap,
                                                    mapcode = "pixelGroup", newRasterCols ="vegType")

      sim$vegTypesRas <- vegTypesRas
      ######NON FORESTED PIXELS#########
      pixelCohortData2 <- LandR::addPixels2CohortData(sim$cohortData, sim$pixelGroupMap)
      ## get the non forested pixelIndices
      nonForestLCC <- data.table(pixelID= 1:ncell(sim$rstLCC),
                                 pg = getValues(sim$pixelGroupMap),
                                 LCC05_Class = getValues(sim$rstLCC))

      nonForestLCC <- nonForestLCC[is.na(pg) & !is.na(LCC05_Class)]

      ## This table has the non forested values and their description
      # nonForestReclassTB <- Cache(prepInputs, url = paste0("https://drive.google.com/file/",
      #                                                     "d/17IGN5vphimjWjIfyF7XLkUeD-ze",
      #                                                     "Kruc1/view?usp=sharing"),
      #                            destinationPath = Paths$inputPath,
      #                            fun = "data.table::fread",
      #                            userTags = "WBnonForest_LCC05")


      nonForestReclassTB <- fread(file.path("inputs", "nonForest_reclass.csv"))

      setkey(nonForestLCC, LCC05_Class)
      setkey(nonForestReclassTB, LCC05_Class)

      ## join both tables based on the LCC05 value
      nonForestPixLCC<-merge(nonForestLCC, nonForestReclassTB)
      ##create an empty raster with the same number of pixels  than pixelGroupMap
      templateRas <- raster(sim$pixelGroupMap)
      ##setup all values to NA
      templateRas[]<- NA
      ## assign nonForested reclassified values
      templateRas[nonForestPixLCC$pixelID]<- nonForestPixLCC$nonForest_Class

      sim$nonForestRas <- templateRas
      ##crop & mask the raster to the BCR selected area.
      nonForestRas <- postProcess( templateRas, studyArea = bcrSelectSA)

      sim$nonForestRas <- nonForestRas


      ######RECLASSIFY AGES#########
      ##subset the pixelCohortData and create a new column with the max age per pixelGroup
      newAgeCD <- vegTypeTable[, list(ageMax = max(age)), by = "pixelGroup"]

      ## create a new  RasterLayer of Age reclassified
      ageRas <- SpaDES.tools::rasterizeReduced(reduced = newAgeCD,
                                               fullRaster = sim$pixelGroupMap,
                                               mapcode = "pixelGroup", newRasterCols = "ageMax")

      ageRas <- postProcess(ageRas, studyArea = bcrSelectSA)
      sim$ageRas <- ageRas

      return(list(sim$vegTypesRas,sim$nonForestRas,
                  sim$ageRas))

      sim <- scheduleEvent(sim, time(sim) + P(sim)$reclassTime, "vegReclass_WBI", "reclass")
   },

# plot = {
# Plot(sim$vegTypeMap) ##I don't think this is necessary ASK STEVE

#sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "vegReclass_WBI", "plot")
# },
warning(paste("Undefined event type: \'", current(sim)[1, "eventType", with = FALSE],
              "\' in module \'", current(sim)[1, "moduleName", with = FALSE], "\'", sep = ""))
  )
  return(invisible(sim))
}

.inputObjects <- function(sim) {

  #cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")
  # if (!suppliedElsewhere("studyArea", sim)) {
  #   message("study area not supplied. Using AB province within WB studyARea")
  #
  #   sim$studyArea <- Cache(prepInputs, url = extractURL(objectName = "studyArea", sim = sim),
  #                          destinationPath = getPaths()$inputPath,
  #                          filename2 = "studyArea", overwrite = TRUE)
  # }
  # if(!suppliedElsewhere("pixelGroupMap", sim)){
  #   message("pixelGrouMap not supplied. Please provide one")
  #   sim$pixelGroupMap <- Cache(LandR::prepInputsLCC,
  #                              year = 2005,
  #                              destinationPath = Paths$inputPath,
  #                              studyArea = sim$studyArea,
  #                              filename2 = "RTM")
  #}

  return(invisible(sim))
}



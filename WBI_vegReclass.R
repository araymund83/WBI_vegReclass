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
  reqdPkgs = list("SpaDES.tools", "googledrive", "data.table", "raster", "reproducible", "LandR"),
  parameters = rbind(
    defineParameter(".plots", "character", "screen", NA, NA,
                    "Used by Plots function, which can be optionally used here"),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA,
                    "Describes the simulation time at which the first plot event should occur."),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "Describes the simulation time interval between plot events."),
    defineParameter("bcr", "character", NA, NA, NA,
                    paste("BCR region used for subdivision of the WB study area when",
                          "using reclassifying land cover classes.Options are:",
                          "bcr4", "bcr6", "bcr7", "bcr8",
                          "bcr4BC", "bcr6BC", "bcr6AB", "bcr6SK", "bcr7SK", "bcr8SK",
                          "bcr6MB", "bcr7MB", "bcr8MB", "bcr4NT", "bcr6NT", "bcr7NT",
                          "bcr7NU", "bcr4YU")),
    defineParameter("reclassTime", "numeric", 1, NA, NA,
                    "Describes the simulation time at which the reclassify event should occur."),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "Describes the simulation time at which the first save event should occur."),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events."),
    # defineParameter("studyAreaName", "character", "bcr6BC", NA, NA,
    #                 paste("study area name for WB project, options are: bcr6BC, bcr6AB",
    #                       "bcr6SK,YK, bcr6NWT,bcr6MB")),
    defineParameter("studyArea", "SpatialPolygonsDataFrame", "BC", NA, NA,
                    paste("study area shapefile for each of the provinces in WB. Options are: BC, AB",
                          "SK", "MB", "NT", "NU", "YU")),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    paste("Should caching of events or module be activated?",
                          "This is generally intended for data-type modules, where stochasticity",
                          "and time are not relevant"))
  ),
  inputObjects = bindrows(
    expectsInput("cohortData",  "data.table",
                 desc = paste("Initial community table, created from available biomass (g/m2"),
                             ("age and species cover data, as well as ecozonation information"),
                             ("Columns: B, pixelGroup, speciesCode")),
    # expectsInput("rstLCC", "RasterLayer",
    #              desc = paste("Initial conditions from LCC 2005 ",
    #                           "XXXXX")),
    expectsInput("sppEquiv", "data.table",
                 desc = "The column in sim$speciesEquivalency data.table to use as a naming convention" ),
    expectsInput("sppEquiv", "data.table",
                 desc = "The column in sim$speciesEquivalency data.table to use as a naming convention" ),
    expectsInput("studyArea", objectClass = "SpatialPolygonsDataFrame",
                  desc = "study area used for REPORTING. This shapefile is created in the WBI preamble module"),
    expectsInput("pixelGroupMap", "RasterLayer",
                 desc = "Initial community map that has mapcodes match initial community table")
    ),

  outputObjects = bindrows(
    createsOutput("vegTypesMap", "RasterLayer",
                  desc = paste("reclassification of cohort data into pre-defined",
                               "vegetation classes for the WBI project"),
                  "nonForestedMap", "RasterLayer",
                  desc = paste("reclassification of non forested pixels"),
                  "ageRas", "RasterLayer",
                  desc = paste("ageMap raster"))
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
     browser()
      bcrSApixels<- raster::extract(sim$pixelGroupMap,
                                    sim$studyArea,
                                    cellnumbers = TRUE, df = TRUE)

      ##assign the BCR number to the ID
      bcrSApixels$ID <- as.factor(bcrSApixels$ID)
      sim$studyArea$BCR <- as.factor(sim$studyArea$BCR)
      levels(bcrSApixels$ID) <- levels(sim$studyArea$BCR)

      ## rename columns for easier reference
      names(bcrSApixels) <- c ("bcrID", "pixelID", "pixelGroup")


      ##merge both DT to have the BCR id in the cohortData

      pixelCohortData <- LandR::addNoPixel2CohortData(cohortData,
                                                      pixelGroupMap,
                                                      doAssertion = getOption("LandR.assertions", TRUE))

      ##merge both DT to have the BCR id in the cohortData
      #bcrpixelCohortData<- merge(pixelCohortData, bcrMBpixels, by = "pixelGroup")
      bcrCohortData<- merge(pixelCohortData, bcrSApixels, by = "pixelGroup")

      ## Add vegetation type column to the bcr cohortData table
      vegTypeTable <- LandR::vegTypeGenerator(bcrCohortData, vegLeadingProportion = 0.8,
                                              mixedType = 2, sppEquiv = sim$sppEquiv,
                                              sppEquivCol = sim$sppEquivCol, pixelGroupColName = "pixelGroup")


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
      bcrSelect <- sim$bcr
      VegTypes<- vegTypes[vegTypes$bcrID %in% bcrSelect, ]

      ##subset pixelGroupMap for an specific bcr
      ## get pixelGroup per bcr polygon
      sim$studyArea <-st_as_sf(sim$studyArea)
      bcrPixelGroupRas <- fasterize(sim$studyArea, sim$pixelGroupMap, field = "BCR")

      bcrSelectSA <- sim$studyArea[sim$studyArea$BCR %in% bcrSelect, ]
      subsetpixelGroupMapRas <- crop(sim$pixelGroupMap, bcrSelectSA)
      bcrPixelGroupRas <- mask(subsetpixelGroupMapRas, mask = bcrSelectSA)


      if (grepl("AB", P(sim)$studyArea)){
          if(grepl("6", P(sim)$bcr)){
      vegTypes[, vegClass:= convertToVegType(.SD, pureCutoff = 0.8,
                                             deciSp = c("Popu_tre", "Popu_bal", "Betu_pap"),
                                             coniSp = c("Pinu_ban", "Pinu_con", "Abie_bal", "Abie_las")),
                   by = pixelGroup, .SDcols = c("speciesCode", "relB", "leading")]
          }
      } else if (grepl ("BC", P(sim)$studyArea)){
        if(grepl("4", P(sim)$bcr)){
          vegTypes[, vegClass:= convertToVegTypeBCR4(.SD, pureCutoff = 0.8,
                                                     deciSp = c("Popu_tre", "Popu_bal", "Betu_pap"),
                                                     coniSp = c("Pinu_ban", "Pinu_con", "Abie_bal", "Abie_las"),
                                                     wetland = c("Pinu_ban", "Pinu_con")),
                   by = pixelGroup, .SDcols = c("speciesCode", "relB", "leading", "bcrID")]

        } else if (grepl ("6", P(sim)$bcr)){
          vegTypes[, vegClass:= convertToVegTypeBCR6MBSK(.SD, pureCutoff = 0.8,
                                                         deciSp = c("Popu_tre", "Popu_bal","Betu_pap"),
                                                         coniSp = c("Pinu_ban", "Pinu_con", "Abie_bal", "Abie_las"),
                                                         wetland = c("Pinu_ban", "Pinu_con")),
                   by = pixelGroup, .SDcols = c("speciesCode", "relB", "leading", "bcrID")]
          }
      }else if (P(sim)$studyArea == "MB" || P(sim)$studyArea == "SK"){
        if(grepl("6", P(sim)$bcr)){
          vegTypes[, vegClass:= convertToVegTypeBCR6MBSK(.SD, pureCutoff = 0.8,
                                                         deciSp = c("Popu_tre", "Popu_bal","Betu_pap"),
                                                         coniSp = c("Pinu_ban", "Pinu_con", "Abie_bal", "Abie_las"),
                                                         wetland = c("Pinu_ban", "Pinu_con")),
                   by = pixelGroup, .SDcols = c("speciesCode", "relB", "leading", "bcrID")]
        }else if (grepl ("7", P(sim)$bcr)){
          vegTypes[, vegClass:= convertToVegTypeBCR7(.SD, pureCutoff = 0.8,
                                                     deciSp = c("Popu_tre", "Popu_bal", "Betu_pap"),
                                                     coniSp = c("Pinu_ban", "Pinu_con", "Abie_bal", "Abie_las"),
                                                     wetland  = c("Pinu_mar", "Lari_lar")),
                   by = pixelGroup, .SDcols = c("speciesCode", "relB", "leading", "bcrID")]
        } else if (grepl ("8", P(sim)$bcr)){
          vegTypes[, vegClass:= convertToVegTypeBCR8(.SD, pureCutoff = 0.8,
                                                     deciSp = c("Popu_tre", "Popu_bal", "Betu_pap"),
                                                     coniSp = c("Pinu_ban", "Pinu_con", "Abie_bal", "Abie_las"),
                                                     wetland  = c("Pinu_mar", "Lari_lar")),
                   by = pixelGroup, .SDcols = c("speciesCode", "relB", "leading", "bcrID")]

        }
        }else if (P(sim)$studyArea == "NT" || P(sim)$studyArea == "NU"){
          if (grepl("4", P(sim)$bcr)){
            vegTypes[, vegClass:= convertToVegTypeBCR4(.SD, pureCutoff = 0.8,
                                                       deciSp = c("Popu_tre", "Popu_bal", "Betu_pap"),
                                                       coniSp = c("Pinu_ban", "Pinu_con", "Abie_bal", "Abie_las"),
                                                       wetland = c("Pinu_ban", "Pinu_con")),
                     by = pixelGroup, .SDcols = c("speciesCode", "relB", "leading", "bcrID")]

          }
        }
      vegTypes$vegClass <- as.factor(vegTypes$vegClass)

      vegTypesCD <- vegTypes[, list(vegType = unique(vegClass)), by = "pixelGroup"]


      vegTypesRas <- SpaDES.tools::rasterizeReduced(reduced = vegTypesCD,
                                      fullRaster = sim$pixelGroupMap,
                                      mapcode = "pixelGroup", newRasterCols ="vegType")

      sim$vegTypesRas <- vegTypesRas


      ######NON FORESTED PIXELS#########
      nonForesReclassTB <- Cache(prepInputs, url = paste0("https://drive.google.com/file/",
                                                          "d/17IGN5vphimjWjIfyF7XLkUeD-ze",
                                                          "Kruc1/view?usp=sharing"),
                                 destinationPath = Paths$inputPath,
                                 fun = "data.table::fread",
                                 userTags = "WBnonForest_LCC05")

      reclassMatrix <- usefulFuns::makeReclassifyMatrix(table = nonForesReclassTB,
                                                        originalCol = "LCC05_Class",
                                                        reclassifiedTo = "nonForest_Class")
      nonForestRas <- raster::reclassify(x = rstLCC, rcl = reclassMatrix[, -1])



      ##subset the pixelCohortData and create a new column with the max age per pixelGroup
      newAgeCD <- vegTypeTable[, list(ageMax = max(age)), by = "pixelGroup"]  ## TODO: BIOMASS WEIGHTED MEAN ?

      ## create a new  RasterLayer of Age reclassified
      ageRas <- SpaDES.tools::rasterizeReduced(reduced = newAgeCD,
                                 fullRaster = sim$pixelGroupMap,
                                 mapcode = "pixelGroup", newRasterCols = "ageMax")
      browser()
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
  #   message("study area not supplied. Using Ecodistrict 644")
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



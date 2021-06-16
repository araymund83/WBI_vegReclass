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
    defineParameter("reclassTime", "numeric", 1, NA, NA,
                    "Describes the simulation time at which the reclassify event should occur."),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "Describes the simulation time at which the first save event should occur."),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events."),
    defineParameter("studyAreaName", "character", "bcr6BC", NA, NA,
                    paste("study area name for WB project, options are: bcr6BC, bcr6AB",
                          "bcr6SK,YK, bcr6NWT,bcr6MB")),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    paste("Should caching of events or module be activated?",
                          "This is generally intended for data-type modules, where stochasticity",
                          "and time are not relevant"))
  ),

  inputObjects = bindrows(
    expectsInput("cohortData",  "data.table",
                 desc = paste0("Initial community table, created from available biomass (g/m2)",
                              "age and species cover data, as well as ecozonation information",
                              "Columns: B, pixelGroup, speciesCode")),
    expectsInput("rstLCC", "RasterLayer",
                 desc = paste("Initial conditions from LCC 2005 ",
                              "XXXXX")),
    expectsInput("sppEquiv", "data.table",
                 desc = "The column in sim$speciesEquivalency data.table to use as a naming convention" ),
    expectsInput("sppEquivCol", "character",
                 desc = "The column in sim$speciesEquivalency data.table to use as a naming convention" ),
    expectsInput("pixelGroupMap", "RasterLayer",
                 desc = "Initial community map that has mapcodes match initial community table")
  ),

  outputObjects = bindrows(
    createsOutput("forestTypeMap", "RasterLayer",
                  desc = paste0("reclassification of cohort data into pre-defined",
                               "vegetation classes for the WBI project")),
    createsOutput("nonForestedMap", "RasterLayer",
                  desc = "reclassification of non forested pixels"),
    createsOutput("ageRas", "RasterLayer",
                   desc = paste0("ageMap raster"))
  )
))


## event types
#   - type `init` is required for initialization

doEvent.WBI_vegReclass = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
   #   browser()
      # schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$reclassTime, "WBI_vegReclass", "reclass")
      #sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "WBI_vegReclass", "plot")
      #sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "WBI_vegReclass", "save")
    },

    reclass = {
      ## Add the pixel number to the cohortData
      pixelCohortData <- LandR::addNoPixel2CohortData(sim$cohortData,
                                                      sim$pixelGroupMap,
                                                      doAssertion = getOption("LandR.assertions", TRUE))
      ## Create leading species column as the species with the highest B
      vegTypeTable <- LandR::vegTypeGenerator(pixelCohortData, vegLeadingProportion = 0,
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
      vegTypes <- unique(vegTypeTable[B > 0, .(pixelGroup,leading, speciesCode, relB, sumB)])

      DT <- data.table::copy(vegTypes)
      data.table::setkeyv(DT, cols = "pixelGroup")

      vegTypes[, vegClass:= convertToVegType(.SD, pureCutoff = 0.8,
                                                 deciSp = c("Popu_tre", "Popu_bal","Betu_pap"),
                                                 coniSp = c("Pinu_ban", "Pinu_con", "Abie_bal", "Abie_las")),
                   by = pixelGroup, .SDcols = c("speciesCode", "relB", "leading")]

      vegTypes$vegClass <- as.factor(vegTypes$vegClass)

      vegTypesCD <- vegTypes[, list(vegType = unique(vegClass)), by = "pixelGroup"]

      vegTypesRas <- SpaDES.tools::rasterizeReduced(reduced = vegTypesCD,
                                      fullRaster = sim$pixelGroupMap,
                                      mapcode = "pixelGroup", newRasterCols ="vegType")

      sim$vegTypesRas <- vegTypesRas
      browser()

      nonForesReclassTB <- Cache(prepInputs, url = paste0("https://drive.google.com/file/",
                                                          "d/17IGN5vphimjWjIfyF7XLkUeD-ze",
                                                          "Kruc1/view?usp=sharing"),
                                 destinationPath = Paths$inputPath,
                                 fun = "data.table::fread",
                                 userTags = "ABnonForest_LCC05")

      reclassMatrix <- usefulFuns::makeReclassifyMatrix(table = nonForesReclassTB,
                                                        originalCol = "LCC05_Class",
                                                        reclassifiedTo = "nonForest_Class")
      nonForestRas <- raster::reclassify(x = sim$rstLCC, rcl = reclassMatrix[, -1])
      sim$nonForestRas <- nonForestRas



      ##subset the pixelCohortData and create a new column with the max age per pixelGroup
      newAgeCD <- vegTypeTable[, list(ageMax = max(age)), by = "pixelGroup"]  ## TODO: BIOMASS WEIGHTED MEAN ?

      ## create a new  RasterLayer of Age reclassified
      ageRas <- SpaDES.tools::rasterizeReduced(reduced = newAgeCD,
                                 fullRaster = sim$pixelGroupMap,
                                 mapcode = "pixelGroup", newRasterCols = "ageMax")

      sim$ageRas <- ageRas

      return(list(sim$vegTypesRas, sim$nonForestRas,
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
  # if(!suppliedElsewhere("rasterToMatch", sim)){
  #   message("rasterToMatch not supplied. Generating from LCC05")
  #   sim$rasterToMatch <- Cache(LandR::prepInputsLCC,
  #                              year = 2005,
  #                              destinationPath = Paths$inputPath,
  #                              studyArea = sim$studyArea,
  #                              filename2 = "RTM")
  #}

  return(invisible(sim))
}



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
  reqdPkgs = list("googledrive", "data.table", "raster", "reproducible", "LandR"),
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
                 desc = paste("Columns: B, pixelGroup, speciesCode, Indicating several features",
                              "about ages and current vegetation of stand")),
    #expectsInput("studyArea", "SpatialPoligonsDataFrame",
     #                         desc = paste("Polygon to use as the study area"),
      #                        sourceURL = "https://drive.google.com/open?id=18XPcOKeQdty102dYHizKH3ZPE187BiYi"),
    expectsInput("rasterToMatch", objectClass = "RasterLayer",
                              desc = "RasterToMatch" ),
    expectsInput("sppEquiv", objectClass = "data.table",
                              desc = "xxxxxxRasterToMatch" ),
    expectsInput("sppEquivCol", objectClass = "character",
                              desc = "xxxxxxRasterToMatch" ),
    expectsInput("pixelGroupMap", objectClass = "RasterLayer",
                 desc = NA, sourceURL = NA)
    ),

  outputObjects = bindrows(
    createsOutput("vegTypeMap", "Rasterlayer",
                  desc = paste("reclassification of biomass map into pre-defined",
                               "vegetation classes for the WBI project"))
  )
))

## event types
#   - type `init` is required for initialization

doEvent.WBI_vegReclass = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      browser()
      sim <- Init(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$reclassTime, "WBI_vegReclass", "reclass")
     # sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "WBI_vegReclass", "plot")
      #sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "WBI_vegReclass", "save")
    },
    reclass = {

      #sim$groupAgeMap  <- ageReclass(cohortData = sim$cohortData)
      ## Add the pixel number to the cohortData
      pixelCohortData <- LandR::addNoPixel2CohortData(sim$cohortData,
                                                      sim$pixelGroupMap,
                                                      doAssertion = getOption("LandR.assertions", TRUE))

      vegTypeTable <- LandR::vegTypeGenerator(pixelCohortData, vegLeadingProportion = 0,
                                              mixedType = 2, sppEquiv = sim$sppEquiv,
                                              sppEquivCol = "WB", pixelGroupColName = "pixelGroup")

      ## subset the pixelCohortData and create a new column with the sum of Biomass(B) & relB
      ## per pixelGroup and RelB per pixelGroup & speciesCode
      vegTypeTable[, sumB := sum(B), by = .(pixelGroup)]
      vegTypeTable[, relB := sum(B)/sumB, by = .(pixelGroup, speciesCode)]
      vegTypeTable[is.na(relB) & sumB == 0, relB := 0]
browser()
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

      vegTypesCD <- vegTypes[, list(vegType = as.factor(vegClass)), by = "pixelGroup"]

      vegTypesRas <- rasterizeReduced(reduced = vegTypesCD,
                                      fullRaster = sim$pixelGroupMap,
                                      mapcode = "pixelGroup", newRasterCols ="vegType")


      sim$vegTypes <- vegTypes

   sim <- scheduleEvent(sim, time(sim) + P(sim)$reclassTime, "vegReclass_WBI", "reclass")
    },

    plot = {
      yearSim <- paste0("Year", time(sim))
      sim$habitatMap[[yearSim]] <- ## TODO: BUILD A RASTERSTACK FUNCTION FOR ALL YEARS

   sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "vegReclass_WBI", "plot")
    },
  warning(paste("Undefined event type: \'", current(sim)[1, "eventType", with = FALSE],
                  "\' in module \'", current(sim)[1, "moduleName", with = FALSE], "\'", sep = ""))
  )
  return(invisible(sim))
}

## event functions
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization
Init <- function(sim) {

  # # ! ----- EDIT BELOW ----- ! #

  # ! ----- STOP EDITING ----- ! #

  return(invisible(sim))
}

### template for save events
Save <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  sim <- saveFiles(sim)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for plot events
plotFun <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  sampleData <- data.frame("TheSample" = sample(1:10, replace = TRUE))
  Plots(sampleData, fn = ggplotFn)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for your event1
Event1 <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # THE NEXT TWO LINES ARE FOR DUMMY UNIT TESTS; CHANGE OR DELETE THEM.
  # sim$event1Test1 <- " this is test for event 1. " # for dummy unit test
  # sim$event1Test2 <- 999 # for dummy unit test

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for your event2
Event2 <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # THE NEXT TWO LINES ARE FOR DUMMY UNIT TESTS; CHANGE OR DELETE THEM.
  # sim$event2Test1 <- " this is test for event 2. " # for dummy unit test
  # sim$event2Test2 <- 777  # for dummy unit test

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  # Any code written here will be run during the simInit for the purpose of creating
  # any objects required by this module and identified in the inputObjects element of defineModule.
  # This is useful if there is something required before simulation to produce the module
  # object dependencies, including such things as downloading default datasets, e.g.,
  # downloadData("LCC2005", modulePath(sim)).
  # Nothing should be created here that does not create a named object in inputObjects.
  # Any other initiation procedures should be put in "init" eventType of the doEvent function.
  # Note: the module developer can check if an object is 'suppliedElsewhere' to
  # selectively skip unnecessary steps because the user has provided those inputObjects in the
  # simInit call, or another module will supply or has supplied it. e.g.,
  # if (!suppliedElsewhere('defaultColor', sim)) {
  #   sim$map <- Cache(prepInputs, extractURL('map')) # download, extract, load file from url in sourceURL
  # }

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



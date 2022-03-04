#This is R script to process development study  data using XCMS version 3.12.0
#call package
library(xcms)

#setting to the number of CPUs available on your system
if (.Platform$OS.type == "unix") {
      register(bpstart(MulticoreParam(1)))
} else {
      register(bpstart(SnowParam(1)))
} 


#######################################pre-processing "Full scan" data, export MS1 features#######################################

MSdata <- paste0("full_scan/",list.files(path = "full_scan/",pattern = ".mzXML$", recursive = TRUE))
pd <- data.frame(sample_name = sub(basename(MSdata), 
                                   pattern = ".mzXML",
                                   replacement = "", 
                                   fixed = TRUE),
                 sample_group = c(rep("Against", 3),
                                  rep("Normal", 3),
                                  rep("Normal", 3),
                                  rep("Against", 3),
                                  rep("QC", 8)),
                 class=c(rep("Against", 3),
                         rep("Normal", 3),
                         rep("Normal", 3),
                         rep("Against", 3),
                         rep("QC", 8)),
                 stringsAsFactors
                 = FALSE)
rawData <- readMSData(MSdata, centroided. = TRUE, mode = "onDisk",
                      pdata = new("NAnnotatedDataFrame", pd))
rawData <- filterEmptySpectra(rawData)
cwp <- CentWaveParam(snthresh = 5, noise = 100, peakwidth = c(5, 30), ppm = 20)
processedData_raw <- findChromPeaks(rawData, param = cwp)
mpp <- MergeNeighboringPeaksParam(expandRt = 3, expandMz = 0, ppm=20)
xdata_pp <- refineChromPeaks(processedData_raw,param=mpp)
processedData <- xdata_pp
processedData$sample_type <- "study"
processedData$sample_type[c(1,2,3)] <- "QC"
pdp_subs <- PeakDensityParam(sampleGroups = processedData$sample_type,
                             minFraction = 0.9)
processedData <- groupChromPeaks(processedData, param = pdp_subs)
pgp_subs <- PeakGroupsParam(minFraction = 0.9,
                            subset = which(processedData$sample_type == "QC"),
                            subsetAdjust = "average", span = 0.2)
processedData <- adjustRtime(processedData, param = pgp_subs)
pdp <- PeakDensityParam(sampleGroups = processedData$sample_group,
                        minFraction = 1, binSize = 0.02)
processedData <- groupChromPeaks(processedData, param = pdp) 
medWidth <- median(chromPeaks(processedData)[, "rtmax"] -
                         chromPeaks(processedData)[, "rtmin"])
processed_Data <- fillChromPeaks(processedData, param = FillChromPeaksParam(fixedRt = medWidth))
featuresDef <- featureDefinitions(processed_Data)
featuresIntensities <- featureValues(processed_Data, value = "into")
dataTable <- merge(featuresDef, featuresIntensities, by=0, all=TRUE)
dataTable <- dataTable[, !(names(dataTable) %in% c("peakidx"))]
head(dataTable)
write.table(dataTable, "xcms_all_fullscan_MS1_features.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#######################################pre-processing "DDA" data, export MS1&MS2 features, MS2 features as MS/MS spectral library#######################################

MSdata <- paste0("DDA/",list.files(path = "DDA/",pattern = ".mzXML$", recursive = TRUE))
pd <- data.frame(sample_name = sub(basename(MSdata), 
                                   pattern = ".mzXML",
                                   replacement = "", 
                                   fixed = TRUE),
                 sample_group = c(rep("study", 15),
                                  rep("QC", 4)),
                 class=c(rep("study", 15),
                         rep("QC", 4)),
                 stringsAsFactors
                 = FALSE)
rawData <- readMSData(MSdata, centroided. = TRUE, mode = "onDisk",
                      pdata = new("NAnnotatedDataFrame", pd))
rawData <- filterEmptySpectra(rawData)
cwp <- CentWaveParam(snthresh = 5, noise = 100, peakwidth = c(5, 30), ppm = 20)
processedData_raw <- findChromPeaks(rawData, param = cwp)
mpp <- MergeNeighboringPeaksParam(expandRt = 3, expandMz = 0, ppm=20)
xdata_pp <- refineChromPeaks(processedData_raw,param=mpp)
processedData <- xdata_pp
processedData$sample_type <- "study"
processedData$sample_type[c(1,2,3)] <- "QC"
pdp_subs <- PeakDensityParam(sampleGroups = processedData$sample_type,
                             minFraction = 0.9)
processedData <- groupChromPeaks(processedData, param = pdp_subs)
pgp_subs <- PeakGroupsParam(minFraction = 0.9,
                            subset = which(processedData$sample_type == "QC"),
                            subsetAdjust = "average", span = 0.2)
processedData <- adjustRtime(processedData, param = pgp_subs)
pdp <- PeakDensityParam(sampleGroups = processedData$sample_group,
                        minFraction = 0.01, binSize = 0.02)
processedData <- groupChromPeaks(processedData, param = pdp) 
medWidth <- median(chromPeaks(processedData)[, "rtmax"] -
                         chromPeaks(processedData)[, "rtmin"])
source("https://raw.githubusercontent.com/jorainer/xcms-gnps-tools/master/customFunctions.R")
filteredMs2Spectra <- featureSpectra(processedData, return.type = c("MSpectra","list"))
filteredMs2Spectra <- clean(filteredMs2Spectra, all = TRUE)
filteredMs2Spectra <- formatSpectraForGNPS(filteredMs2Spectra)
writeMgfData(filteredMs2Spectra, "DDA_spectral_library.mgf")
processedData <- adjustRtime(processedData, param = pgp_subs)
pdp <- PeakDensityParam(sampleGroups = processedData$sample_group,
                        minFraction = 1, binSize = 0.02)
processedData <- groupChromPeaks(processedData, param = pdp) 
medWidth <- median(chromPeaks(processedData)[, "rtmax"] -
                         chromPeaks(processedData)[, "rtmin"])
featuresDef <- featureDefinitions(processedData)
featuresIntensities <- featureValues(processedData, value = "into")
dataTable <- merge(featuresDef, featuresIntensities, by=0, all=TRUE)
dataTable <- dataTable[, !(names(dataTable) %in% c("peakidx"))]
head(dataTable)
write.table(dataTable, "xcms_all_DDA_MS1_features.txt", sep = "\t", quote = FALSE, row.names = FALSE)


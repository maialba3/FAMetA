

################################################################################
## FAMetA: example script
##
## Maribel Alcoriza Balaguer: maribel_alcoriza@iislafe.es
################################################################################

library(FAMetA)
library(tools)

# Example files can be found at 
# https://drive.google.com/drive/folders/1-mhqBd6W8VJkIYwVuIOr2Gn4Z9vBUN9-?usp=sharing.
# This dataset contains 12 samples of linfocites grown in culture media 
# supplemented with 13Cglc, 13Cglc+iFASN, 13Cglc+iSCD and 13Cglc+iFADS2 
# (3 replicates of each condition), 2 blank samples and 3 injections of a 
# standards mix.  


# If any other external software is used for processing, data can be loaded from
# a csv file using the following function:
# fadata <- readfadatafile("externalfadata.csv", sep=",", dec=".")

# In this case, go directly to data correction step.

#==============================================================================#
# Data pre-processing using LipidMS package
#==============================================================================#

###################
# Read metadata
metadata <- read.csv("samples.csv", header = T, sep=",")

# check file names (they must include .mzXML)
if (!all(file_ext(metadata$sample) == "mzXML")){
  metadata$sample[file_ext(metadata$sample) != "mzXML"] <- 
    paste(metadata$sample[file_ext(metadata$sample) != "mzXML"], ".mzXML", sep="")
}


###################
# Set processing parameters

# Peak-picking
polarity <- "negative"
dmzagglom <- 15                       # dmz and drt to generate bins/partitions for peak-picking
drtagglom <- 200                      # max rt window for bins
drtclust <- 100                       # drt window for clustering (redefines previous bins)
minpeak <- 8                          # min number of points to define a peak (MS1, MS2)
minint <- 100000                      # at least minpeak points must have minint intensity
drtgap <- 5                           # max rt gap to fill missing points in a peak
drtminpeak <- 8                       # min width of a peak when there are more than 1 peak in a EIC
drtmaxpeak <- 30                      # max rt window for a EIC
recurs <- 10                          # max number of peaks in a EIC
sb <- 5                               # signal-to-baseline ratio (MS1, MS2)
sn <- 5                               # signal-to-noise ratio
weight <- 2                           # weight to assign new peaks
dmzIso <- 5                           # dmz for isotopes search
drtIso <- 5                           # drt for isotopes search

# Batch processing
dmzalign <- 10                        # max dmz and rt to group peaks for alignment
drtalign <- 60                        # max rt window for clustering in alignment
span <- 0.2                           # span for alignment
minsamplesfracalign <- 0.50           # min fraction of samples represented in a peak group to be used for alignment
dmzgroup <- 10                        # max dmz and rt to group peaks for grouping
drtagglomgroup <- 50                  # max rt window for clustering in grouping
drtgroup <- 10                        # max rt difference within a peak group
minsamplesfracgroup <- 0.20           # min fraction of samples represented in a peak group to be kept


parallel <- TRUE                      # parallel processing 
ncores <- 4                           # number of cores


###################
# Peak-picking
msbatch <- batchdataProcessing(metadata = metadata,
                               polarity = polarity,
                               dmzagglom = dmzagglom,
                               drtagglom = drtagglom,
                               drtclust = drtclust,
                               minpeak = minpeak,
                               drtgap = drtgap,
                               drtminpeak = drtminpeak,
                               drtmaxpeak = drtmaxpeak,
                               recurs = recurs,
                               sb = sb,
                               sn = sn,
                               minint = minint,
                               weight = weight,
                               dmzIso = dmzIso,
                               drtIso = drtIso,
                               parallel = parallel,
                               ncores = ncores)

ncores <- 2 
###################
# Alignment
msbatch <- alignmsbatch(msbatch, dmz = dmzalign, drt = drtalign, span = span,
                        minsamplesfrac = minsamplesfracalign,
                        parallel = parallel, ncores = ncores)

###################
# Grouping
msbatch <- groupmsbatch(msbatch, dmz = dmzgroup, drtagglom = drtagglomgroup,
                        drt = drtgroup, minsamplesfrac = minsamplesfracgroup,
                        parallel = parallel, ncores = ncores)

###################
# Save msbatch
save(msbatch, file="msbatch.rda.gz", compress = TRUE)
gc()


#==============================================================================#
# FA annotation
#==============================================================================#

# remove samples without 13C-tracer
msbatch$msobjects <- msbatch$msobjects[msbatch$metaData$sampletype != "mixFA"]
msbatch$metaData <- msbatch$metaData[msbatch$metaData$sampletype != "mixFA",]
msbatch$features <- msbatch$features[,!grepl("mixFA", colnames(msbatch$features))]


###################
# Annotate FA
msbatch <- annotateFA(msbatch, dmz = 5)

###################
# plot peaks from identified FAs to check them
plots <- plotFA(msbatch, dmz = 10)

pdf("FAs.pdf")
for (p in 1:length(plots)){
  print(plots[[p]])
}
dev.off()

###################
# export annotations for curation
write.csv(msbatch$fas, file="faid.csv", row.names=FALSE)

# FA annotations may be modified by removing rows of unwanted FA, modifying 
# initial and end retention times or adding new rows with missing compounds. 
# Here, unique compound names with the nomenclature "FA(16:1)n7", where n7 
# (omega-7) indicates the last double bound position, are required to differentiate 
# between FA isomers. In case of unknown positions, x, y and z letters are allowed 
# (i.e. FA(16:1)nx). In addition, internal standards for later normalization can 
# also be added at this point in anew row indicating IS in the compound name 
# column. Check faid.csv and faid_curated.csv files to see an example.
# Changes can also be performed directly on the msbatch using addFA, removeFA, 
# changeFArt and searchIS functions. See documentation for details. 

# Search IS within the msbatch features to add it to the annotations file
View(msbatch$features)

###################
# read csv file with modified annotations
faid <- read.csv("faid_curated.csv", sep=",", dec=".")

###################
# change FA annotations
msbatch <- curateFAannotations(msbatch, faid)

###################
# check new identifications
plots <- plotFA(msbatch, dmz = 10)

pdf("FAs_curated.pdf")
for (p in 1:length(plots)){
  print(plots[[p]])
}
dev.off()


#==============================================================================#
# Search FA isotopes and get the fadata object
#==============================================================================#
fadata <- searchFAisotopes(msbatch, dmz = 10, coelCutoff = 0.4)

save(fadata, file="fadata.rda")


# if you want to save fadata in a csv to subset it for example:
df <- cbind(rbind(fadata$fattyacids, data.frame(Compound="IS", Label="")), 
            rbind(fadata$intensities, fadata$IS))
df <- rbind(c("", "sampletype", fadata$metadata$sampletype),
            c("Compound", "Label", colnames(fadata$intensities)), df)
write.table(df, file="fadata.csv", sep=",", col.names = FALSE, row.names = FALSE)

# and then, you could read it using:
fadata <- readfadatafile("fadata.csv", sep=",", dec=".")

#==============================================================================#
# Data correction
#==============================================================================#
fadata <- dataCorrection(fadata, blankgroup = "Blank")

# Alternatively, to add external normalization:
# fadata <- dataCorrection(fadata, blankgroup = "blank", 
# externalnormalization = "protein")
# It requires a "protein" column at the fadata$metadata data frame.


#==============================================================================#
# Metabolic analysis
#==============================================================================#

###################
# Synthesis analysis
fadata <- synthesisAnalysis(fadata=fadata, R2Thr = 0.95, maxiter = 1e3,
                            maxconvergence = 100, startpoints = 5)

# If inhibitors have been used, make sure D2 has not been underestimated. If so,
# D2 could be set as the one calculated for 13-Glc Control samples to improve 
# the results:

# D2 <- fadata$synthesis$results$D2[fadata$synthesis$results$FA == "FA(16:0)"]
# fadata$synthesis$results$Group[fadata$synthesis$results$FA == "FA(16:0)"]
 
# D2[4:12] <- rep(mean(D2[1:3]))
 
# relaunch synthesis analysis fixing D2
# fadata <- synthesisAnalysis(fadata=fadata, R2Thr = 0.95, maxiter = 1e3,
#                             maxconvergence = 100, startpoints = 5, D2 = D2)


# Explore results
View(fadata$synthesis$results)
View(fadata$synthesis$predictedValues)
pdf("plotsSynthesis.pdf")
for (f in 1:length(fadata$synthesis$plots)){
  for (s in 1:length(fadata$synthesis$plots[[f]])){
    print(fadata$synthesis$plots[[f]][[s]])
  }
}
dev.off()

###################
# Elongation analysis
fadata <- elongationAnalysis(fadata, R2Thr = 0.95, maxiter = 1e4,
                             maxconvergence=100, startpoints = 5, DThr = 0.1)


# Explore results
View(fadata$elongation$results)
View(fadata$elongation$predictedValues)
pdf("plotsElongation.pdf")
for (f in 1:length(fadata$elongation$plots)){
  for (s in 1:length(fadata$elongation$plots[[f]])){
    print(fadata$elongation$plots[[f]][[s]])
  }
}
dev.off()

###################
# Desaturation analysis
fadata <- desaturationAnalysis(fadata)

# Explore results
View(fadata$desaturations$results)

###################
# Summarize results
fadata <- summarizeResults(fadata, controlgroup = "Control13Cglc")

###################
# Save fadata
save(fadata, file="fadata.rda")


###################
# Export results
write.csv(fadata$results$results, file = "results.csv", row.names=FALSE)
write.csv(fadata$results$summary, file = "summary.csv")
write.table(rbind(fadata$metadata$sampletype, 
                  colnames(fadata$mid),
                  fadata$mid),
            file = "mid.csv", sep=",", col.names = FALSE)
write.table(rbind(fadata$metadata$sampletype, 
                  colnames(fadata$synthesis$predictedValues),
                  fadata$synthesis$predictedValues),
            file = "predictedmid.csv", sep=",", col.names = FALSE)


pdf("relativepoolsizeRaw.pdf")
print(fadata$results$heatmaps$relativepoolsize$raw$plot)
dev.off()

write.table(rbind(fadata$metadata$sampletype, 
                  colnames(fadata$results$heatmaps$relativepoolsize$raw$values),
                  fadata$results$heatmaps$relativepoolsize$raw$values),
            file = "relativepoolsizeRaw.csv", sep=",", col.names = FALSE)


pdf("relativepoolsizeZscore.pdf")
print(fadata$results$heatmaps$relativepoolsize$zscore$plot)
dev.off()

write.table(rbind(fadata$metadata$sampletype, 
                  colnames(fadata$results$heatmaps$relativepoolsize$raw$values),
                  fadata$results$heatmaps$relativepoolsize$zscore$values),
            file = "relativepoolsizeZscore.csv", sep=",", col.names = FALSE)

if ("log2FC" %in% names(fadata$results$heatmaps$relativepoolsize)){
  pdf("relativepoolsizeLog2FC.pdf")
  print(fadata$results$heatmaps$relativepoolsize$log2FC$plot)
  dev.off()
  
  write.table(rbind(fadata$metadata$sampletype, 
                    colnames(fadata$results$heatmaps$relativepoolsize$log2FC$values),
                    fadata$results$heatmaps$relativepoolsize$log2FC$values),
              file = "relativepoolsizeLog2FC.csv", sep=",", col.names = FALSE)
}

pdf("resultsRaw_endogenouslySynthesized.pdf")
print(fadata$results$heatmaps$synthesized$raw$plot)
dev.off()

write.table(rbind(fadata$metadata$sampletype, 
                  colnames(fadata$results$heatmaps$synthesized$raw$values),
                  fadata$results$heatmaps$synthesized$raw$values),
            file = "resultsRaw_endogenouslySynthesized.csv", sep=",", 
            col.names = FALSE)


if ("log2FC" %in% names(fadata$results$heatmaps$synthesized)){
  pdf("resultsLog2FC_endogenouslySynthesized.pdf")
  print(fadata$results$heatmaps$synthesized$log2FC$plot)
  dev.off()
  
  write.table(rbind(fadata$metadata$sampletype, 
                    colnames(fadata$results$heatmaps$synthesized$log2FC$values),
                    fadata$results$heatmaps$synthesized$log2FC$values),
              file = "resultsLog2FC_endogenouslySynthesized.csv", sep=",", 
              col.names = FALSE)
}


# Isotopologue distributions: observed vs predicted

pdf("isotopologueDistributions.pdf")
for (f in 1:length(fadata$synthesis$plots)){
  for (s in 1:length(fadata$synthesis$plots[[f]])){
    print(fadata$synthesis$plots[[f]][[s]])
  }
}
for (f in 1:length(fadata$elongation$plots)){
  for (s in 1:length(fadata$elongation$plots[[f]])){
    print(fadata$elongation$plots[[f]][[s]])
  }
}
dev.off()



# Reorganized tables for synthesis and elongation parameters (S16, E1, E2, E3, 
# E4 and E5)

write.table(rbind(fadata$metadata$sampletype, 
                  colnames(fadata$results$allparameters$S16),
                  fadata$results$allparameters$S16),
            file = "S16.csv", sep=",", col.names = FALSE)

write.table(rbind(fadata$metadata$sampletype, 
                  colnames(fadata$results$allparameters$E1),
                  fadata$results$allparameters$E1),
            file = "E1.csv", sep=",", col.names = FALSE)

write.table(rbind(fadata$metadata$sampletype, 
                  colnames(fadata$results$allparameters$E2),
                  fadata$results$allparameters$E2),
            file = "E2.csv", sep=",", col.names = FALSE)

write.table(rbind(fadata$metadata$sampletype, 
                  colnames(fadata$results$allparameters$E3),
                  fadata$results$allparameters$E3),
            file = "E3.csv", sep=",", col.names = FALSE)

write.table(rbind(fadata$metadata$sampletype, 
                  colnames(fadata$results$allparameters$E4),
                  fadata$results$allparameters$E4),
            file = "E4.csv", sep=",", col.names = FALSE)

write.table(rbind(fadata$metadata$sampletype, 
                  colnames(fadata$results$allparameters$E5),
                  fadata$results$allparameters$E5),
            file = "E5.csv", sep=",", col.names = FALSE)

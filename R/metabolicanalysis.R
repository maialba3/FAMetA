# synthesisAnalysis
#' De novo synthesis analysis of fatty acids until 16 carbons.
#'
#' De novo synthesis analysis of fatty acids until 16 carbons.
#'
#' @param fadata fadata obtained from the msbatch with \link{searchFAisotopes}
#' function or read from csv file with \link{readfadatafile} function.
#' @param R2Thr positive numeric between 0 and 1 specifying the minimum R2
#' allowed for fits.
#' @param maxiter parameter passed to \link{nls.control}. Positive integer
#' specifying the maximum number of iterations allowed.
#' @param maxconvergence positive integer specifying the maximum number of
#' successes before choosing the winning model.
#' @param D1 positive numeric between 0 and 1 specifying the contribution of
#' acetate M+1. If NA it is estimated.
#' @param D2 positive numeric between 0 and 1 specifying the contribution of
#' acetate M+2. If NA it is estimated.
#' @param P overdispersion parameter. If NA it is estimated (quasi-multinomial
#' distribution). If set to 0, no overdispersion is assumed (multinomial
#' distribution).
#' @param startpoints positive integer specifying the number of starting points
#' for each parameter to be estimated.
#' @param parameters parameters to be estimated for each fatty acid. It can be
#' modified to change them or to add new fatty acids.
#' @param propagateD logical. If TRUE, unsaturated fatty acids use estimated D0,
#' D1,D2 and P values for saturated fatty acids (14:0 for FA shorter than 16C
#' and 16:0 for FA with 16C.).
#' @param verbose print information messages.
#'
#' @details Synthesis analysis will model FA data for FA up to 16 carbons to
#' estimate 13C-tracer contribution to the acetyl-CoA pool for FA synthesis (D)
#' and the FA fraction that has been synthesized de novo. D0, D1 and D2
#' represent the contribution of M+0, M+1 and M+2 acetate, respectively, and P
#' (phi) is the overdispersion parameter of the quasi-multinomial distribution.
#' D0, D1, D2 can also be fixed if they are known. This is particularly useful
#' in case inhibitors have been used as they could reduce S below the confidence
#' interval and thus, S and D parameters could be misestimated.
#'
#' @return fadata list. Synthesis analysis results will be saved at the 
#' synthesis element of the fa list.
#'
#' @examples
#' \donttest{
#' ssdata <- dataCorrection(ssexamplefadata, blankgroup="Blank")
#' ssdata <- synthesisAnalysis(ssdata, R2Thr = 0.95, maxiter = 1e3,
#' maxconvergence = 100, startpoints = 5)
#' }
#' 
#' \dontrun{
#' fadata <- dataCorrection(examplefadata, blankgroup = "Blank")
#' fadata <- synthesisAnalysis(fadata, R2Thr = 0.95, maxiter = 1e3,
#' maxconvergence = 100, startpoints = 5)
#'
#' # If inhibitors have been used, make sure D2 has not been underestimated. If so,
#' # D2 could be set as the one calculated for 13-Glc Control samples to improve
#' # the results:
#'
#' # D2 <- fadata$synthesis$results$D2[fadata$synthesis$results$FA == "FA(16:0)"]
#' # fadata$synthesis$results$Group[fadata$synthesis$results$FA == "FA(16:0)"]
#'
#' # D2[4:12] <- rep(mean(D2[1:3]))
#'
#' # relaunch synthesis analysis fixing D2
#' # fadata <- synthesisAnalysis(fadata, R2Thr = 0.95, maxiter = 1e3,
#' #                             maxconvergence = 100, startpoints = 5, D2 = D2)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
synthesisAnalysis <- function(fadata,
                        R2Thr = 0.98,
                        maxiter = 1e3,
                        maxconvergence = 100,
                        D1 = NA,
                        D2 = NA,
                        P = NA,
                        startpoints = 5,
                        parameters = FAMetA::parameters,
                        propagateD = TRUE, 
                        verbose = TRUE){

  #============================================================================#
  # Check arguments
  #============================================================================#
  if (!is.list(fadata) | length(fadata) < 6 |
      !all(c("metadata", "fattyacids", "intensities", "mid", "poolsize",
             "relativepoolsize") %in% names(fadata))){
    stop("fadata must be a list with at least 5 elements: metadata, fattyacids
         (data.frame with two columns: Compound and Label), intensities, mid,
         poolsize and relativepoolsize")
  }
  if (!length(D1) %in% c(1, length(fadata$metadata$sample))){
    stop("D1 must have a length of 1 or equal to the number of samples")
  }
  if (!length(D2) %in% c(1, length(fadata$metadata$sample))){
    stop("D2 must have a length of 1 or equal to the number of samples")
  }
  if (!length(P) %in% c(1, length(fadata$metadata$sample))){
    stop("P must have a length of 1 or equal to the number of samples")
  }
  
  if (length(D1) == 1){inputD1 <- rep(D1, nrow(fadata$metadata))} else {inputD1 <- D1}
  if (length(D2) == 1){inputD2 <- rep(D2, nrow(fadata$metadata))} else {inputD2 <- D2}
  if (length(P) == 1){inputP <- rep(P, nrow(fadata$metadata))} else {inputP <- P}
  #============================================================================#
  # Run DNSIA algorithm for each compound and sample
  #============================================================================#
  compounds <- fadata$fattyacids$Compound
  toDo <- parameters[which(parameters$FattyAcid %in% compounds &
                                     parameters[, "M"]  <= 16), "FattyAcid"]

  if (length(toDo) > 0){
    toDo1 <- toDo[grep("FA\\(.*:0\\)", toDo)]
    toDo2 <- toDo[grep("FA\\(.*:0\\)", toDo, invert = TRUE)]


    if (length(toDo1) > 0){
      results1 <- runSynthesisAnalysis(fadata = fadata, toDo = toDo1,
                                       R2Thr = R2Thr, maxiter = maxiter,
                                       maxconvergence = maxconvergence,
                                       D1 = inputD1, D2 = inputD2, P = inputP,
                                       startpoints = startpoints,
                                       verbose = verbose)
      if (length(toDo2) > 0){
        toDo14 <- toDo2[grep("FA\\(16", toDo2, invert = TRUE)]
        if (length(toDo14) > 0){
          if (propagateD){
            if ("FA(14:0)" %in% toDo1){
              D1 <- results1$results$D1[results1$results$FA == "FA(14:0)"]
              D2 <- results1$results$D2[results1$results$FA == "FA(14:0)"]
              P <- results1$results$P[results1$results$FA == "FA(14:0)"]
            } else {
              D1 <- results1$results$D1[results1$results$FA == "FA(16:0)"]
              D2 <- results1$results$D2[results1$results$FA == "FA(16:0)"]
              P <- results1$results$P[results1$results$FA == "FA(16:0)"]*8/7
            }
          } else {
            D1 <- inputD1
            D2 <- inputD2
            P <- inputP
          }
          synt_14 <- runSynthesisAnalysis(fadata = fadata, toDo = toDo14,
                                          R2Thr = R2Thr, maxiter = maxiter,
                                          maxconvergence = maxconvergence,
                                          D1 = D1, D2 = D2, P = P,
                                          startpoints = startpoints,
                                          verbose = verbose)
        } else {
          synt_14 <- list()
        }

        toDo16 <- toDo2[grep("FA\\(16", toDo2)]
        if (length(toDo16) > 0){
          if (propagateD){
            if ("FA(16:0)" %in% toDo1){
              D1 <- results1$results$D1[results1$results$FA == "FA(16:0)"]
              D2 <- results1$results$D2[results1$results$FA == "FA(16:0)"]
              P <- results1$results$P[results1$results$FA == "FA(16:0)"]
            } else {
              D1 <- results1$results$D1[results1$results$FA == "FA(14:0)"]
              D2 <- results1$results$D2[results1$results$FA == "FA(14:0)"]
              P <- results1$results$P[results1$results$FA == "FA(14:0)"]*7/8
            }
          } else {
            D1 <- inputD1
            D2 <- inputD2
            P <- inputP
          }
          synt_16 <- runSynthesisAnalysis(fadata = fadata, toDo = toDo16,
                                         R2Thr = R2Thr, maxiter = maxiter,
                                         maxconvergence = maxconvergence,
                                         D1 = D1, D2 = D2, P = P,
                                         startpoints = startpoints,
                                         parameters = parameters,
                                         verbose = verbose)
        } else {
          synt_16 <- list()
        }

        if (length(synt_14) > 0 & length(synt_16) > 0){
          results <- list()
          results$fas <- c(results1$fas, synt_14$fas, synt_16$fas)
          results$results <- data.frame(rbind(results1$results,
                                              synt_14$results,
                                              synt_16$results))
          results$predictedValues <- data.frame(rbind(results1$predictedValues,
                                                      synt_14$predictedValues,
                                                      synt_16$predictedValues))
          results$plots <- c(results1$plots, synt_14$plots, synt_16$plots)
          results$details <- c(results1$details, synt_14$details, synt_16$details)
        } else {
          if (length(synt_14) > 0){
            results <- list()
            results$fas <- c(results1$fas, synt_14$fas)
            results$results <- data.frame(rbind(results1$results,
                                                synt_14$results))
            results$predictedValues <- data.frame(rbind(results1$predictedValues,
                                                        synt_14$predictedValues))
            results$plots <- c(results1$plots, synt_14$plots)
            results$details <- c(results1$details, synt_14$details)
          } else {
            results <- list()
            results$fas <- c(results1$fas, synt_16$fas)
            results$results <- data.frame(rbind(results1$results,
                                                synt_16$results))
            results$predictedValues <- data.frame(rbind(results1$predictedValues,
                                                    synt_16$predictedValues))
            results$plots <- c(results1$plots, synt_16$plots)

            results$details <- c(results1$details, synt_16$details)
          }
        }
      } else {
        results <- results1
      }
    } else {
      warning("At least one of myristic or palmitic acid should be included in the data set.")
      results <- runSynthesisAnalysis(fadata = fadata, toDo = toDo2,
                                      R2Thr = R2Thr, maxiter = maxiter,
                                      maxconvergence = maxconvergence,
                                      D1 = inputD1, D2 = inputD2, P = inputP, 
                                      startpoints = startpoints)
    }
    fadata$synthesis <- results
  } else {
    stop("No FA available for synthesis analysis.")
  }

  return(fadata)
}

# elongationAnalysis
#' Elongation analysis of fatty acids longer than 16 carbons.
#'
#' Elongation analysis of fatty acids longer than 16 carbons.
#'
#' @param fadata fadata containing synthesis results.
#' @param R2Thr positive numeric between 0 and 1 specifying the minimum R2
#' allowed for fits.
#' @param maxiter parameter passed to \link{nls.control}. Positive integer
#' specifying the maximum number of iterations allowed.
#' @param maxconvergence positive integer specifying the maximum number of
#' successes before choosing the winning model.
#' @param startpoints positive integer specifying the number of starting points
#' for each parameter to be estimated.
#' @param D2Thr minimum D2 value allowed to perform the elongation analysis.
#' @param parameters parameters to be estimated for each fatty acid. It can be
#' modified to change them or to add new fatty acids (adding new rows).
#' @param verbose print information messages.
#'
#' @details Main route of de novo synthesis plus elongation starts at 16 carbons
#' and then adds blocks of 2 carbons. Therefore, isotopologue distributions for
#' FA longer than 16 carbons will be modeled taking into account de novo
#' synthesis until FA(16:0), followed by single and independent elongation steps
#' (E1, E2 …, En). Parameters D0, D1 and D2 are imported from FA(16:0) or FA(14:0)
#' and thus, the only relevant parameters to be estimated in the elongation
#' analysis are Ei and I. For n6 and n3 series, elongation is expected from
#' FA(18:2)n6 and FA(18:3)n3 so that synthesis (S16:0) and first elongation step
#' (E1) are set to 0.
#'
#' @return fadata list. Elongation analysis results will be saved at the 
#' elongation element of the fa list.
#'
#' @examples
#' \donttest{
#' ssdata <- dataCorrection(ssexamplefadata, blankgroup="Blank")
#' ssdata <- synthesisAnalysis(ssdata, R2Thr = 0.95, maxiter = 1e3,
#' maxconvergence = 100, startpoints = 5)
#' ssdata <- elongationAnalysis(ssdata, R2Thr = 0.95, maxiter = 1e4,
#' maxconvergence=100, startpoints = 5, D2Thr = 0.1)
#' }
#' 
#' \dontrun{
#' fadata <- dataCorrection(examplefadata, blankgroup = "Blank")
#' fadata <- synthesisAnalysis(fadata, R2Thr = 0.95, maxiter = 1e3,
#' maxconvergence = 100, startpoints = 5)
#' fadata <- elongationAnalysis(fadata, R2Thr = 0.95, maxiter = 1e4,
#' maxconvergence=100, startpoints = 5, D2Thr = 0.1)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
elongationAnalysis <- function(fadata, R2Thr = 0.98, maxiter = 1e4,
                         maxconvergence=100, startpoints=5, D2Thr = 0.1,
                         parameters = FAMetA::parameters, verbose = TRUE){

  #============================================================================#
  # Check arguments
  #============================================================================#
  if (!is.list(fadata) | length(fadata) < 7 |
      !all(c("metadata", "fattyacids", "intensities", "mid", "poolsize",
             "relativepoolsize", "synthesis") %in% names(fadata))){
    stop("fadata must be a list with at least 5 elements: metadata, fattyacids
         (data.frame with two columns: Compound and Label), intensities, mid and
         poolsize, relativepoolsize and synthesis")
  }

  if (missing(parameters)){
    parameters = FAMetA::parameters
  }

  #============================================================================#
  # Run EIA algorithm for each compound and sample
  #============================================================================#
  fas <- data.frame(cbind(fadata$fattyacids[,c("Compound", "Label")],
                          fadata$mid), stringsAsFactors = FALSE)
  fas <- split(fas, fas$Compound)

  compounds <- names(fas)
  toDo <- parameters[which(parameters$FattyAcid %in% compounds &
                                      parameters[, "M"]  > 16), "FattyAcid"]

  # get D and P parameters from synthesis analysis results
  # preference order: 16:0, 14:0, any 16:X, any 14:X
  if ("FA(16:0)" %in% fadata$synthesis$fas){
    D1 <- fadata$synthesis$results$D1[fadata$synthesis$results$FA == "FA(16:0)"]
    D2 <- fadata$synthesis$results$D2[fadata$synthesis$results$FA == "FA(16:0)"]
    P <- fadata$synthesis$results$P[fadata$synthesis$results$FA == "FA(16:0)"]
  } else {
    D1 <- D2 <- P <- rep(NA, ncol(fadata$mid))
  }
  if (any(is.na(D2)) & "FA(14:0)" %in% fadata$synthesis$fas){
    D1[is.na(D1)] <- fadata$synthesis$results$D1[fadata$synthesis$results$FA == "FA(14:0)"][is.na(D1)]
    D2[is.na(D2)] <- fadata$synthesis$results$D2[fadata$synthesis$results$FA == "FA(14:0)"][is.na(D2)]
    P[is.na(P)] <- fadata$synthesis$results$P[fadata$synthesis$results$FA == "FA(14:0)"][is.na(P)]*7/8
  }
  if (any(is.na(D2)) & any(grepl("FA\\(16", fadata$synthesis$fas))){
    index <- fadata$synthesis$results$Sample[grepl("FA\\(16", fadata$synthesis$results$FA)]
    means16 <- apply(fadata$synthesis$results[grepl("FA\\(16", fadata$synthesis$results$FA),4:10],
                     2, function(x) tapply(x, index, mean, na.rm=TRUE))
    D1[is.na(D1)] <- means16[is.na(D1),"D1"]
    D2[is.na(D2)] <- means16[is.na(D2),"D2"]
    P[is.na(P)] <- means16[is.na(P),"P"]
  }
  if (any(is.na(D2)) & any(grepl("FA\\(14", fadata$synthesis$fas))){
    index <- fadata$synthesis$results$Sample[grepl("FA\\(14", fadata$synthesis$results$FA)]
    means16 <- apply(fadata$synthesis$results[grepl("FA\\(14", fadata$synthesis$results$FA),4:10],
                     2, function(x) tapply(x, index, mean, na.rm=TRUE))
    D1[is.na(D1)] <- means16[is.na(D1),"D1"]
    D2[is.na(D2)] <- means16[is.na(D2),"D2"]
    P[is.na(P)] <- means16[is.na(P),"P"]*7/8
  }
  
  if (length(toDo) > 0){
    resultsElong <- list()
    for (c in toDo){
      if(verbose){cat(paste("\n", c, "...", sep=""))}
      fa <- fas[[c]]
      fa <- fa[order(fa$Label),, drop = FALSE]
      M <- as.numeric(parameters[parameters$FattyAcid == c, "M"])
      S16 <- as.numeric(parameters[parameters$FattyAcid == c, "S16"])
      E1 <- as.numeric(parameters[parameters$FattyAcid == c, "E1"])
      E2 <- as.numeric(parameters[parameters$FattyAcid == c, "E2"])
      E3 <- as.numeric(parameters[parameters$FattyAcid == c, "E3"])
      E4 <- as.numeric(parameters[parameters$FattyAcid == c, "E4"])
      E5 <- as.numeric(parameters[parameters$FattyAcid == c, "E5"])
      if (S16 == 1){S16 <- NA}else{S16 <- 0}
      if (E1 == 1){E1 <- NA}else{E1 <- 0}
      if (E2 == 1){E2 <- NA}else{E2 <- 0}
      if (E3 == 1){E3 <- NA}else{E3 <- 0}
      if (E4 == 1){E4 <- NA}else{E4 <- 0}
      if (E5 == 1){E5 <- NA}else{E5 <- 0}
      resultsElong[[c]] <- runElongationAnalysis(fa = fa, M = M, D1 = D1,
                                                 D2 = D2, P = P, S16 = S16,
                                                 E1 = E1, E2 = E2, E3 = E3,
                                                 E4 = E4, E5 = E5, R2Thr = R2Thr,
                                                 maxiter = maxiter,
                                                 maxconvergence = maxconvergence,
                                                 startpoints = startpoints,
                                                 D2Thr = D2Thr)
      if(verbose){cat("OK")}
    }
    #==========================================================================#
    # Summarize results
    #==========================================================================#
    results <- cbind(FA = rep(toDo, each = ncol(fadata$intensities)),
                            Sample = colnames(fadata$intensities),
                            Group = fadata$metadata$sampletype,
                            do.call(rbind, lapply(resultsElong, function(x) x[[1]])))

    predictedValues <- data.frame(do.call(rbind, lapply(resultsElong, function(x)
      do.call(cbind, lapply(x$models, function(y) if (!any(is.na(y))){predict(y)}else{rep(NA, nrow(x$mid))})))))
    compounds <- do.call(rbind, lapply(resultsElong, function(x) x[[2]][,c("Compound", "Label")]))
    predictedValues <- cbind(compounds, predictedValues)
    rownames(predictedValues) <- paste(predictedValues$Compound, "_M+", 
                                       predictedValues$Label, sep="")

  } else {
    stop("Unable to run EIA")
  }

  #============================================================================#
  # Plot observed vs predicted distributions based on the estimated parameters
  #============================================================================#
  plots <- plotResults(resultsElong, fadata$metadata$sampletype, toDo)

  fadata$elongation <- list(fas = toDo, results = results,
                            predictedValues = predictedValues, plots = plots,
                            details = resultsElong)

  return(fadata)
}

# desaturationAnalysis
#' Desaturation analysis of fatty acids.
#'
#' Desaturation analysis of fatty acids.
#'
#' @param fadata fadata containing synthesis and elongation results.
#' @param desaturationsdb desaturation reactions considered. It can be
#' modified to change them or to add new reactions.
#' @param SEThr minimum S or E value allowed to perform estimate desaturations.
#'
#' @details Once synthesis and elongation parameters have been estimated, these
#' results can be used to calculate the FA fraction that comes from desaturation
#' in unsaturated FA. For a given unsaturated FA (e.g. FA(18:1n9) we can
#' conceptually consider a one-step elongation-desaturation reaction (in this
#' example directly from FA(16:0) to FA(18:1n9) (E1') or a two-step elongation
#' followed by desaturation process (in this example FA(16:0) is elongated to
#' FA(18:0) (E1) and then desaturated to FA(18:1n9) (Des). Therefore,
#' desaturation can be estimated based on the fraction of E1', which is E1 from
#' FA(18:1)n9, and E1, which is E1 from FA(18:0). This same model can be used
#' for all known desaturation steps (see FAMetA::desaturationsdb) as long as
#' precursor and product FA isomers have been correctly and uniquely identified
#' and stationary state has been reached.
#'
#' @return fadata list. Desaturation analysis results will be saved at the 
#' desaturation element of the fa list.
#'
#' @examples
#' \donttest{
#' ssdata <- dataCorrection(ssexamplefadata, blankgroup="Blank")
#' ssdata <- synthesisAnalysis(ssdata, R2Thr = 0.95, maxiter = 1e3,
#' maxconvergence = 100, startpoints = 5)
#' ssdata <- elongationAnalysis(ssdata, R2Thr = 0.95, maxiter = 1e4,
#' maxconvergence=100, startpoints = 5, D2Thr = 0.1)
#' ssdata <- desaturationAnalysis(ssdata, SEThr = 0.05)
#' }
#' 
#' \dontrun{
#' fadata <- dataCorrection(examplefadata, blankgroup = "Blank")
#' fadata <- synthesisAnalysis(fadata, R2Thr = 0.95, maxiter = 1e3,
#' maxconvergence = 100, startpoints = 5)
#' fadata <- elongationAnalysis(fadata, R2Thr = 0.95, maxiter = 1e4,
#' maxconvergence=100, startpoints = 5, D2Thr = 0.1)
#' fadata <- desaturationAnalysis(fadata, SEThr = 0.05)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
desaturationAnalysis <- function(fadata,
                                 desaturationsdb = FAMetA::desaturationsdb,
                                 SEThr = 0.05){

  #============================================================================#
  # Check arguments
  #============================================================================#
  if (!is.list(fadata) | length(fadata) < 7 |
      !all(c("metadata", "fattyacids", "intensities", "mid", "poolsize",
             "relativepoolsize", "synthesis") %in% names(fadata))){
    stop("fadata must be a list with at least 5 elements: metadata, fattyacids
         (data.frame with two columns: Compound and Label), intensities, mid and
         poolsize, relativepoolsize and synthesis")
  }

  #============================================================================#
  # Run desaturation analysis for synthesized and elongated fatty acids
  #============================================================================#
  desaturations <- c()
  for (d in 1:nrow(desaturationsdb)){
    prec <- desaturationsdb[d,1]
    prod <- desaturationsdb[d,2]
    param <- desaturationsdb[d,3]

    if (param == "S"){
      if (prec %in% fadata$synthesis$fas &
          prod %in% fadata$synthesis$fas){
        desat <- fadata$synthesis$results[fadata$synthesis$results$FA == prod,"S"]/
          fadata$synthesis$results[fadata$synthesis$results$FA == prec,"S"]
        observ <- rep("", length(desat))
        observ[fadata$synthesis$results[fadata$synthesis$results$FA == prod,"S"] < SEThr |
                 fadata$synthesis$results[fadata$synthesis$results$FA == prec,"S"] < SEThr] <-
          "Warning: One S parameter out of confidence interval (SEThr).
        Desaturation level could be innacurate."
        observ[fadata$synthesis$results[fadata$synthesis$results$FA == prod,"S"] < SEThr &
                 fadata$synthesis$results[fadata$synthesis$results$FA == prec,"S"] < SEThr] <-
          "Warning: Both S parameters out of confidence interval (SEThr).
        Desaturation level could not be estimated."
        desat[fadata$synthesis$results[fadata$synthesis$results$FA == prod,"S"] < SEThr &
                fadata$synthesis$results[fadata$synthesis$results$FA == prec,"S"] < SEThr] <- NA
        desat <- data.frame(FA = prod, Sample = colnames(fadata$intensities),
                            Group = fadata$metadata$sampletype, Des = desat, Observations = observ)
        desaturations <- data.frame(rbind(desaturations, desat))
      }
    } else if (param %in% c("E1", "E2", "E3", "E4", "E5")){
      if (prec %in% fadata$elongation$fas &
          prod %in% fadata$elongation$fas){
        desat <- fadata$elongation$results[fadata$elongation$results$FA == prod,param]/
          fadata$elongation$results[fadata$elongation$results$FA == prec,param]
        observ <- rep("", length(desat))
        observ[fadata$elongation$results[fadata$elongation$results$FA == prod,param] < SEThr |
                 fadata$elongation$results[fadata$elongation$results$FA == prec,param] < SEThr] <-
          "Warning: One E parameter out of confidence interval (SEThr).
        Desaturation level could be innacurate."
        observ[fadata$elongation$results[fadata$elongation$results$FA == prod,param] < SEThr &
                 fadata$elongation$results[fadata$elongation$results$FA == prec,param] < SEThr] <-
          "Warning: Both E parameters out of confidence interval (SEThr).
        Desaturation level could not be estimated."
        desat[fadata$elongation$results[fadata$elongation$results$FA == prod,param] < SEThr &
                fadata$elongation$results[fadata$elongation$results$FA == prec,param] < SEThr] <- NA
        desat <- data.frame(FA = prod, Sample = colnames(fadata$intensities),
                            Group = fadata$metadata$sampletype, Des = desat, Observations = observ)
        desaturations <- data.frame(rbind(desaturations, desat))
      }
    }
  }
  if (length(desaturations) > 0){
    desaturations$Des <- as.numeric(desaturations$Des)
    desaturations$I <- as.numeric(1 - desaturations$Des)
    if (any(desaturations$Des[!is.na(desaturations$Des)]  > 1)){
      desaturations$Observations[desaturations$Des > 1] <-
        "Warning: desaturation parameter > 1, make sure your samples have reach
    the steady state and check isomers annotations."
      desaturations$I[desaturations$Des > 1] <- 0
      desaturations$Des[desaturations$Des > 1] <- 1
    }
    desaturations <- desaturations[,c("FA", "Sample", "Group", "Des", "I", "Observations")]
    fas <- unique(desaturations$FA)
  } else {
    fas <- vector()
    desaturations <- data.frame()
  }
  fadata$desaturation <- list(fas = fas, results = desaturations)
  return(fadata)
}

# summarizeResults
#' Obtain result tables and heatmaps that help interpreting your results.
#'
#' Obtain result tables and heatmaps that help interpreting your results.
#'
#' @param fadata fadata containing synthesis, elongation and desaturation results.
#' @param controlgroup name of the control group to compare the results.
#' @param parameters parameters to be estimated for each fatty acid. It can be
#' modified to change them or to add new fatty acids.
#'
#' @return fadata list with a results element which contains: results data frame
#' (results for the main parameters for each fatty acid and sample), summary 
#' data frame (mean and sd by sample groups for each parameter and fatty acids 
#' from the results table), different heatmaps representing pool size and 
#' results (values represented are also saved in data frames) and tables 
#' summarizing all parameters values for synthesis and elongation (S16, E1, E2, 
#' E3, E4 and E5).
#'
#' @examples
#' \donttest{
#' ssdata <- dataCorrection(ssexamplefadata, blankgroup="Blank")
#' ssdata <- synthesisAnalysis(ssdata, R2Thr = 0.95, maxiter = 1e3,
#' maxconvergence = 100, startpoints = 5)
#' ssdata <- elongationAnalysis(ssdata, R2Thr = 0.95, maxiter = 1e4,
#' maxconvergence=100, startpoints = 5, D2Thr = 0.1)
#' ssdata <- desaturationAnalysis(ssdata, SEThr = 0.05)
#' ssdata <- summarizeResults(ssdata)
#' }
#' 
#' \dontrun{
#' fadata <- dataCorrection(examplefadata, blankgroup = "Blank")
#' fadata <- synthesisAnalysis(fadata, R2Thr = 0.95, maxiter = 1e3,
#' maxconvergence = 100, startpoints = 5)
#' fadata <- elongationAnalysis(fadata, R2Thr = 0.95, maxiter = 1e4,
#' maxconvergence=100, startpoints = 5, D2Thr = 0.1)
#' fadata <- desaturationAnalysis(fadata, SEThr = 0.05)
#' fadata <- summarizeResults(fadata, controlgroup = "Control13Cglc")
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
summarizeResults <- function(fadata, controlgroup = NA,
                             parameters = FAMetA::parameters){

  #============================================================================#
  # Check arguments
  #============================================================================#
  #============================================================================#
  # Check arguments
  #============================================================================#
  if (!is.list(fadata) | length(fadata) < 7 |
      !all(c("metadata", "fattyacids", "intensities", "mid", "poolsize",
             "relativepoolsize", "synthesis") %in% names(fadata))){
    stop("fadata must be a list with at least 5 elements: metadata, fattyacids
         (data.frame with two columns: Compound and Label), intensities, mid and
         poolsize, relativepoolsize and synthesis")
  }
  if (missing(parameters)){
    parameters = FAMetA::parameters
  }

  results <- list()
  #============================================================================#
  # Create results table
  #============================================================================#
  results$results <- getResultsTable(fadata, parameters = parameters)

  #============================================================================#
  # Create summary table
  #============================================================================#
  results$summary <- getSummaryTable(fadata, resultstable = results$results,
                                     parameters = parameters)

  #============================================================================#
  # Create graphical output: pool size heatmap
  #============================================================================#
  palette <- c("#42858C", "#FE9300", "#870E75", "#3E71A8",
               "#7F8E39", "#5F3659", "#E5C616",
               "#16A08CFF", "#FE6900", "#628395", "#C5D86D", "#969696FF",
               "#358359FF", "#9F4D23FF", "#D86C4FFF", "#170C2EFF",
               "#473B75FF", "#F19C1FFF",
               "#117733", "#DDCC77", "#CC6677", "#88CCEE",
               "#44AA99", "#332288", "#AA4499", "#999933",
               "#882255", "#661100", "#6699CC", "#888888")
  if (length(unique(fadata$metadata$sampletype)) > length(palette)){
    set.seed(19580811)
    colors <- grDevices::colors()[grep('gr(a|e)y|white|light', grDevices::colors(), invert = T)]
    colors <- sample(colors, size = (length(unique(fadata$metadata$sampletype)) - length(palette)))
    palette <- c(palette, colors)
  }
  cols <- palette[as.factor(fadata$metadata$sampletype)]

  #============================================================================#
  # pool size heatmaps
  #============================================================================#
  results$heatmaps <- list()
  results$heatmaps$poolsize <- list()
  toPlot <- fadata$poolsize

  results$heatmaps$poolsize$raw <- list()
  results$heatmaps$poolsize$raw$plot <- heatmapResults(toPlot, cols, scale = "none")
  results$heatmaps$poolsize$raw$values <- toPlot

  results$heatmaps$poolsize$zscore <- list()
  results$heatmaps$poolsize$zscore$plot <- heatmapResults(toPlot, cols, scale = "row")
  if (ncol(toPlot) > 1){
    results$heatmaps$poolsize$zscore$values <- t(apply(toPlot, 1, scale, center=TRUE))
  } else {
    results$heatmaps$poolsize$zscore$values <- data.frame()
  }

  if (controlgroup %in% fadata$metadata$sampletype){
    toPlot <- log2(t(apply(fadata$poolsize, 1, function(x){
      x/mean(x[toupper(fadata$metadata$sampletype) == toupper(controlgroup)],
             na.rm = TRUE)
      })))

    results$heatmaps$poolsize$log2FC <- list()
    results$heatmaps$poolsize$log2FC$plot <- heatmapResults(toPlot, cols, scale = "none",
                                                       breaks = seq(-2 , 2, length=20))
    results$heatmaps$poolsize$log2FC$values <- toPlot
  }


  #============================================================================#
  # relative pool size heatmaps
  #============================================================================#
  results$heatmaps$relativepoolsize <- list()
  toPlot <- fadata$relativepoolsize

  results$heatmaps$relativepoolsize$raw <- list()
  results$heatmaps$relativepoolsize$raw$plot <- heatmapResults(toPlot, cols, scale = "none")
  results$heatmaps$relativepoolsize$raw$values <- toPlot

  results$heatmaps$relativepoolsize$zscore <- list()
  results$heatmaps$relativepoolsize$zscore$plot <- heatmapResults(toPlot, cols, scale = "row")
  if (ncol(toPlot) > 1){
    results$heatmaps$relativepoolsize$zscore$values <- t(apply(toPlot, 1, scale, center=TRUE))
  } else {
    results$heatmaps$relativepoolsize$zscore$values <- data.frame()
  }
  

  if (controlgroup %in% fadata$metadata$sampletype){
    toPlot <- log2(t(apply(fadata$relativepoolsize, 1, function(x){
      x/mean(x[toupper(fadata$metadata$sampletype) == toupper(controlgroup)], na.rm = TRUE)
    })))

    results$heatmaps$relativepoolsize$log2FC <- list()
    results$heatmaps$relativepoolsize$log2FC$plot <- heatmapResults(toPlot, cols, scale = "none",
                                                       breaks = seq(-2 , 2, length=20))
    results$heatmaps$relativepoolsize$log2FC$values <- toPlot
  }

  #============================================================================#
  # endogenously synthesized
  #============================================================================#
  results$heatmaps$synthesized <- list()
  toPlot <- data.frame(Sample = results$results$Sample, FA = results$results$FA,
                       S = apply(results$results[,c("S", "E", "Des")], 1, function(x){
                         if (!all(is.na(x))){max(x, na.rm=TRUE)} else {NA}
                       }))
  toPlot <- as.data.frame(pivot_wider(toPlot, names_from = "Sample",
                                      values_from = "S"))
  rownames(toPlot) <- toPlot$FA
  toPlot <- toPlot[,colnames(toPlot) != "FA", drop = FALSE]

  results$heatmaps$synthesized$raw <- list()
  results$heatmaps$synthesized$raw$plot <- heatmapResults(toPlot, cols, scale = "none",
                                                     breaks = seq(0 , 1, length=20))
  results$heatmaps$synthesized$raw$values <- toPlot

  if (controlgroup %in% fadata$metadata$sampletype){
    toPlot <- log2(t(apply(toPlot, 1, function(x){
      x/mean(x[toupper(fadata$metadata$sampletype) == toupper(controlgroup)],
             na.rm = TRUE)
    })))

    results$heatmaps$synthesized$log2FC <- list()
    results$heatmaps$synthesized$log2FC$plot <- heatmapResults(toPlot, cols, scale = "none",
                                                          breaks = seq(-2 , 2, length=20))
    results$heatmaps$synthesized$log2FC$values <- toPlot
  }
  
  #============================================================================#
  # intermediate parameters
  #============================================================================#
  results$allparameters <- list()
  
  # S
  synth <- fadata$synthesis$results
  df <- data.frame(Sample = synth$Sample, FA = synth$FA, S = synth$S)
  df <- as.data.frame(tidyr::pivot_wider(df, names_from = "Sample", 
                                         values_from = "S"))
  
  if ("elongation"%in% names(fadata)){
    elo <- fadata$elongation$results
    
    df2 <- data.frame(Sample = elo$Sample, FA = elo$FA,
                          S = elo$S16)
    df2 <- as.data.frame(tidyr::pivot_wider(df2, names_from = "Sample",
                                         values_from = "S"))
    
    df <- rbind(df, df2)
  }
  
  S16 <- df
  rownames(S16) <- df$FA
  S16 <- S16[,colnames(S16) != "FA", drop=FALSE]
  results$allparameters$S16 <- S16
  
  
  if ("elongation" %in% names(fadata)){
    
    # E1
    df <- data.frame(Sample = elo$Sample, FA = elo$FA, 
                     E1 = elo$E1)
    df <- as.data.frame(pivot_wider(df, names_from = "Sample",
                                         values_from = "E1"))
    E1 <- df
    rownames(E1) <- df$FA
    E1 <- E1[,colnames(E1) != "FA", drop=FALSE]
    results$allparameters$E1 <- E1
    
    
    # E2
    df <- data.frame(Sample = elo$Sample, FA = elo$FA, 
                     E2 = elo$E2)
    df <- as.data.frame(pivot_wider(df, names_from = "Sample",
                                    values_from = "E2"))
    E2 <- df
    rownames(E2) <- df$FA
    E2 <- E2[,colnames(E2) != "FA", drop=FALSE]
    results$allparameters$E2 <- E2
    
    
    # E3
    df <- data.frame(Sample = elo$Sample, FA = elo$FA, 
                     E3 = elo$E3)
    df <- as.data.frame(pivot_wider(df, names_from = "Sample",
                                    values_from = "E3"))
    E3 <- df
    rownames(E3) <- df$FA
    E3 <- E3[,colnames(E3) != "FA", drop=FALSE]
    results$allparameters$E3 <- E3
    
    
    # E4
    df <- data.frame(Sample = elo$Sample, FA = elo$FA, 
                     E4 = elo$E4)
    df <- as.data.frame(pivot_wider(df, names_from = "Sample",
                                    values_from = "E4"))
    E4 <- df
    rownames(E4) <- df$FA
    E4 <- E4[,colnames(E4) != "FA", drop=FALSE]
    results$allparameters$E4 <- E4
    
    
    # E5
    df <- data.frame(Sample = elo$Sample, FA = elo$FA, 
                     E5 = elo$E5)
    df <- as.data.frame(pivot_wider(df, names_from = "Sample",
                                    values_from = "E5"))
    E5 <- df
    rownames(E5) <- df$FA
    E5 <- E5[,colnames(E5) != "FA", drop=FALSE]
    results$allparameters$E5 <- E5
    
  }
  
  fadata$results <- results

  return(fadata)
}



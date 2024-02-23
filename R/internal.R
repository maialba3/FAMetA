
# msbatch2fadata
#' Extract FA data from an annotated msbatch.
#'
#' Extract FA data from an annotated msbatch.
#'
#' @param msbatch annotated msbatch.
#' @param faid data frame with two columns (ID and Compound) specifying FA ids
#' and FA names. FA names must be unique and omega series must be indicated
#' (i.e. FA(20:4)n3, FA(24:1)n9, FA(16:0)). Unknown FA series can be named as nx,
#' ny, nz to differentiate between isomers.
#'
#' @return fadata.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
msbatch2fadata <- function(msbatch, faid){

  #============================================================================#
  # Check arguments
  #============================================================================#
  ##############################################################################
  # check msbatch structure
  if (!"fas" %in% names(msbatch)){
    stop("FA must be annotated. Use annotateFA(msbatch).")
  }
  if (!is.list(msbatch) | !all(c("metaData", "msobjects", "alignment", "grouping", "features") %in% names(msbatch)) |
      !is.data.frame(msbatch$metaData) | !is.list(msbatch$msobjects) | !is.list(msbatch$alignment) |
      !is.list(msbatch$grouping) | !is.data.frame(msbatch$features)){
    stop("Wrong msbatch format")
  }
  ##############################################################################
  # check faid structure
  if (!is.data.frame(faid)){
    stop("faid must be a data.frame with two columns: ID and FAid.")
  }
  if (!all(c("ID", "FAid") %in% colnames(faid))){
    stop("faid must be a data.frame with two columns: ID and FAid.")
  }
  if (any(duplicated(faid$FAid))){
    stop("Compound names must be unique.")
  }
  if (any(duplicated(faid$ID))){
    stop("ID must be unique.")
  }

  #============================================================================#
  # Change data structure
  #============================================================================#
  fas <- merge(faid[,c("ID"), drop=FALSE], msbatch$isotopes, by="ID")

  if ("Isotope" %in% colnames(fas)){
    label <- fas$Isotope
    label <- as.numeric(unlist(sapply(label, function(x) {
      x <- gsub("\\[M\\+", "", x)
      x <- gsub("]", "", x)
      x
    })))
    fas$Label <- label
  }
  fas <- fas[order(fas$FAid, fas$Label),]


  fattyacids <- fas[,c("ID", "FAid", "Label")]
  colnames(fattyacids) <- c("ID", "Compound", "Label")
  if (length(unique(fattyacids$ID)) != length(unique(fattyacids$Compound))){
    stop("Unique Compound names must be provided for each ID")
  }
  fattyacids <- fattyacids[,c("Compound", "Label")]


  data <- fas[,colnames(fas) %in% make.names(msbatch$metaData$sample)]
  data <- data[order(fattyacids$Compound, fattyacids$Label),]
  fattyacids <- fattyacids[order(fattyacids$Compound, fattyacids$Label),]
  rownames(data) <- paste(fattyacids$Compound, "_M+", fattyacids$Label, sep="")

  if ("IS" %in% names(msbatch)){
    IS <- msbatch$IS
  } else {
    IS <- t(data.frame(IS = rep(1, ncol(data))))
    colnames(IS) <- colnames(data)
  }

  # create fadata
  fadata <- list(metadata = msbatch$metaData, fattyacids = fattyacids, IS = IS,
                 intensities = data)
  fadata$metadata <- fadata$metadata[order(colnames(fadata$intensities)),]
  if ("IS" %in% names(fadata)){
    fadata$IS <- fadata$IS[,order(colnames(fadata$intensities)), drop=FALSE]
  }
  fadata$intensities <- fadata$intensities[order(fadata$fattyacids$Compound, fadata$fattyacids$Label),
                                           order(colnames(fadata$intensities))]
  fadata$fattyacids <- fadata$fattyacids[order(fadata$fattyacids$Compound, fadata$fattyacids$Label),]

  return(fadata)
}

# combAcetate
#' acetate combinations for M+0, M+1 and M+2.
#'
#' acetate combinations for M+0, M+1 and M+2.
#'
#' @param M total number of carbons for the FA.
#'
#' @return acetate combination for M carbon atoms.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
combAcetate <- function(M){
  if (M%%2 != 0){
    stop("M must be an even number")
  }

  if (M > 0){
    N <- M/2

    perm <- gtools::permutations(n=3, r=N, v=c(0,1,2), repeats.allowed=TRUE)
    sumperm <- apply(perm, 1, sum)

    tperm <- t(apply(perm, 1, function(x){
      zero <- sum(x == 0)
      one <- sum(x ==1)
      two <- sum(x == 2)
      return(c(zero = zero, one = one, two = two))
    }))

    comb <- sapply(c(0:M), function(x){
      c <- unique(tperm[which(sumperm == x),, drop=FALSE])
      return(c)
    }, simplify = FALSE)
  } else {
    comb <- as.matrix(data.frame(zero = 0, one = 0, two = 0))
  }

  return(comb)
}

# synthesisqmult
#' calculate FA isotope distribution for newly synthesized FAs using the quasi
#' multinomial distribution.
#'
#' calculate FA isotope distribution for newly synthesized FAs using the quasi
#' multinomial distribution.
#'
#' @param D1 tracer contribution to M+1 acetate pool.
#' @param D2 tracer contribution to M+2 acetate pool.
#' @param P overdispersion parameter. If different to 0, quasi-multinomial
#' distribution is obtained, while if set to 0, it is simplified to a multinomial
#' distribution.
#' @param S fraction of newly synthesized FA.
#' @param M total number of carbons for the FA.
#' @param vcomb list of acetate combinations obtained with \link{combAcetate}
#' function.
#'
#' @return numeric vector with the FA isotope distribution.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
synthesisqmult <- function(D1, D2, P, S, M, vcomb){

  N <- M/2             # N molecules of acetate
  D0 <- 1 - D1 - D2    # contributions from M+0, M+1 and M+2 to the acetate pool
  if (D0 == 0){D0 <- 1e-20}
  if (D1 == 0){D1 <- 1e-20}
  if (D2 == 0){D2 <- 1e-20}

  # quasi-multinomial distribution with 2 terms: Synthesis and Importation (1-S)
  dist <- S * unlist(lapply(vcomb, function(x){
    nc <- factorial(N)/(factorial(x[,"zero"]) * factorial(x[,"one"]) * factorial(x[,"two"]))
    sum(nc * (1+N*P) *
          D0/(1+N*P) * ((D0 + x[,"zero"]*P)/(1+N*P))^(x[,"zero"]-1) *
          D1/(1+N*P) * ((D1 + x[,"one"]*P)/(1+N*P))^(x[,"one"]-1) *
          D2/(1+N*P) * ((D2 + x[,"two"]*P)/(1+N*P))^(x[,"two"]-1))
  })) + (1-S) * 0^c(0:M)

  # check results
  if (any(is.na(dist)) | any(dist < 0) |  round(sum(dist), 4) != 1 | D0 < 0 |
      P > ((1-max(D0,D1,D2))/N)){dist <- rep(0, M+1)}
  return(dist)
}

# elongationqmult
#' calculate FA isotope distribution for elongated FAs using the quasi
#' multinomial distribution.
#'
#' calculate FA isotope distribution for elongated FAs using the quasi
#' multinomial distribution.
#'
#' @param S16 fraction of newly synthesized palmitate.
#' @param D1 tracer contribution to M+1 acetate pool.
#' @param D2 tracer contribution to M+2 acetate pool.
#' @param P overdispersion parameter. If different to 0, quasi-multinomial
#' distribution is obtained, while if set to 0, it is simplified to a multinomial
#' distribution.
#' @param E1 fraction of elongated C18 FA from C16.
#' @param E2 fraction of elongated C20 FA from C18.
#' @param E3 fraction of elongated C22 FA from C20.
#' @param E4 fraction of elongated C24 FA from C22.
#' @param E5 fraction of elongated C26 FA from C24.
#' @param M total number of carbons for the FA.
#' @param vcomb16 list of acetate combinations for C16 synthesis obtained with
#' combAcetate(16) function.
#' @param mcombe list of acetate combinations for each elongation step obtained
#' with combAcetate(2) function.
#'
#' @return FA isotope distribution.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
elongationqmult <- function(S16, D1, D2, P,
                E1, E2, E3, E4, E5,
                M, vcomb16, mcombe){

  D0 <- 1 - D1 - D2
  if (D0 == 0){D0 <- 1e-20}
  if (D1 == 0){D1 <- 1e-20}
  if (D2 == 0){D2 <- 1e-20}
  E <- c(E1, E2, E3, E4, E5)
  n <- (M/2)-8

  # De novo synthesis and importation for FA(16:0)
  dist <- S16 * unlist(lapply(vcomb16, function(x){
    nc <- factorial(8)/(factorial(x[,"zero"]) * factorial(x[,"one"]) * factorial(x[,"two"]))
    sum(nc * (1+8*P) *
          D0/(1+8*P) * ((D0 + x[,"zero"]*P)/(1+8*P))^(x[,"zero"]-1) *
          D1/(1+8*P) * ((D1 + x[,"one"]*P)/(1+8*P))^(x[,"one"]-1) *
          D2/(1+8*P) * ((D2 + x[,"two"]*P)/(1+8*P))^(x[,"two"]-1))
  })) + (1-S16) * 0^c(0:16)

  # elongation steps from FA(16:0): each step accounts for elongation of
  # synthesized and imported precursors (Ix = 1-Ex)
  if (n > 0){
    for (e in 1:n){
      dist <- c(dist, rep(0, 2))
      dist <- E[e] * apply(t(sapply(c(1:(length(dist)-2)), function(i) {
        c(rep(0, (i-1)), dist [i] * sapply(1:length(mcombe), function (c) {
          sum(1 * D0^mcombe[[c]][,"zero"] * (D1)^mcombe[[c]][,"one"] * D2^mcombe[[c]][,"two"])
        }), rep(0,(length(dist)-i-2)))
      })), 2, sum) +
        c(1-E[e], rep(0, length(dist)-1))
    }
  }

  # check results
  if (any(E > 1) | any(E < 0) | any(dist < 0) |  round(sum(dist), 4) != 1 |
      D0 < 0){dist <- rep(0, M+1)}

  return(dist)
}

# runSynthesisAnalysis
#' run synthesis analysis.
#'
#' run synthesis analysis.
#'
#' @param fadata fadata containing synthesis results.
#' @param toDo fatty acids to analyse.
#' @param R2Thr positive numeric between 0 and 1 specifying the minimum R2
#' allowed for fits.
#' @param maxiter parameter passed to \link{nls.control}. Positive integer
#' specifying the maximum number of iterations allowed.
#' @param maxconvergence positive integer specifying the maximum number of
#' successes before choosing the winning model.
#' @param D1 positive numeric vector with values between 0 and 1 specifying the
#' contribution of acetate M+1. If NA it is estimated.
#' @param D2 positive numeric vector with values between 0 and 1 specifying the
#' contribution of acetate M+2. If NA it is estimated.
#' @param P overdispersion parameter. If NA it is estimated (quasi-multinomial
#' distribution). If set to 0, no overdispersion is assumed (multinomial
#' distribution).
#' @param startpoints positive integer specifying the number of starting points
#' for each parameter to be estimated.
#' @param parameters parameters to be estimated for each fatty acid. It can be
#' modified to change them or to add new fatty acids (adding new rows).
#' @param verbose print information messages.
#'
#' @return De novo-synthesis analysis results.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
runSynthesisAnalysis <- function(fadata, toDo, R2Thr = R2Thr,
                                 maxiter = maxiter, maxconvergence = maxconvergence,
                                 D1 = D1, D2 = D2, P = P, startpoints = startpoints,
                                 parameters = FAMetA::parameters,
                                 verbose = TRUE){

  fas <- data.frame(cbind(fadata$fattyacids[,c("Compound", "Label")],
                          fadata$mid), stringsAsFactors = FALSE)
  fas <- split(fas, fas$Compound)

  if(length(toDo) > 0){
    resultsSynth <- list()
    for (c in toDo){
      if(verbose){cat(paste("\n", c, "...", sep=""))}
      fa <- fas[[c]]
      fa <- fa[order(fa$Label),, drop = FALSE]
      M <- as.numeric(parameters[parameters$FattyAcid == c, "M"])
      resultsSynth[[c]] <- estimateSynthesis(fa=fa, M=M, R2Thr = R2Thr,
                                             maxiter = maxiter,
                                             maxconvergence = maxconvergence,
                                             D1 = D1, D2 = D2, P = P,
                                             startpoints = startpoints)
      if(verbose){cat("OK")}
    }
    #==========================================================================#
    # Summarize results
    #==========================================================================#
    results <- cbind(FA = rep(toDo, each = ncol(fadata$intensities)),
                             Sample = colnames(fadata$intensities),
                             Group = fadata$metadata$sampletype,
                             do.call(plyr::rbind.fill, lapply(resultsSynth, function(x) x[[1]])))
    predictedValues <- data.frame(do.call(rbind, lapply(resultsSynth, function(x)
      do.call(cbind, lapply(x$models, function(y) if (!any(is.na(y))){predict(y)}else{rep(NA, nrow(x$mid))})))))
    compounds <- do.call(rbind, lapply(resultsSynth, function(x) x[[2]][,c("Compound", "Label")]))
    predictedValues <- cbind(compounds, predictedValues)
    rownames(predictedValues) <- paste(predictedValues$Compound, "_M+", 
                                       predictedValues$Label, sep = "")

    plots <- plotResults(resultsSynth, fadata$metadata$sampletype, toDo)

    results <- list(fas = toDo, results = results, predictedValues = predictedValues,
                    plots = plots, details = resultsSynth)

  } else {
    stop("Unable to run synthesis analysis.")
  }

  return(results)
}

# estimateSynthesis
#' estimate synthesis parameters.
#'
#' estimate synthesis parameters.
#'
#' @param fa data frame with isotope intensities for a FA. First two columns
#' correspond to Compound and Label information.
#' @param M total number of carbons for the FA.
#' @param R2Thr positive numeric between 0 and 1 specifying the minimum R2
#' allowed for fits.
#' @param maxiter parameter passed to \link{nls.control}. Positive integer
#' specifying the maximum number of iterations allowed.
#' @param maxconvergence positive integer specifying the maximum number of
#' successes before choosing the winning model.
#' @param D1 positive numeric vector with values between 0 and 1 specifying the
#' contribution of acetate M+1. If NA it is estimated.
#' @param D2 positive numeric vector with values between 0 and 1 specifying the
#' contribution of acetate M+2. If NA it is estimated.
#' @param P overdispersion parameter. If NA it is estimated (quasi-multinomial
#' distribution). If set to 0, no overdispersion is assumed (multinomial
#' distribution).
#' @param startpoints positive integer specifying the number of starting points
#' for each parameter to be estimated.
#'
#' @return De novo-synthesis and elongation analysis results.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
estimateSynthesis <- function(fa, M, R2Thr = 0.9, maxiter = 1e3,
                                 maxconvergence = 100, D1 = NA, D2 = NA,
                                 P = NA, startpoints = 5){

  # Combinations of acetate molecules for each isotopologue
  vcomb <- combAcetate(M)

  # check distributions
  values <- apply(fa[,3:ncol(fa), drop = FALSE], 2, function(x) x/sum(x))
  values <- data.frame(apply(values, 2, function(x){
    if(all(is.na(x))){
      x <- rep(0, length(x))
    }
    return(x)
  }))

  # formula for synthesis distributions
  formula <- as.formula(resp ~ synthesisqmult(D1, D2, P, S, M, vcomb))

  # fit distributions
  results <- list()
  models <- list()

  rowsResults <- c("D0", "D1", "D2", "P", "S", "I", "R2")
  for (s in 1:ncol(values)){
    resp <- as.vector(values[,s])

    gridStart <- unique(gridSynth(D1=D1[s], D2=D2[s], P=P[s], M=M,
                                  startpoints=startpoints))
    datanls <- dataSynth(D1=D1[s], D2=D2[s], P=P[s], M=M, vcomb=vcomb)
    datanls$resp <- resp

    if (!(all(round(resp,4) == c(1, rep(0, M))) | all(resp == rep(0, M+1)))){
      model <- tryCatch({estimatePars(formula, gridStart, datanls,
                                      maxiter = maxiter,
                                      maxconvergence=maxconvergence)},
                        error = function(e){NULL})

      if (!is.null(model)){
        if (is.na(P[s]) & is.na(D2[s])){
          if (coef(model)["D2"] >= 0.75 & coef(model)["P"] >= 0.5*(1-coef(model)["D2"])/(M/2)){
            model <- tryCatch({estimatePars(formula, gridStart, datanls,
                                            maxiter = maxiter,
                                            maxconvergence=maxconvergence,
                                            limitPhi = 0.25*(1-coef(model)["D2"])/(M/2))}, 
                              error = function(e){NULL})
          }
        }
      }
    } else {
      model <- NULL
    }

    # save results
    if (!is.null(model)){
      R2 <- c(1 - (sum((resp - predict(model))^2)/sum((resp - mean(predict(model)))^2)))
      if (R2 < R2Thr){
        results[[s]] <- NA[rowsResults]
        models[[s]] <- NA
      } else {
        results[[s]] <- c(coef(model), I = 1-as.numeric(coef(model)["S"]), R2 = R2)
        if (!is.na(D1[s])){results[[s]]["D1"] <- D1[s]}
        if (!is.na(D2[s])){results[[s]]["D2"] <- D2[s]}
        if (!is.na(P[s])){results[[s]]["P"] <- P[s]}
        if (results[[s]]["I"] == 1){
          results[[s]]["S"] <- results[[s]]["D1"] <- results[[s]]["D2"] <- results[[s]]["P"] <- 0
        }
        results[[s]]["D0"] <- 1 - results[[s]]["D1"] - results[[s]]["D2"]
        results[[s]] <- results[[s]][rowsResults]
        models[[s]] <- model
      }
    } else {
      results[[s]] <- NA[rowsResults]
      models[[s]] <- NA
    }
  }
  results <- data.frame(do.call(rbind, results))
  colnames(results) <- rowsResults
  rownames(results) <-  names(models) <- colnames(fa)[3:ncol(fa)]
  values <- data.frame(cbind(fa[,1:2], values), stringsAsFactors = F)

  return(list(results = results, mid = values, models = models))
}

# runElongationAnalysis
#' run the elongation analysis.
#'
#' run the elongation analysis.
#'
#' @param fa data frame with isotope intensities for a FA. First two columns
#' correspond to Compound and Label information.
#' @param M total number of carbons for the FA.
#' @param D1 positive numeric between 0 and 1 specifying the contribution of
#' acetate M+1. Estimated with \link{synthesisAnalysis}.
#' @param D2 positive numeric between 0 and 1 specifying the contribution of
#' acetate M+2. Estimated with \link{synthesisAnalysis}.
#' @param P overdispersion parameter. Estimated with \link{synthesisAnalysis}.
#' @param S16 fraction of newly synthesized C16 FA. If NA it is estimated.
#' It is set to 0 for n3 and n6 FA series.
#' @param E1 fraction of elongated C18 FA from C16. If NA it is estimated.
#' It is set to 0 for n3 and n6 FA series.
#' @param E2 fraction of elongated C20 FA from C18. If NA it is estimated.
#' @param E3 fraction of elongated C22 FA from C20. If NA it is estimated.
#' @param E4 fraction of elongated C24 FA from C22. If NA it is estimated.
#' @param E5 fraction of elongated C26 FA from C24. If NA it is estimated.
#' @param R2Thr positive numeric between 0 and 1 specifying the minimum R2
#' allowed for fits.
#' @param maxiter parameter passed to \link{nls.control}. Positive integer
#' specifying the maximum number of iterations allowed.
#' @param maxconvergence positive integer specifying the maximum number of
#' successes before choosing the winning model.
#' @param startpoints positive integer specifying the number of starting points
#' for each parameter to be estimated.
#' @param D2Thr minimum D2 value allowed to perform the elongation analysis.
#'
#' @return Elongation and importation analysis results.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
runElongationAnalysis <- function(fa, M, D1, D2, P, S16, E1, E2, E3, E4, E5, R2Thr = 0.9,
                   maxiter = 1e4, maxconvergence = 100, startpoints = 5, D2Thr = D2Thr){

  vcomb16 <- combAcetate(16)
  mcombe <- combAcetate(2) # acetate combination for elongation steps

  SEini <- list(S16, E1, E2, E3, E4, E5)

  # distributions
  values <- apply(fa[,3:ncol(fa), drop = FALSE], 2, function(x) x/sum(x))
  values <- data.frame(apply(values, 2, function(x){
    if(all(is.na(x))){
      x <- rep(0, length(x))
    }
    return(x)
  }))

  formula <- as.formula(resp ~ elongationqmult(S16, D1, D2, P, E1, E2, E3, E4, E5,
                                   M, vcomb16, mcombe))

  # Run FASA
  results <- list()
  models <- list()

  rowsResults <- c("D0", "D1", "D2", "P", "S16", "E1", "E2", "E3", "E4", "E5", "I", "R2")
  for (s in 1:ncol(values)){
    resp <- as.vector(values[,s])
    
    # check parameters based on distribution
    n <- (M-16)/2
    if (!is.na(S16)) {imported <- TRUE} else {imported <- FALSE}
    SE <- checkparameters(resp, M, n, SEini, imported, D2=D2[s])
    
    gridStart <- gridElong(S16=SE[[1]], E1 = SE[[2]], E2=SE[[3]], E3=SE[[4]],
                           E4=SE[[5]], E5=SE[[6]], startpoints=startpoints)

    datanls <- dataElong(S16=SE[[1]], D1=D1[s], D2=D2[s], P=P[s], E1 = SE[[2]], 
                         E2=SE[[3]], E3=SE[[4]], E4=SE[[5]], E5=SE[[6]],
                       M=M, vcomb16=vcomb16, mcombe=mcombe)
    datanls$resp <- resp
    
    if (!any(is.na(c(D2[s], D1[s], P[s])))){
      if (D2[s]+D1[s]*0.5 >= D2Thr){
        if (!all(resp == rep(0, M+1)) | any(unlist(lapply(SE[2:6], is.na)))){
          model <- tryCatch({estimatePars(formula, gridStart, datanls,
                                          maxiter = maxiter, maxconvergence=maxconvergence)},
                            error = function(e){NULL})
        } else {
          model <- NULL
        }
      } else {
        model <- NULL
      }
    } else {
      model <- NULL
    }

    if (!is.null(model) & !any(is.na(model))){
      R2 <- c(1 - (sum((resp - predict(model))^2)/sum((resp - mean(predict(model)))^2)))
      results[[s]] <- c(coef(model), P = as.numeric(P[s]), R2 = R2)
      if (R2 < R2Thr){
        results[[s]] <- NA[rowsResults]
        models[[s]] <- NA
      } else {
        if (!is.na(D1[s])){results[[s]]["D1"] <- D1[s]}
        if (!is.na(D2[s])){results[[s]]["D2"] <- D2[s]}
        if (!is.na(P[s])){results[[s]]["P"] <- P[s]}
        results[[s]]["D0"] <- 1 - results[[s]]["D1"] - results[[s]]["D2"]
        if (!is.na(SE[[1]])){results[[s]]["S16"] <- SE[[1]]}
        results[[s]][c("E1", "E2", "E3", "E4", "E5")[!unlist(lapply(SE[2:6], is.na))]] <- 0
        results[[s]]["I"] <- 1 - results[[s]][c("E1", "E2", "E3", "E4", "E5")[(M/2)-8]]
        results[[s]] <- results[[s]][rowsResults]
        models[[s]] <- model
      }
    } else {
      results[[s]] <- NA[rowsResults]
      models[[s]] <- NA
    }
  }
  results <- data.frame(do.call(rbind, results))
  colnames(results) <- rowsResults
  rownames(results) <-  names(models) <- colnames(fa)[3:ncol(fa)]
  values <- data.frame(cbind(fa[,1:2], values),
                       stringsAsFactors = F)
  return(list(results = results, mid = values, models = models))
}

# checkparameters
#' check S and E steps based on isotopologue distribution.
#'
#' check S and E steps based on isotopologue distribution.
#'
#' @param resp isotopologue distribution
#' @param M total number of carbons of the fatty acid.
#' @param n maximum number of elongation steps.
#' @param SE list with S and E parameters. NA indicates they must be estimated 
#' while 0 indicates it does not occur.
#' @param imported logical. TRUE if S16 is predefined as 0 (n3 or n6 series). 
#' @param D2 numeric between 0 and 1. Only if D2 >= 0.4 parameters are checked 
#' based on distribution.
#'
#' @return SE list corrected.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
checkparameters <- function(resp, M, n, SE, imported, D2){
  
  if (!is.na(D2)){
    if (D2 >= 0.4){
      if (!imported & (sum(resp[seq(2*n+3, M+1, 2)] > 0.0001) >= 3 | 
                       sum(resp[seq(2*n+3, M+1, 2)] > 0.001) >= 2)){  
        #floor(0.33*length(resp[seq(2*n+3, M+1, 2)]))
        SE[[1]] <- NA
      } else {
        SE[[1]] <- 0
      }
      if (!is.na(SE[[1]])){
        if (resp[2*n+1] > 0.001){
          SE[[2]] <- NA
        } else {
          SE[[2]] <- 0
        }
      }
      if (!is.na(SE[[2]]) & n > 1){
        if (resp[2*(n-1)+1] > 0.001){
          SE[[3]] <- NA
        } else {
          SE[[3]] <- 0
        }
      }
      if (!is.na(SE[[3]]) & n > 2){
        if (resp[2*(n-2)+1] > 0.001){
          SE[[4]] <- NA
        } else {
          SE[[4]] <- 0
        }
      }
      if (!is.na(SE[[4]]) & n > 3){
        if (resp[2*(n-3)+1] > 0.001){
          SE[[5]] <- NA
        } else {
          SE[[5]] <- 0
        }
      }
      if (!is.na(SE[[5]]) & n > 4){
        if (resp[2*(n-4)+1] > 0.001){
          SE[[6]] <- NA
        } else {
          SE[[6]] <- 0
        }
      }
    }
  }
  
  return(SE)
}

# gridSynth
#' generate starting points grid for  synthesis analysis.
#'
#' generate starting points grid for synthesis analysis.
#'
#' @param D1 positive numeric between 0 and 1 specifying the contribution of
#' acetate M+1. If NA it is included in the starting points grid.
#' @param D2 positive numeric between 0 and 1 specifying the contribution of
#' acetate M+2. If NA it is included in the starting points grid.
#' @param P overdispersion parameter. If NA, it is included in the starting
#' points grid.
#' @param S fraction of newly synthesized C16 FA. If NA it is included in the
#' starting points grid.
#' @param M total number of carbons for the FA.
#' @param startpoints positive integer specifying the number of starting points
#' for each parameter to be estimated.
#'
#' @return Starting points grid for synthesis analysis.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
gridSynth <- function(D1=NA, D2=NA, P=NA, S=NA, M, startpoints){
  toestimate <- c()

  if (missing(D1) | is.na(D1)){
    D1 <- (0:(startpoints-1))/(startpoints-1)
    toestimate <- c(toestimate, "D1")
  } else {
    D1 <- D1
  }
  if (missing(D2) | is.na(D2)){
    D2 <- (0:(startpoints-1))/(startpoints-1)
    toestimate <- c(toestimate, "D2")
  } else {
    D2 <- D2
  }

  gridContrib <- expand.grid(D1 = D1, D2 = D2)
  gridContrib <- gridContrib[apply(gridContrib, 1, sum) <= 1,]

  if (missing(P) | is.na(P)){
    iniPhi <- 0 # -1/(M/2)
    endPhi <- 0.1 # 1/(M/2)

    gridOD <- seq(iniPhi, endPhi, length=startpoints)
    gridMult <- as.data.frame(tidyr::crossing(gridContrib, gridOD))
    colnames(gridMult) <- c("D1", "D2", "P")
    gridMult <- gridMult[apply(gridMult, 1, function(x) x["P"] >= -min(x[c("D1", "D2")])/(M/2) &
                                 x["P"] <= (1-max(x[c("D1", "D2")]))/(M/2)),]
    toestimate <- c(toestimate, "P")
  } else {
    gridMult <- gridContrib
    gridMult$P <- P
  }

  gridSE <- expand.grid(S = 0:(startpoints-1))/(startpoints-1)
  toestimate <- c(toestimate, "S")

  gridFinal <- as.data.frame(tidyr::crossing(gridMult, gridSE))
  colnames(gridFinal) <- c("D1", "D2", "P", "S")
  gridFinal <- gridFinal[,toestimate, drop = FALSE]

  return(gridFinal)
}

# gridElong
#' generate starting points grid for elongation analysis.
#'
#' generate starting points grid for elongation analysis.
#'
#' @param S16 fraction of newly synthesized C16 FA. If NA it is included in
#' the starting points grid.
#' @param E1 fraction of elongated C18 FA from C16. If NA it is included in
#' the starting points grid.
#' @param E2 fraction of elongated C20 FA from C18. If NA it is included in
#' the starting points grid.
#' @param E3 fraction of elongated C22 FA from C20. If NA it is included in
#' the starting points grid.
#' @param E4 fraction of elongated C24 FA from C22. If NA it is included in
#' the starting points grid.
#' @param E5 fraction of elongated C26 FA from C124. If NA it is included in
#' the starting points grid.
#' @param startpoints positive integer specifying the number of starting points
#' for each parameter to be estimated.
#'
#' @return Starting points grid for elongation analysis.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
gridElong <- function(S16=NA, E1=NA, E2=NA, E3=NA, E4=NA, E5=NA, startpoints){
  toestimate <- c()

  if (is.na(S16)){
    S16 <- (0:(startpoints-1))/(startpoints-1)
    toestimate <- c(toestimate, "S16")
  }else{
    S16 <- S16
  }
  if (is.na(E1)){
    E1 <- (0:(startpoints-1))/(startpoints-1)
    toestimate <- c(toestimate, "E1")
  }else{
    E1 <- E1
  }
  if (is.na(E2)){
    E2 <- (0:(startpoints-1))/(startpoints-1)
    toestimate <- c(toestimate, "E2")
  }else{
    E2 <- E2
  }
  if (is.na(E3)){
    E3 <- (0:(startpoints-1))/(startpoints-1)
    toestimate <- c(toestimate, "E3")
  }else{
    E3 <- E3
  }
  if (is.na(E4)){
    E4 <- (0:(startpoints-1))/(startpoints-1)
    toestimate <- c(toestimate, "E4")
  }else{
    E4 <- E4
  }
  if (is.na(E5)){
    E5 <- (0:(startpoints-1))/(startpoints-1)
    toestimate <- c(toestimate, "E5")
  }else{
    E5 <- E5
  }

  gridSE <- expand.grid(S16 = S16,
                        E1 = E1,
                        E2 = E2,
                        E3 = E3,
                        E4 = E4,
                        E5 = E5)

  colnames(gridSE) <- c("S16", "E1", "E2", "E3", "E4", "E5")
  gridFinal <- gridSE[,toestimate, drop = FALSE]

  return(gridFinal)
}

# dataSynth
#' data list for synthesis analysis.
#'
#' data list for synthesis analysis.
#'
#' @param D1 positive numeric between 0 and 1 specifying the contribution of
#' acetate M+1. If different to NULL it is included in the data list.
#' @param D2 positive numeric between 0 and 1 specifying the contribution of
#' acetate M+2. If different to NULL it is included in the data list.
#' @param P overdispersion parameter. If different to NULL it is included in the
#' data list.
#' @param M total number of carbons for the FA.
#' @param vcomb list of acetate combinations for C16 synthesis.
#'
#' @return Data list for synthesis analysis.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
dataSynth <- function(D1, D2, P, M, vcomb){
  datanls <- list(D1 = D1,
                  D2 = D2,
                  P = P,
                  M = M,
                  vcomb = vcomb)
  datanls <- datanls[!unlist(lapply(datanls, function(x) any(is.na(x))))]
  return(datanls)
}

# dataElong
#' data list for elongation analysis.
#'
#' data list for elongation analysis.
#'
#' @param S16 fraction of newly synthesized C16 FA. If different to NULL it is
#' included in the data list.
#' @param D1 positive numeric between 0 and 1 specifying the contribution of
#' acetate M+1.
#' @param D2 positive numeric between 0 and 1 specifying the contribution of
#' acetate M+2.
#' @param P overdispersion parameter.
#' @param E1 fraction of elongated C18 FA from C16. If different to NULL it is
#' included in the data list.
#' @param E2 fraction of elongated C20 FA from C18. If different to NULL it is
#' included in the data list.
#' @param E3 fraction of elongated C22 FA from C20. If different to NULL it is
#' included in the data list.
#' @param E4 fraction of elongated C24 FA from C22. If different to NULL it is
#' included in the data list.
#' @param E5 fraction of elongated C26 FA from C124. If different to NULL it is
#' included in the data list.
#' @param M total number of carbons for the FA.
#' @param vcomb16 list of acetate combinations for C16 synthesis.
#' @param mcombe list of acetate combinations for each elongation step.
#'
#' @return data list for elongation analysis.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
dataElong <- function(S16, D1, D2, P, E1, E2, E3, E4, E5,
                      M, vcomb16, mcombe){
  datanls <- list(S16 = S16,
                  D1 = D1,
                  D2 = D2,
                  P = P,
                  E1 = E1,
                  E2 = E2,
                  E3 = E3,
                  E4 = E4,
                  E5 = E5,
                  M = M,
                  vcomb16 = vcomb16,
                  mcombe = mcombe)
  datanls <- datanls[!unlist(lapply(datanls, function(x) any(is.na(x))))]
  return(datanls)
}

# estimatePars
#' Fit FA isotope distribution with non-linear regression using multiple
#' starting points to find the best fit.
#'
#' Fit FA isotope distribution with non-linear regression using multiple
#' starting points to find the best fit.
#'
#' @param formula formula.
#' @param gridStart starting points grid.
#' @param datanls data list.
#' @param maxiter parameter passed to \link{nls.control}. Positive integer
#' specifying the maximum number of iterations allowed.
#' @param maxconvergence positive integer specifying the maximum number of
#' successes before choosing the winning model.
#' @param limitPhi upper bound for overdispersion parameter.
#'
#' @references Daniel Padfield and Granville Matheson (2020). nls.multstart: Robust Non-Linear Regression using AIC Scores. R package version 1.2.0. <https://CRAN.R-project.org/package=nls.multstart>
#'
#' @return model.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
estimatePars <- function(formula, gridStart, datanls, maxiter, maxconvergence,
                         limitPhi = 0.1){
  control <- minpack.lm::nls.lm.control(maxiter = maxiter)

  allfits <- c()
  count <- 0
  for (t in 1:nrow(gridStart)){
    if (count == maxconvergence){break}
    f <- Inf
    f <- fitSSE(formula, startpars=gridStart[t,, drop=FALSE], datanls, control,
                limitPhi = limitPhi)
    if (!is.infinite(f) & f < 0.1){
      count <- count + 1
    }
    allfits <- c(allfits, f)
  }
  orderallfits <- order(allfits, decreasing = FALSE)
  head(allfits[orderallfits])

  fit_best <- NULL
  lower <- rep(0, ncol(gridStart))
  upper <- rep(1, ncol(gridStart))
  if("P" %in% colnames(gridStart)){
    lower[colnames(gridStart) == "P"] <- 0
    upper[colnames(gridStart) == "P"] <- limitPhi
  }

  try(
    fit_best <- minpack.lm::nlsLM(
      formula,
      start = gridStart[orderallfits[1],, drop=FALSE],
      control = control,
      data = datanls,
      lower = lower,
      upper = upper
    ),
    silent <- TRUE
  )
  return(fit_best)
}

# fitSSE
#' Fit FA isotope distribution with non-linear regression and calculate the sum
#' of squared estimate of errors.
#'
#' Fit FA isotope distribution with non-linear regression and calculate the sum
#' of squared estimate of errors.
#'
#' @param formula formula.
#' @param startpars starting points for parameters to estimate.
#' @param datanls data list.
#' @param control nls control.
#'
#' @references Daniel Padfield and Granville Matheson (2020). nls.multstart: Robust Non-Linear Regression using AIC Scores. R package version 1.2.0. <https://CRAN.R-project.org/package=nls.multstart>
#'
#' @return sum of squared estimate errors.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
fitSSE <- function(formula, startpars, datanls, control, limitPhi = 0.1) {
  fit <- NULL
  lower <- rep(0, ncol(startpars))
  upper <- rep(1, ncol(startpars))
  if("P" %in% colnames(startpars)){
    lower[colnames(startpars) == "P"] <- 0
    upper[colnames(startpars) == "P"] <- limitPhi
  }

  try(
    fit <- minpack.lm::nlsLM(
      formula=formula,
      start = as.list(startpars),
      control = control,
      data = datanls,
      lower = lower,
      upper = upper
    ),
    silent <- TRUE
  )

  sse <- ifelse(!is.null(fit), sum((predict(fit) - datanls$resp)^2), Inf)

  return(sse)
}

# plotDistributions
#' Plot observed and predicted isotope distributions.
#'
#' Plot observed and predicted isotope distributions.
#'
#' @param results results.
#' @param groups character vector specifying sample groups.
#' @param title title.
#'
#' @return plot.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
plotDistributions <- function(results, groups, title = ""){
  
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar, new = FALSE))
  
  plots <- list()
  for (i in 1:length(results$models)){
    grDevices::pdf(NULL) # use a pdf NULL device to save plots to an object
    grDevices::dev.control(displaylist = "enable")
    graphics::par(mar = c(3,4,4,1), mgp=c(2,1,0), bg = "white")
    if(!is.na(results$models[i])){
      resp <- results$mid[,i+2]
      isotope <- paste("M+", seq(0, length(resp)-1), sep="")
      if(!is.null(results$models[[i]])){
        prediction <- predict(results$models[[i]])
      } else {
        prediction <- rep(0, length(resp))
      }
    } else {
      resp <- results$mid[,i+2]
      isotope <- paste("M+", seq(0, length(resp)-1), sep="")
      prediction <- rep(0, length(resp))
    }

    coul <- c("#66C2A5", "#FC8D62")
    barplot(t(cbind(resp, prediction)), col = coul,
            beside = T, names = rep("", length(isotope)),
            main = paste(paste(title,
                               paste(paste(colnames(results$results),
                                           round(results$results[i,],4),
                                           sep=":"), collapse="; "),
                               sep = "\n"),
                         paste(colnames(results$mid)[i+2], groups[i], sep="; "), sep="\n"),
            cex.main = 0.6,
            cex.axis = 0.6,
            cex.lab = 0.6,
            ylab = "Relative Intensity",
            ylim = c(0,1),
            las = 1)
    lab <- isotope
    at <- (1:(length(isotope)*3))[seq(2, length(isotope)*3, 3)]
    axis(1, at = at, labels = lab, cex.axis = 0.6, las = 2)
    legend("topright",
           legend = c("Observed", "Predicted"),
           fill = coul[1:2], cex = 0.6)

    plots[[i]] <- grDevices::recordPlot() # save plot
    invisible(grDevices::dev.off()) # close pdf NULL device
  }

  return(plots)
}

# plotResults
#' Plot observed and predicted isotope distributions from the results obtained.
#'
#' Plot observed and predicted isotope distributions from the results obtained.
#'
#' @param results results.
#' @param groups character vector specifying sample groups.
#' @param fas character vector specifying the FAs to be plotted.
#'
#' @return plot.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
plotResults <- function(results, groups, fas){
  plots <- list()
  for(f in fas){
    plots[[f]] <- plotDistributions(results[[f]], groups, title = f)
  }
  return(plots)
}

# getResultsTable
#' Obtain the final results table from synthesis, elongation and desaturation
#' analysis.
#'
#' Obtain the final results table from synthesis, elongation and desaturation
#' analysis.
#'
#' @param fadata fadata list.
#' @param parameters parameters to be estimated for each fatty acid. It can be
#' modified to change them or to add new fatty acids.
#'
#' @return results data frame.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
getResultsTable <- function(fadata, parameters = FAMetA::parameters){

  cnames <- c("FA", "Sample", "Group", "D0", "D1", "D2", "P", "S", "E", "Des", "I")
  resultstable <- fadata$synthesis$results
  resultstable$E <- resultstable$Des <- NA
  resultstable <- resultstable[,cnames]

  if ("elongation" %in% names(fadata)){
    for (fa in fadata$elongation$fas){
      param <- parameters[parameters$FattyAcid == fa,]
      params <- c("FA", "Sample", "Group", "D0", "D1", "D2", "P",
                    tail(names(param)[param == 1], n=1), "I")
      newresult <- fadata$elongation$results[fadata$elongation$results$FA == fa,
                                             params, drop=FALSE]
      if (fa %in% c("FA(18:2)n6", "FA(18:3)n3", "FA(18:3)n6")){ # these FA are never considered endogenously synthesized
        if ("E1" %in% params){
          newresult["E1"] <- NA
          newresult["I"] <- 1
        }
      }
      if (ncol(newresult) == 9){
        colnames(newresult) <- c("FA", "Sample", "Group", "D0", "D1", "D2", "P",
                                 "E", "I")
      }
      resultstable <- plyr::rbind.fill(resultstable, newresult)
    }
  }

  if ("desaturation" %in% names(fadata)){
    for (fa in fadata$desaturation$fas){
      newresult <- fadata$desaturation$results[fadata$desaturation$results$FA == fa,
                                                c("Des", "I")]
      resultstable[resultstable$FA == fa, "Des"] <- newresult$Des
      resultstable[resultstable$FA == fa, "I"] <- newresult$I
      resultstable[resultstable$FA == fa, c("S", "E")] <- NA
    }
  }
  return(resultstable)
}


# getSummaryTable
#' Obtain a summary results table.
#'
#' Obtain a summary results table containing means and standard deviation for
#' each relevant parameter for each fatty acids by sample groups.
#'
#' @param fadata fadata list.
#' @param resultstable results data frame obtained with \link{getResultsTable}.
#' @param parameters parameters to be estimated for each fatty acid. It can be
#' modified to change them or to add new fatty acids.
#'
#' @return summary data frame.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
getSummaryTable <- function(fadata, resultstable,
                            parameters = FAMetA::parameters){

  cnames <- c("D0", "D1", "D2", "P", "S", "E", "Des", "I")
  means <- apply(resultstable[cnames], 2, function(x)
    tapply(x, paste(resultstable$FA, resultstable$Group), mean, na.rm=TRUE))
  sds <- apply(resultstable[cnames], 2, function(x)
    tapply(x, paste(resultstable$FA, resultstable$Group), sd, na.rm=TRUE))

  colnames(means) <- paste("Mean_", colnames(means), sep="")
  colnames(sds) <- paste("SD_", colnames(sds), sep="")

  summary <- c()
  cnames <- c()
  for (c in 1:ncol(means)){
    summary <- cbind(summary, means[,c], sds[,c])
    cnames <- c(cnames, colnames(means)[c], colnames(sds)[c])
  }
  colnames(summary) <- cnames
  summary[is.nan(summary)] <- NA

  fas <- sapply(rownames(summary), function(x) unlist(strsplit(x, " "))[1])
  groups <- unique(sapply(rownames(summary), function(x) unlist(strsplit(x, " "))[2]))

  summaryArranged <- c()
  cnames <- c()
  for (fa in unique(fas)){
    if (fa %in% fadata$desaturation$fas){
      columns <- c("Mean_Des", "SD_Des", "Mean_I", "SD_I")
    } else if (parameters$M[parameters$FattyAcid == fa] <= 16){
      columns <- c("Mean_D0", "SD_D0", "Mean_D1", "SD_D1", "Mean_D2", "SD_D2",
                   "Mean_P", "SD_P", "Mean_S", "SD_S", "Mean_I", "SD_I")
    } else if (parameters$M[parameters$FattyAcid == fa] > 16){
      columns <- c("Mean_E", "SD_E", "Mean_I", "SD_I")
    }
    summaryArranged <- data.frame(cbind(summaryArranged, summary[fas == fa, columns, drop=FALSE]))
    cnames <- c(cnames, paste(fa, columns, sep="_"))
  }
  colnames(summaryArranged) <- cnames
  rownames(summaryArranged) <- groups

  return(summaryArranged)

}

# heatmapResults
#' Plot results on a heatmap.
#'
#' Plot results on a heatmap.
#'
#' @param toPlot data frame with data to be plotted.
#' @param cols colors for the side bar (by group).
#' @param scale "none", "row" or "column".
#' @param breaks numeric vector with breaks for colouring. Optional.
#' @param nacolor color for NA values. Grey by default.
#'
#' @return heatmap.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
heatmapResults <- function(toPlot, cols, scale = "none", breaks,
                           nacolor = "grey"){
  
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar, new = FALSE))
  
  grDevices::pdf(NULL) # use a pdf NULL device to save plots to an object
  grDevices::dev.control(displaylist = "enable")
  graphics::par(mar = c(3,4,4,1), mgp=c(2,1,0), bg = "white")
  cexRow <- 1 - nrow(toPlot)*0.005
  if (nrow(toPlot) > 1 & ncol(toPlot) > 1){
    if (missing(breaks)){
      gplots::heatmap.2(as.matrix(toPlot),
                        col=colorRampPalette(c("blue1", "white", "red"))(19),
                        ColSideColors = cols,
                        na.color = nacolor,
                        density.info = "none",
                        Rowv=F, Colv = F, trace = "none", cexCol = 0.6, cexRow = cexRow,
                        scale = scale,
                        dendrogram = "none", key=T, lhei=c(1, 4, 25), sepcolor="black",
                        margins = c(15,10),
                        colsep=0:ncol(toPlot),
                        rowsep=0:nrow(toPlot),
                        sepwidth=c(0.01,0.01))
    } else {
      gplots::heatmap.2(as.matrix(toPlot),
                        col=colorRampPalette(c("blue1", "white", "red"))(19),
                        ColSideColors = cols,
                        na.color = nacolor,
                        breaks = breaks,
                        symkey = F,
                        density.info = "none",
                        Rowv=F, Colv = F, trace = "none", cexCol = 0.6, cexRow = cexRow,
                        scale = scale,
                        dendrogram = "none", key=T, lhei=c(1, 4, 25), sepcolor="black",
                        margins = c(15,10),
                        colsep=0:ncol(toPlot),
                        rowsep=0:nrow(toPlot),
                        sepwidth=c(0.01,0.01))
    }
  } else {
    plot(0,type='n',axes=FALSE,ann=FALSE)
  }
  hm <- grDevices::recordPlot() # save plot
  invisible(grDevices::dev.off()) # close pdf NULL device

  return(hm)
}

# mzMatch
#' mz match withing a vector of mz values
#'
#' This function searches marches between a given mz and a vector of mz values
#' with certain mass  tolerance and returns the index of the matched values. It
#' is used by identification functions to find candidates of each class of lipid
#' based on full MS information.
#'
#' @param mz mz value to be matched
#' @param mzvector vector of mz values
#' @param ppm mass error tolerance
#'
#' @return Numeric vector indicating the index of matched mz values and ppms for
#' each one of those matches (match1, ppm1, match2, ppm2, etc.)
#'
#' @keywords internal
#'
#' @references M Isabel Alcoriza-Balaguer (2021). LipidMS: Lipid Annotation for
#' LC-MS/MS DDA or DIA Data. R package version 3.0.1.
#' <https://CRAN.R-project.org/package=LipidMS>
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
mzMatch <- function(mz, mzvector, ppm){
  matches_ppm <- vector()
  for (i in 1:length(mzvector)){
    ppm_observed <- abs(((mz-mzvector[i])/mz)*1000000)
    if (ppm_observed <= ppm){
      matches_ppm <- append(matches_ppm, c(i, ppm_observed))
    }
  }
  return(matches_ppm)
}

# getFormula
#' Get formula and neutral mass for annotated compounds
#'
#' Get formula and neutral mass for annotated compounds.
#'
#' @param df data frame with the input results
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#'
#' @return Data frame
#'
#' @keywords internal
#'
#' @references M Isabel Alcoriza-Balaguer (2021). LipidMS: Lipid Annotation for
#' LC-MS/MS DDA or DIA Data. R package version 3.0.1.
#' <https://CRAN.R-project.org/package=LipidMS>
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
getFormula <- function(df, dbs){
  lipidClass <- df["Class"]
  cdb <- df["CDB"]
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (lipidClass == "BA"){
    db <- dbs$badb
  } else if (lipidClass == "Carnitine"){
    db <- dbs$carnitinesdb
  } else if (lipidClass == "Cer"){
    db <- dbs$cerdb
  } else if (lipidClass == "CL"){
    db <- dbs$cldb
  } else if (lipidClass == "DG"){
    db <- dbs$dgdb
  } else if (lipidClass == "CE"){
    db <- dbs$CEdb
  } else if (lipidClass == "FA"){
    db <- dbs$fadb
  } else if (lipidClass == "FAHFA"){
    db <- dbs$fahfadb
  } else if (lipidClass == "HFA"){
    db <- dbs$hfadb
  } else if (lipidClass == "LPC"){
    db <- dbs$lysopcdb
  } else if (lipidClass == "LPE"){
    db <- dbs$lysopedb
  } else if (lipidClass == "LPG"){
    db <- dbs$lysopgdb
  } else if (lipidClass == "LPI"){
    db <- dbs$lysopidb
  } else if (lipidClass == "LPS"){
    db <- dbs$lysopsdb
  } else if (lipidClass == "MG"){
    db <- dbs$mgdb
  } else if (lipidClass == "PC"){
    db <- dbs$pcdb
  } else if (lipidClass == "PE"){
    db <- dbs$pedb
  } else if (lipidClass == "PG"){
    db <- dbs$pgdb
  } else if (lipidClass == "PI"){
    db <- dbs$pidb
  } else if (lipidClass == "PS"){
    db <- dbs$psdb
  } else if (lipidClass == "SM"){
    db <- dbs$smdb
  } else if (lipidClass == "Sph"){
    db <- dbs$sphdb
  } else if (lipidClass == "SphP"){
    db <- dbs$sphPdb
  } else if (lipidClass == "TG"){
    db <- dbs$tgdb
  } else if (lipidClass == "TG"){
    db <- dbs$tgdb
  } else if (lipidClass == "PCo"){
    db <- dbs$pcodb
  } else if (lipidClass == "PCp"){
    db <- dbs$pcpdb
  } else if (lipidClass == "PEo"){
    db <- dbs$peodb
  } else if (lipidClass == "PEp"){
    db <- dbs$pepdb
  } else if (lipidClass == "LPCo"){
    db <- dbs$lysopcodb
  } else if (lipidClass == "LPCp"){
    db <- dbs$lysopcpdb
  } else if (lipidClass == "LPEo"){
    db <- dbs$lysopeodb
  } else if (lipidClass == "LPEp"){
    db <- dbs$lysopepdb
  } else if (lipidClass == "CerP"){
    db <- dbs$cerPdb
  } else if (lipidClass == "AcylCer"){
    db <- dbs$acylcerdb
  }
  index <- which(db$total == as.character(cdb))
  if(length(index) > 0){
    formula <- db$formula[index]
    Mn <- db$Mass[index]
  } else {
    tempdb <- LipidMS::createLipidDB(lipidClass,
                            chains = c(as.character(cdb),
                                       "0:0", "0:0", "0:0"),
                            chains2 = "0:0")[[1]]
    formula <- tempdb[tempdb$total == as.character(cdb), 1]
    Mn <- tempdb[tempdb$total == as.character(cdb), 3]
  }
  return(data.frame(Formula = as.character(formula),
                    Mn = as.numeric(Mn), stringsAsFactors = FALSE))
}

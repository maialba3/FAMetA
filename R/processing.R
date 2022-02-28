# annotateFA
#' FA annotation
#'
#' FA annotation
#'
#' @param msbatch msbatch obtained from LipidMS package.
#' @param dmz mz tolerance in ppm.
#' @param rt Optional. Numeric vector of length two specifying the rt range to
#' search for FA.
#' @param adducts character vector specifying adducts.
#' @param db FA database. Data frame with three columns: formula, total (number
#' of carbons and double bounds, i.e. "18:1") and Mass.
#'
#' @return annotated msbatch.
#'
#' @examples
#' \dontrun{
#' msbatch <- annotateFA(msbatch, dmz = 5)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
annotateFA <- function(msbatch,
                       dmz = 5,
                       rt,
                       adducts = c("M-H"),
                       db){

  #============================================================================#
  # Check arguments
  #============================================================================#
  ##############################################################################
  # check msbatch structure
  if (!msbatch$grouping$grouped){
    stop("msbatch must be grouped.")
  }
  if (!is.list(msbatch) | !all(c("metaData", "msobjects", "alignment", "grouping", "features") %in% names(msbatch)) |
      !is.data.frame(msbatch$metaData) | !is.list(msbatch$msobjects) | !is.list(msbatch$alignment) |
      !is.list(msbatch$grouping) | !is.data.frame(msbatch$features)){
    stop("Wrong msbatch format")
  }
  ##############################################################################
  # check that all msobjects have an mslevel 1
  whichmslevel1 <- which(unlist(lapply(msbatch$msobjects, function(x)
    1 %in% unique(x$metaData$scansMetadata$msLevel))))
  if (length(whichmslevel1) != nrow(msbatch$metaData)){
    warning("Removing samples with no MS1 level for alignment")
    msbatch$metaData <- msbatch$metaData[whichmslevel1,]
    msbatch$msobjects <- msbatch$msobjects[whichmslevel1]
  }

  if (missing(rt)){
    rt <- c(min(msbatch$features$RT),
            max(msbatch$features$RT))
  }
  dbs <- LipidMS::assignDB()
  if (missing(db)){
    dbs$fadb <- FAMetA::fattyacidsdb
  } else {
    dbs$fadb <- db
  }

  #============================================================================#
  # Annotate FA
  #============================================================================#
  fas <- LipidMS::findCandidates(msbatch$features[msbatch$features$isotope == "[M+0]",],
                                 db = dbs$fadb,
                                 ppm = dmz,
                                 rt = rt, adducts = adducts, rttol = 5,
                                 dbs = dbs)
  fas <- fas[order(fas$cb, fas$RT),]

  # add an identifier to trace FAs
  fas$ID <- 1:nrow(fas)

  #============================================================================#
  # Save results
  #============================================================================#
  fas <- fas[,unlist(sapply(c("ID", "cb", "adducts", "ppms",
                              colnames(msbatch$features)), match, colnames(fas)))]
  fas$cb <- paste("FA(", fas$cb, ")", sep="")
  colnames(fas)[1:4] <- c("ID", "FAid", "Adducts", "ppm")
  msbatch$fas <- fas[,c("ID", "FAid", "Adducts", "mz", "RT", "iniRT", "endRT")]

  return(msbatch)
}

# plot FAs
#' Plot FA EICs
#'
#' Plot FA EICs
#'
#' @param msbatch annotated msbatch.
#' @param dmz mz tolerance in ppm for EIC extraction.
#' @param verbose print information messages.
#'
#' @return annotated msbatch with saved plots.
#'
#' @examples
#' \dontrun{
#' msbatch <- annotateFA(msbatch, dmz = 5)
#'
#' plots <- plotFA(msbatch, dmz = 10)
#'
#' pdf("FAs.pdf")
#' for (p in 1:length(plots)){
#'   print(plots[[p]])
#' }
#' dev.off()
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
plotFA <- function(msbatch, dmz, verbose = TRUE){

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar, new = FALSE))
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

  if (missing(dmz)){
    dmz <- msbatch$grouping$parameters$dmz
  }

  plots <- list()
  counter <- 1
  ##################################################
  # for each fa (unique molecular formula)
  for (fa in unique(msbatch$fas$FAid)){
    # subset all isomers for that FA
    isomers <- msbatch$fas[msbatch$fas$FAid == fa,]
    mz <- mean(isomers$mz)
    # define RT range
    iniRT <- min(isomers$iniRT)-20
    if (iniRT < 0){iniRT <- 0}
    endRT <- max(isomers$endRT)+20

    #============================================================================#
    # Extract EIC
    #============================================================================#
    eics <- list()
    for (s in 1:length(msbatch$msobjects)){
      e <- msbatch$msobjects[[s]]$rawData$MS1[abs(msbatch$msobjects[[s]]$rawData$MS1$mz - mz)*1e6/mz <= dmz,, drop = FALSE]
      e <- e[e$RT >= iniRT & e$RT <= endRT,]
      if (nrow(e) > 0){
        eics[[s]] <- e[order(e$RT), c("RT", "int")]
      } else {
        eics[[s]] <- data.frame()
      }
    }
    maxint <- max(unlist(lapply(eics, function(x) if(nrow(x) > 0){max(x$int, na.rm = TRUE)} else {0})), na.rm = TRUE)

    #============================================================================#
    # Extract peak for each isomer
    #============================================================================#
    peaks <- list()
    for (p in 1:nrow(isomers)){
      peaks[[p]] <- list()
    }
    for (p in 1:nrow(isomers)){
      iniRT2 <- isomers$iniRT[p]
      endRT2 <- isomers$endRT[p]
      for (s in 1:length(msbatch$msobjects)){
        e <- eics[[s]]
        e <- e[e$RT >= iniRT2 & e$RT <= endRT2,]
        if (nrow(e) > 0){
          peaks[[p]][[s]] <- e[order(e$RT), c("RT", "int")]
        } else {
          peaks[[p]][[s]] <- data.frame()
        }
      }
    }

    #============================================================================#
    # Plot each isomer
    #============================================================================#
    palette <- c("#42858C", "#FE9300", "#870E75", "#3E71A8",
                 "#7F8E39", "#5F3659", "#E5C616",
                 "#16A08CFF", "#FE6900", "#628395", "#C5D86D", "#969696FF",
                 "#358359FF", "#9F4D23FF", "#D86C4FFF", "#170C2EFF",
                 "#473B75FF", "#F19C1FFF",
                 "#117733", "#DDCC77", "#CC6677", "#88CCEE",
                 "#44AA99", "#332288", "#AA4499", "#999933",
                 "#882255", "#661100", "#6699CC", "#888888")
    if (nrow(isomers) > length(palette)){
      set.seed(19580811)
      colors <- grDevices::colors()[grep('gr(a|e)y|white|light', grDevices::colors(), invert = T)]
      colors <- sample(colors, size = (nrow(isomers) - length(palette)))
      palette <- c(palette, colors)
    }

    for (p in 1:nrow(isomers)){
      grDevices::pdf(NULL) # use a pdf NULL device to save plots to an object
      grDevices::dev.control(displaylist = "enable")
      graphics::par(mar = c(3,4,4,1), mgp=c(2,1,0), bg = "white")
      startat <- 1
      start <- FALSE
      while (!start){
        if (nrow(eics[[startat]]) > 0){
          smt <- tryCatch({stats::predict(stats::smooth.spline(eics[[startat]]$RT,
                                                               eics[[startat]]$int,
                                                               spar=0.05),
                                          x = eics[[startat]]$RT)},
                          error = function(e) {return(list(x = eics[[startat]]$RT,
                                                           y = eics[[startat]]$int))})


          smt$y[smt$y < 0] <- 0

          plot(smt$x, smt$y, type = "l", lwd = 2, ylim = c(0, maxint),
               xlim = c(iniRT, endRT),
               col = scales::alpha("grey", 0.7),
               main = paste("ID: ", isomers$ID[p], "; ", fa, "; EIC: ",
                            round(isomers$mz[p],3), "; RT: ",
                            round(isomers$iniRT[p]), "-",
                            round(isomers$endRT[p]),
                            sep=""),
               xlab = "RT (sec)", ylab = "Intensity")
          start <- TRUE
        } else if (nrow(eics[[startat]]) == 0 & startat < length(msbatch$msobjects)) {
          startat <- startat + 1
        } else {
          if(verbose){cat("No peaks found")}
          plot(0, xlim = c(iniRT, endRT), ylim = c(0, 100), type = "l", lwd = 2,
               col = scales::alpha("grey", 0.7),
               main = paste(fa, "; EIC: ", round(isomers$mz[p],3), "; RT: ",
                            round(isomers$iniRT[p]), "-", round(isomers$endRT[p]),
                            sep=""),
               xlab = "RT (sec)", ylab = "Intensity")
          break
        }
      }
      if (length(eics) > startat){
        if (start){
          for (i in (startat+1):length(eics)){
            if (nrow(eics[[i]]) > 0){
              smt <- tryCatch({stats::predict(stats::smooth.spline(eics[[i]]$RT,
                                                                   eics[[i]]$int, spar=0.05),
                                              x = eics[[i]]$RT)},
                              error = function(e) {return(list(x = eics[[i]]$RT,
                                                               y = eics[[i]]$int))})

              smt$y[smt$y < 0] <- 0
              lines(smt$x, smt$y, type = "l", lwd = 2,
                    col = scales::alpha("grey", 0.7))
            }
          }
        }
      }

      for (i in 1:length(peaks[[p]])){
        if (nrow(peaks[[p]][[i]]) > 0){
          smt <- tryCatch({stats::predict(stats::smooth.spline(peaks[[p]][[i]]$RT,
                                                               peaks[[p]][[i]]$int,
                                                               spar=0.05),
                                          x = peaks[[p]][[i]]$RT)},
                          error = function(e) {return(list(x = peaks[[p]][[i]]$RT,
                                                           y = peaks[[p]][[i]]$int))})

          smt$y[smt$y < 0] <- NA
          lines(smt$x, smt$y, type = "l", lwd = 2,
                col = scales::alpha(palette[p], 0.7))
        }
      }
      plots[[counter]] <- grDevices::recordPlot() # save plot
      counter <- counter + 1
      invisible(grDevices::dev.off()) # close pdf NULL device
    }
  }

  return(plots)
}

# curateFAannotations
#' Modify FA annotations
#'
#' Modify FA annotations
#'
#' @param msbatch annotated msbatch.
#' @param faid data frame with 7 columns (ID, FAid, Adducts, mz, RT, iniRT and
#' endRT) containing curated FAs.
#' @param dmz mz tolerance in ppm.
#'
#' @description after FA annotation using \link{annotateFA}, the resulting
#' data frame can be modified to remove rows with unwanted annotation, iniRT and
#' endRT can be changed to redefine peak limits and extra rows may be written to
#' add new annotations. FAid should also be modified to contain unique names
#' such as "FA(16:1)n7" and "FA(16:1)n10" instead of generic "FA(16:1)". For
#' unknown fatty acids use FA(16:1)nx (nx, ny and nz are availables for all FA).
#'
#' Internal standards can also be added to normalize data later. Leave ID and
#' Adducts columns empty, write "IS" at the FAid column and add mz, RT, iniRT
#' and endRT information.
#'
#' @return annotated msbatch.
#'
#' @examples
#' \dontrun{
#' msbatch <- annotateFA(msbatch, dmz = 5)
#'
#' plots <- plotFA(msbatch, dmz = 10)
#'
#' pdf("FAs.pdf")
#' for (p in 1:length(plots)){
#'   print(plots[[p]])
#' }
#' dev.off()
#'
#' write.csv(msbatch$fas, file="faid.csv", row.names=FALSE)
#'
#' faid <- read.csv("faid_curated.csv", sep=",", dec=".")
#'
#' msbatch <- curateFAannotations(msbatch, faid)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
curateFAannotations <- function(msbatch, faid, dmz = 10){


  #============================================================================#
  # check arguments (por acabar)
  #============================================================================#
  faid$FAid <- gsub("[ \\.]", "", faid$FAid)

  msbatch$fas$iniRT <- round(msbatch$fas$iniRT, 4)
  msbatch$fas$endRT <- round(msbatch$fas$endRT, 4)
  faid$iniRT <- round(faid$iniRT, 4)
  faid$endRT <- round(faid$endRT, 4)
  rownames(msbatch$fas) <- msbatch$fas$ID


  if (any(duplicated(faid$FAid))){
    message(paste(unique(faid$FAid[duplicated(faid$FAid)]), collapse=", "), "duplicated.")
    stop("Compound names (FAid column) must be unique.")
  }
  if (any(duplicated(faid$ID) & !is.na(faid$ID))){
    stop("ID must be unique.")
  }

  faid <- faid[order(faid$ID),]
  #============================================================================#
  # remove FA
  #============================================================================#
  toremove <- msbatch$fas$ID[!msbatch$fas$ID %in% faid$ID]
  if (length(toremove) > 0){
    msbatch <- removeFA(msbatch, ids = toremove)
  }

  msbatch$fas <- msbatch$fas[unlist(sapply(msbatch$fas$ID, match, faid$ID)),]
  msbatch$fas <- msbatch$fas[!is.na(msbatch$fas$FAid),]
  #============================================================================#
  # change FA
  #============================================================================#
  tochange <- faid$ID[which(faid$iniRT[faid$ID %in% msbatch$fas$ID] != msbatch$fas$iniRT |
                      faid$endRT[faid$ID %in% msbatch$fas$ID] != msbatch$fas$endRT)]
  if (length(tochange) > 0){
    msbatch <- changeFArt(msbatch, id = tochange,
                          from = faid$iniRT[unlist(sapply(tochange, match, faid$ID))],
                          to = faid$endRT[unlist(sapply(tochange, match, faid$ID))])
  }


  #============================================================================#
  # add new FA
  #============================================================================#
  toadd <- which(is.na(faid$ID) & faid$FAid != "IS")
  if (length(toadd) > 0){
    msbatch <- addFA(msbatch, dmz=dmz, faid = faid$FAid[toadd],
                     adducts = faid$Adducts[toadd],
                     mz = faid$mz[toadd],
                     from = faid$iniRT[toadd],
                     to = faid$endRT[toadd])
  }

  # update FAid
  msbatch$fas$FAid[msbatch$fas$ID %in% faid$ID] <-
    faid$FAid[unlist(sapply(faid$ID[faid$ID %in% msbatch$fas$ID], match, msbatch$fas$ID))]

  msbatch$fas <- msbatch$fas[order(msbatch$fas$FAid),]

  #============================================================================#
  # search IS
  #============================================================================#
  if (any(faid$FAid == "IS")){
    is <- which(faid$FAid == "IS")
    msbatch <- searchIS(msbatch, mz = faid$mz[is], rt = faid$RT[is],
                        minRT = faid$iniRT[is], maxRT = faid$endRT[is],
                        dmz = dmz)
  }

  return(msbatch)
}

# remove incorrect FA
#' Remove incorrect FA annotations
#'
#' Remove incorrect FA annotations
#'
#' @param msbatch annotated msbatch.
#' @param ids integer vector specifying FA ids to be removed.
#'
#' @return annotated msbatch.
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
removeFA <- function(msbatch, ids){

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
  # check FA info
  if (missing(ids)){
    stop("FA IDs are required to remove them from annotation results")
  }

  #============================================================================#
  # Remove FA
  #============================================================================#
  msbatch$fas <- msbatch$fas[!msbatch$fas$ID %in% ids,]

  return(msbatch)
}

# changeFArt
#' Modify rt peak limits of annotated FAs
#'
#' Modify rt peak limits of annotated FAs
#'
#' @param msbatch annotated msbatch.
#' @param id integer vector specifying FA ids to be modified.
#' @param from numeric vector specifying the peak start.
#' @param to numeric vector specifying the peak end.
#'
#' @return annotated msbatch.
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
changeFArt <- function(msbatch, id, from, to){

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
  # check FA info
  if (missing(id)){
    stop("FA ID is required")
  }
  if (missing(from)){
    stop("from argument is required. When does the peak start?")
  }
  if (missing(to)){
    stop("to argument is required. When does the peak end?")
  }
  if (length(id) != length(from) || length(id) != length(to)){
    stop(" id, from and to vectors must have the same length")
  }

  #============================================================================#
  # Change FA RTs
  #============================================================================#
  for (i in 1:length(id)){
    msbatch$fas$iniRT[msbatch$fas$ID == id[i]] <- from[i]
    msbatch$fas$endRT[msbatch$fas$ID == id[i]] <- to[i]
  }

  return(msbatch)
}

# addFA
#' Add missing FA annotations
#'
#' Add missing FA annotations
#'
#' @param msbatch annotated msbatch.
#' @param dmz mz tolerance in ppm.
#' @param faid character vector specifying FA names (i.e. "FA(16:1)").
#' @param adducts character vector specifying adducts.
#' @param mz numeric vector specifying FA mz.
#' @param from numeric vector specifying the peak start.
#' @param to numeric vector specifying the peak end.
#'
#' @return annotated msbatch.
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
addFA <- function(msbatch, dmz = 5, faid, adducts = "M-H", mz, from, to){

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
  # check FA info
  if (missing(faid)){
    stop("FA ID is required")
  }
  if (length(adducts) == 1){
    adducts <- rep(adducts, length(faid))
  }
  if (missing(mz)){
    stop("mz is required")
  }
  if (missing(from)){
    stop("from argument is required. When does the peak start?")
  }
  if (missing(to)){
    stop("to argument is required. When does the peak end?")
  }
  if (length(faid) != length(from) || length(faid) != length(to)){
    stop("faid, adducts, from and to vectors must have the same length")
  }

  #============================================================================#
  # add FA
  #============================================================================#
  idstart <- max(msbatch$fas$ID)+1
  id <- idstart
  for (f in 1:length(faid)){
    newrow <- c(ID = id,
                FAid = faid[f],
                Adducts = adducts[f],
                mz = mz[f],
                RT = mean(from[f], to[f]),
                iniRT = from[f],
                endRT = to[f])
    msbatch$fas <- rbind(msbatch$fas, newrow)
    id <- id + 1
  }
  # tonumeric <- which(!colnames(msbatch$fas) %in% c("FAid", "Adducts", "isotope"))
  tonumeric <- which(!colnames(msbatch$fas) %in% c("FAid", "Adducts"))
  msbatch$fas[,tonumeric] <- apply(msbatch$fas[,tonumeric], 2, as.numeric)

  return(msbatch)
}

# searchIS
#' Search internal stardards.
#'
#' Search internal stardards.
#'
#' @param msbatch annotated msbatch.
#' @param mz numeric vector specifying IS mz.
#' @param rt numeric vector specifying IS rt.
#' @param minRT numeric vector specifying lower limits for IS rt.
#' @param maxRT numeric vector specifying upper limits for IS rt.
#' @param dmz mz tolerance in ppm.
#'
#' @return annotated msbatch.
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
searchIS <- function(msbatch, mz, rt, minRT, maxRT, dmz = 10){

  if(length(mz) != length(rt)){
    stop("mz and rt vectors must have the same length")
  }

  is <- list()
  for (i in 1:length(mz)){
    matches <- mzMatch(mz[i], msbatch$features$mz, ppm = dmz)
    matches <- matches[seq(1, length(matches), 2)]
    rts <- msbatch$features$RT[matches]
    if (length(rts) > 0){
      rts <- rts[rts <= maxRT[i] & rts >= minRT[i]]
    }

    if (length(rts) > 0){
      matches <- matches[which.min(abs(rts - rt[i]))]
      is[[i]] <- as.numeric(msbatch$features[matches, colnames(msbatch$features) %in%
                                               make.names(msbatch$metaData$sample)])
      is[[i]][is[[i]] == 0] <- 1
    } else {
      is[[i]] <- rep(1, nrow(msbatch$metaData))
    }
  }
  is <- data.frame(do.call(rbind, is))
  colnames(is) <- make.names(msbatch$metaData$sample)
  rownames(is) <- make.names(mz)
  msbatch$IS <- is

  return(msbatch)
}

# searchFAisotopes
#' Search FA isotopes
#'
#' Search FA isotopes
#'
#' @param msbatch annotated msbatch.
#' @param dmz mz tolerance in ppm.
#' @param coelCutoff coelution score threshold between parent and isotope peaks.
#'
#' @return fadata list.
#'
#' @examples
#' \dontrun{
#' fadata <- searchFAisotopes(msbatch, dmz = 10, coelCutoff = 0.4)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
searchFAisotopes <- function(msbatch,
                             dmz = 5,
                             coelCutoff = 0.7){

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
  # check FA info

  #============================================================================#
  # Change msbatch structure to make it compatible with LipidMS::searchIsotopes
  #============================================================================#
  features <- msbatch$features
  msbatch$features <- msbatch$fas
  msbatch$features$LipidMSid <- gsub("n[a-z]|n[0-9]*", "", msbatch$features$FAid)
  intmatrix <- data.frame(matrix(NA, nrow=nrow(msbatch$fas),
                                 ncol=length(make.names(msbatch$metaData$sample))))
  colnames(intmatrix) <- make.names(msbatch$metaData$sample)
  msbatch$features <- data.frame(cbind(msbatch$features, intmatrix))
  msbatch$features$minmz <- msbatch$features$maxmz <- msbatch$features$minRT <-
    msbatch$features$maxRT <- msbatch$features$npeaks <- msbatch$features$group <-
    msbatch$features$isotope <- msbatch$features$isoGroup <- NA
  msbatch$features <- msbatch$features[,c("ID", "FAid", "LipidMSid", "Adducts", colnames(features))]

  #============================================================================#
  # Search isotopes for all identified FA
  #============================================================================#
  isotopes <- searchIsotopesmsbatch(msbatch,
                                    label = "13C",
                                    adductsTable = LipidMS::adductsTable,
                                    ppm = dmz,
                                    coelCutoff = coelCutoff)
  msbatch$features <- msbatch$features[,c("ID", "FAid", "Adducts", "mz", "RT",
                                          "iniRT", "endRT")]

  # add identifiers to isotopes
  ids <- msbatch$features$ID
  m0 <- which(isotopes$Isotope == "[M+0]")
  m0 <- c(m0, nrow(isotopes)+1)

  newids <- rep(0, nrow(isotopes))
  for (r in 1:(length(m0)-1)){
    newids[m0[r]:(m0[r+1]-1)] <- ids[r]
  }
  isotopes <- data.frame(ID = newids, isotopes)

  isotopes <- isotopes[order(isotopes$ID, isotopes$mz),]

  #============================================================================#
  # Save results
  #============================================================================#
  # undo structure changes
  msbatch$fas <- msbatch$features
  msbatch$features <- features
  # save results
  isotopes$LipidMSid <- msbatch$fas$FAid[unlist(sapply(isotopes$ID, match, msbatch$fas$ID))]
  colnames(isotopes)[2] <- "FAid"
  msbatch$isotopes <- isotopes

  #============================================================================#
  # get fadata list
  #============================================================================#
  fadata <- msbatch2fadata(msbatch, msbatch$fas)

  return(fadata)
}

# readfadatafile
#' read FA data from a csv file.
#'
#' read FA data from a csv file.
#'
#' @param file csv file name.
#' @param sep column delimiter.
#' @param dec character used for decimal points.
#'
#' @description First rows must contain metadata information such
#' as sample groups (row named sampletype) and any other extra information like
#' protein levels for external normalization. Then, the following row must
#' contain a Compound and Label columns followed by all sample names.
#' FA names must be unique and omega series must be indicated
#' (i.e. FA(20:4)n3, FA(24:1)n9, FA(16:0)). Unknown FA series can be named as nx,
#' ny, nz to differentiate between isomers. Labels must be specified with integer
#' numbers for 0 to maximum number of carbons.
#'
#' @return fadata.
#'
#' @examples
#' \dontrun{
#' fadata <- readfadatafile("externafadata.csv", sep=",", dec=".")
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
readfadatafile <- function(file, sep = ",", dec = "."){

  # read fadata file
  if (tools::file_ext(file) != "csv"){
    stop("Input file must be a csv")
  }
  data <- read.csv(file, header = FALSE, sep = sep, dec = dec)
  if (!all(c("sampletype", "Label") %in% data[,2])){
    stop("Please check input format")
  }

  # read metadata
  toptable <- 1:(which(data[,2] == "Label")-1)
  metadata <-  t(data[toptable, 2:ncol(data), drop=FALSE])
  colnames(metadata) <- metadata[1,]
  metadata <- metadata[2:nrow(metadata),, drop=FALSE]

  # Clean input data
  cnames <- data[which(data[,2] == "Label"),]
  if (any(!cnames[1:2] %in% c("Compound", "Label"))){
    stop("Two first columns must be Compound and Label")
  }
  if (nrow(metadata) != ncol(data) - 2){
    stop("Check groups")
  }

  # read data
  data <- data[(which(data[,2] == "Label")+1):nrow(data), , drop = FALSE]
  colnames(data) <- as.character(cnames)

  fattyacids <- data[,c("Compound", "Label")]
  fattyacids <- fattyacids[!apply(fattyacids, 1, function(x) all(x == "")),]
  fattyacids$Label <- as.numeric(fattyacids$Label)
  data <- data[,3:ncol(data) , drop = FALSE]
  data <- data[!apply(data, 1, function(x) all(x == "")),, drop=FALSE]
  data <- data.frame(apply(data, 2, as.numeric))

  data <- data[order(fattyacids$Compound, fattyacids$Label),, drop = FALSE]
  fattyacids <- fattyacids[order(fattyacids$Compound, fattyacids$Label),]
  rownames(data) <- paste(fattyacids$Compound, "_M+", fattyacids$Label, sep="")

  metadata <- data.frame(sample = colnames(data), metadata)

  # extract IS info
  if(any(fattyacids$Compound == "IS")){
    IS <- data[fattyacids$Compound == "IS",, drop=FALSE]
  } else {
    IS <- t(data.frame(IS = rep(1, nrow(metadata))))
  }
  data <- data[fattyacids$Compound != "IS", , drop=FALSE]
  fattyacids <- fattyacids[fattyacids$Compound != "IS", , drop=FALSE]


  # tidy data
  M <- as.numeric(sapply(unique(fattyacids$Compound), function(x){
    unlist(strsplit(x, "[():]"))[2]
  }))

  compounds <- unique(fattyacids$Compound)
  newdata <- c()
  newfattyacids <- c()
  for (c in 1:length(compounds)){
    comp <- data[fattyacids$Compound == compounds[c], , drop = FALSE]
    ss <- sapply(0:M[c], match, fattyacids$Label[fattyacids$Compound == compounds[c]])
    ss <- unlist(lapply(ss, function(x) if(length(x) == 0){NA}else{x}))
    comp <- comp[ss, , drop = FALSE]
    comp[is.na(comp)] <- 0
    mdata <- data.frame(Compound = compounds[c], Label = 0:M[c])
    newdata <- rbind(newdata, comp)
    newfattyacids <- rbind(newfattyacids, mdata)
  }

  colnames(IS) <- colnames(newdata)
  # create fadata
  fadata <- list(metadata = metadata, fattyacids = newfattyacids, IS = IS,
                 intensities = newdata)
  fadata$metadata <- fadata$metadata[order(colnames(fadata$intensities)),]
  if ("IS" %in% names(fadata)){
    fadata$IS <- fadata$IS[,order(colnames(fadata$intensities)), drop=FALSE]
  }
  fadata$intensities <- fadata$intensities[order(fadata$fattyacids$Compound,
                                                 fadata$fattyacids$Label),
                             order(colnames(fadata$intensities)), drop = FALSE]
  fadata$fattyacids <- fadata$fattyacids[order(fadata$fattyacids$Compound,
                                               fadata$fattyacids$Label),, drop = FALSE]

  return(fadata)
}

# dataCorrection
#' Data correction for natural abundance of 13C and data normalization using
#' internal standards followed by blank substraction.
#'
#' Data correction for natural abundance of 13C and data normalization using
#' internal standards followed by blank substraction.
#'
#' @param fadata fadata list.
#' @param correct13C logical. If TRUE, data is corrected for natural abundance
#' of 13C. Set to FALSE if data has been already been corrected.
#' @param blankgroup name used to define blank samples group.
#' @param externalnormalization column name at the metadata data frame of any
#' additional measure that must be used to normalize data (i.e. protein).
#' @param resolution resolution of the mass spectrometer.
#' @param purity13C purity of the tracer employed.
#' @param verbose print information messages.
#'
#' @return corrected fadata.
#' 
#' @references Su X, Lu W, Rabinowitz J (2017). Metabolite Spectral Accuracy on 
#' Orbitraps. Analytical Chemistry, 89(11), 5940-5948, PMID: 28471646, R package 
#' version 0.2.4 (2021), <https://doi.org/10.1021/acs.analchem.7b00396>. 
#'
#' @examples
#' \donttest{
#' ssdata <- dataCorrection(ssexamplefadata, blankgroup="Blank")
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
dataCorrection <-function(fadata,
                          correct13C = TRUE,
                          blankgroup = "blank", 
                          externalnormalization = c(),
                          resolution = 140000, 
                          purity13C = 0.99,
                          verbose = TRUE){

  #============================================================================#
  # Check arguments
  #============================================================================#
  if (!is.list(fadata) | length(fadata) < 3 |
      !all(c("metadata", "fattyacids", "intensities") %in% names(fadata))){
    stop("fadata must be a list with at least 3 elements: metadata, fattyacids
         (data.frame with two columns: Compound and Label) and intensities.")
  }

  #============================================================================#
  # Correct data
  #============================================================================#
  if (correct13C){
    if(verbose){cat("Data correction for natural abundance of 13C...")}
    fadata <- correctNatAb13C(fadata, resolution = resolution, purity = purity13C)
    if(verbose){cat("OK\n")}
  }
  if ("IS" %in% names(fadata)){
    if(verbose){cat("Data normalization with internal standards...")}
    fadata <- normalizeIS(fadata, verbose = verbose)
    if(verbose){cat("OK\n")}
  }
  if (length(blankgroup) > 0){
    if(verbose){cat("Blank substraction...")}
    fadata <- blankSubstraction(fadata, blankgroup = blankgroup, verbose = verbose)
    if(verbose){cat("OK\n")}
  }
  if (length(externalnormalization) > 0){
    if(verbose){cat("External normalization...")}
    fadata <- externalNormalization(fadata, 
                                    externalnormalization = externalnormalization,
                                    verbose = verbose)
    if(verbose){cat("OK\n")}
  }
  

  # update fadata
  falist <- lapply(unique(fadata$fattyacids$Compound), function(x)
    fadata$intensities[fadata$fattyacids$Compound==x, ,drop=FALSE])
  mid <- do.call(rbind, lapply(falist, function(x) apply(x, 2, function(y) y/sum(y))))
  rownames(mid) <- paste(fadata$fattyacids$Compound, "_M+", fadata$fattyacids$Label, sep="")
  poolsize <- do.call(rbind, lapply(falist, function(x) apply(x, 2, sum)))
  rownames(poolsize) <- unique(fadata$fattyacids$Compound)
  relativepoolsize <- apply(poolsize, 2, function(x) x/sum(x))
  if (nrow(poolsize) == 1){
    relativepoolsize <- t(as.data.frame(relativepoolsize))
  }
  rownames(relativepoolsize) <- rownames(poolsize)
  
  fadata$mid <- mid
  fadata$poolsize <- poolsize
  fadata$relativepoolsize <- relativepoolsize

  return(fadata)
}


# normalizeIS
#' Data normalization using internal stardards.
#'
#' Data normalization using internal stardards.
#'
#' @param fadata fadata list.
#' @param verbose print information messages.
#'
#' @return normalised fadata by IS.
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
normalizeIS <-function(fadata, verbose = TRUE){

  #============================================================================#
  # Check arguments
  #============================================================================#
  if (!is.list(fadata) | length(fadata) < 3 |
      !all(c("metadata", "fattyacids", "intensities") %in% names(fadata))){
    stop("fadata must be a list with at least 3 elements: metadata, fattyacids
         (data.frame with two columns: Compound and Label) and intensities.")
  }

  if (!"IS" %in% names(fadata)){
    if(verbose){cat("No IS data available.")}
  } else {
    if (nrow(fadata$IS) == 0){
      if(verbose){cat("No IS data available.")}
    }
    if (all(fadata$IS[1,] == 1)){
      if(verbose){cat("No IS data available.")}
    }
  }


  #============================================================================#
  # Normalize by IS
  #============================================================================#
  if ("IS" %in% names(fadata)){
    if (nrow(fadata$IS) > 0){
      fadata$intensities <- sweep(fadata$intensities, 2,
                                  apply(fadata$IS, 2, prod), `/`)
      fadata$IS <- sweep(fadata$IS, 2,
                         apply(fadata$IS, 2, prod), `/`)
    }
  }

  return(fadata)
}

# correctNatAb13C
#' correct data for natural abundance of 13C using accucor algorithm.
#'
#' correct data for natural abundance of 13C using accucor algorithm.
#'
#' @param fadata fadata.
#' @param resolution resolution of the mass spectrometer.
#' @param purity purity of the tracer employed.
#'
#' @return corrected fadata.
#'
#' @references Su X, Lu W, Rabinowitz J (2017). “Metabolite Spectral Accuracy on 
#' Orbitraps.” Analytical Chemistry, 89(11), 5940-5948, PMID: 28471646, R package 
#' version 0.2.4 (2021), <https://doi.org/10.1021/acs.analchem.7b00396>.
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
correctNatAb13C <- function(fadata, resolution = 140000, purity = 0.99){

  #============================================================================#
  # Check arguments
  #============================================================================#
  if (!is.list(fadata) | length(fadata) < 3 |
      !all(c("metadata", "fattyacids", "intensities") %in% names(fadata))){
    stop("fadata must be a list with at least 3 elements: metadata, fattyacids
         (data.frame with two columns: Compound and Label) and intensities.")
  }

  #============================================================================#
  # Prepare data: order each compound by label and get composition
  #============================================================================#
  fadata$intensities <- fadata$intensities[order(fadata$fattyacids$Compound,
                                   fadata$fattyacids$Label),, drop = FALSE]

  fadata$fattyacids <- fadata$fattyacids[order(fadata$fattyacids$Compound,
                                           fadata$fattyacids$Label),, drop = FALSE]

  fas <- fadata$fattyacids$Compound
  fas <- data.frame(t(sapply(fas, function(x){
    x <- unlist(strsplit(x, "[()]"))
    class <- x[1]
    cdb <- x[2]
    return(c(class, cdb))
  })))
  colnames(fas) <- c("Class", "CDB")
  formulas <- do.call(rbind, apply(fas, 1, getFormula,
                                   dbs = list(fadb = FAMetA::fattyacidsdb)))[,1]

  #============================================================================#
  # Correct natural abundance of 13C using accucor package
  #============================================================================#
  CorrMatrix <- matrix(0, ncol=ncol(fadata$intensities), nrow=0)
  for (f in unique(fadata$fattyacids$Compound)){

    ss <- which(fadata$fattyacids$Compound == f)
    form <- unique(formulas[ss[order(fadata$fattyacids$Label[ss])]])
    fa <- fadata$intensities[ss[order(fadata$fattyacids$Label[ss])],, drop = FALSE]
    label <- fadata$fattyacids$Label[ss[order(fadata$fattyacids$Label[ss])]]
    corrected <- accucor::carbon_isotope_correction(formula = form,
                                                     datamatrix = as.matrix(fa),
                                                     label = label,
                                                     Resolution = resolution,
                                                     purity = purity,
                                                     ReportPoolSize = FALSE)
    CorrMatrix <- rbind(CorrMatrix, corrected)
  }
  CorrMatrix <- data.frame(CorrMatrix)
  colnames(CorrMatrix) <- colnames(fadata$intensities)

  # save results
  fadata$intensities <- CorrMatrix

  return(fadata)
}

# blankSubstraction
#' substract blank samples.
#'
#' substract blank samples.
#'
#' @param fadata fadata.
#' @param blankgroup name used to define blank samples group.
#' @param verbose print information messages.
#'
#' @return blank substracted fadata.
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
blankSubstraction <- function(fadata, blankgroup = "blank", verbose = TRUE){

  #============================================================================#
  # Check arguments
  #============================================================================#
  if (!is.list(fadata) | length(fadata) < 3 |
      !all(c("metadata", "fattyacids", "intensities") %in% names(fadata))){
    stop("fadata must be a list with at least 3 elements: metadata, fattyacids
         (data.frame with two columns: Compound and Label) and intensities.")
  }

  #============================================================================#
  # Substract mean of blank samples intensities and remove them
  #============================================================================#
  if (any(toupper(fadata$metadata$sampletype) == toupper(blankgroup))){
    blank <- apply(fadata$intensities[,toupper(fadata$metadata$sampletype) == toupper(blankgroup), drop=FALSE],
                   1, mean, na.rm = T)
    fadata$intensities <- fadata$intensities - blank
    fadata$intensities[fadata$intensities < 0] <- 0
    fadata$intensities <- fadata$intensities[,toupper(fadata$metadata$sampletype) != toupper(blankgroup),
                               drop = FALSE]
    fadata$IS <- fadata$IS[,toupper(fadata$metadata$sampletype) != toupper(blankgroup),
                                             drop = FALSE]
    fadata$metadata <- fadata$metadata[toupper(fadata$metadata$sampletype) != toupper(blankgroup),,
                                   drop = FALSE]
  } else {
    if(verbose){cat("No blank samples to substract. Use blankgroup argument to 
    indicate the group name for blank samples.")}
  }

  return(fadata)
}

# externalNormalization
#' External normalization using additional measures (i.e. protein levels).
#'
#' External normalization using additional measures (i.e. protein levels).
#'
#' @param fadata fadata list.
#' @param externalnormalization column names of metadata data frame used to
#' define external measures.
#' @param verbose print information messages.
#'
#' @return normalised fadata by external measures.
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
externalNormalization <-function(fadata, externalnormalization, verbose = TRUE){

  #============================================================================#
  # Check arguments
  #============================================================================#
  if (!is.list(fadata) | length(fadata) < 3 |
      !all(c("metadata", "fattyacids", "intensities") %in% names(fadata))){
    stop("fadata must be a list with at least 3 elements: metadata, fattyacids
         (data.frame with two columns: Compound and Label) and intensities.")
  }

  if (length(externalnormalization) == 0){
    if(verbose){cat("No external measures available")}
  } else {
    fadata$metadata[,externalnormalization] <-
      apply(fadata$metadata[,externalnormalization, drop=FALSE], 2, as.numeric)
    vec <- apply(fadata$metadata[,externalnormalization, drop=FALSE], 1, prod)
    fadata$intensities <- sweep(fadata$intensities, 2, vec, `/`)
  }

  return(fadata)
}

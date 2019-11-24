#' Fill in vector sequence
#'
#' Takes a vector of integers provided by segment.cpp and returns a list of
#' integer sequences. The start/end points of each sequence are the values of 
#' the inputed integer vector.
#' of integers 
#'
#' @param x Integer. A vector of integers.
#' 
#' @param overlap Integer. How much overlap to add to each sequence
#'
#' @return A list of integer vectors
#'
segment_index <- function(x, overlap=7)
{
  if (min(x) == 0)
    x <- x + 1L
  from <- pmax.int(x[-length(x)] - overlap, min(x))
  to <- pmin.int(x[-1] + overlap, max(x))
  Map(seq, from=from, to=to)
}



#' Merge list of sequences that contain too few observations
#'
#' Takes a list of sequences and merges adjacent sequences if one or more has
#' fewer than a specified minimum number of entries.
#'
#' @param lst List. A list of sequences.
#' 
#' @param mindays Integer. The lower threshhold. Sequences with fewer than the
#' threshold are merged withand adjacent sequence in the list.
#' 
#' @param combine Character. Combine offending list with preceding ('last') or
#' the following ('next') sequence in the list. This doesn't apply to the first
#' or last sequence in the list for obvious reasons. 
#' 
#' @return A list of integer sequences
#'
regroup_segments <- function(lst, mindays=30, combine=c("last", "next"))
{
  combine <- match.arg(combine)
  neighbor <- switch(combine, "last"=`-`, "next"=`+`)
  newlst <- lst
  nlst <- length(lst)
  nx <- vapply(lst, length, integer(1L))
  idx <- which(nx < mindays)
  if (length(idx) == 0)
    return(lst)
  midx <- min(idx)
  if (midx == 1) {
    newlst[[2]] <- unique(c(lst[[1]], lst[[2]]))
    newlst[[1]] <- NULL
  } else if (midx == nlst) {
    newlst[[nlst - 1]] <- unique(c(lst[[nlst - 1]], lst[[nlst]]))
    newlst[[nlst]] <- NULL
  } else {
    jj <- neighbor(midx, 1)
    newlst[[jj]] <- sort(unique(c(lst[[midx]], lst[[jj]])))
    newlst[[midx]] <- NULL
  }
  regroup_segments(newlst, mindays, combine)
}



#' Perform parametric outlier detection for count data
#'
#' This function fits three different parametric count data distributions: poisson,
#' negative binomial, and geometric. Selection between the three is carried out
#' by choosing which returned the largest maximum likelihood.
#'
#' @param x Integer. A vector of integers.
#' 
#' @param sigma Numeric. The significance level for outlier detection.
#' 
#' @return Index vector of detected outliers
#'
#' @importFrom MASS fitdistr
#'
#' @importFrom stats ppois pnbinom pgeom
#' 
count_dist_outlier <- function(x, sigma=0.01)
{
  sigma <- sigma/2
  fit <- list(fitdistr(x, "poisson"), fitdistr(x, "negative binomial"),
              fitdistr(x, "geometric"))
  idx <- which.max(sapply(fit, `[[`, 'loglik'))
  fit <- fit[[idx]]
  rx <- min(x):max(x)
  if (idx==1) {
    dtf <- data.frame(rx, F=ppois(rx, lambda=fit$estimate[1]))
  } else if (idx==2) {
    dtf <- data.frame(rx, F=pnbinom(rx, size=fit$estimate['size'], mu=fit$estimate['mu']))
  } else {
    dtf <- data.frame(rx, F=pgeom(rx, prob=fit$estimate['prob']))
  }
  low <- dtf[dtf$F < sigma, "rx"]
  if (length(low) == 0) {
    low <- -Inf
  } else {
    low <- max(low)
  }
  high <- dtf[dtf$F > 1 - sigma, "rx"]
  if (length(high) == 0) {
    high <- Inf
  } else {
    high <- min(high)
  }
  which(x <= low | x >= high)
}



#' Perform Isolation Forest outlier detection
#'
#' Run an isolation forest 
#'
#' @param x Integer. A vector of integers.
#' 
#' @param STL Logical. Should the series be seasonally detrended before isolation
#' forest is run.
#' 
#' @param ccF Numeric. The significance level for outlier detection. See Details.
#' 
#' @details The result of the isolation forest assigns an anomaly score for each
#' observations. These observations are transformed by a z-score and then those
#' z-scores above \code{ccF} are considered outliers.
#' 
#' @return Index vector of detected outliers
#'
#' @importFrom isofor iForest
#'
#' @importFrom stats stl sd ts 
#' 
iso_forest <- function(x, STL=TRUE, ccF=2.5)
{
  N <- length(x)
  if (STL) {
    series <- stl(ts(x, frequency=7), "periodic")
    x <- series$time.series[,3,drop=F]
  }
  isof <- iForest(cbind(1, x), nt=200, phi=4*round(N/5))
  p <- predict(isof, cbind(1, x))
  zp <- (p - mean(p)) / sd(p)
  bool <- zp > ccF
  which(bool)
}




#' Perform mean-shift detection
#'
#' Perform mean-shift detection using a rolling window
#'
#' @param x Integer. A vector of integers.
#' 
#' @param width Integer. Size of window
#' 
#' @param ccB Numeric. Threshold level for mean-shift detection.
#' 
#' @param idx.outlier Integer. A vector of series indexes specifying outliers to
#' control for before applying algorithm.
#' 
#' @details Two adjacent blocks (each of size \code{width}) of the series create
#' a rolling window. At each index, we test for statistically significant
#' differences in the mean of the two series blocks.
#' 
#' @return Index vector of detected outliers
#' 
#' @seealso Algorithm details on Hammer at ~/Documents/aF_products/ts_anomaly
#' 
#' @importFrom zoo rollapply
#'
#' @importFrom stats var
#;
block_average <- function(x, width=14, ccB = 2, idx.outlier=NULL)
{
  nx <- length(x)
  if (!is.null(idx.outlier)) {
    weighted_mean <- function(ii)
    {
      jj <- min(ii + width, nx)
      mean(x[seq(ii + 1, jj)])
    }
    x[idx.outlier] <- vapply(idx.outlier, weighted_mean, numeric(1L))
  }
  rollsum <- zoo::rollapply(x, width=width, sum)
  rmn <- length(rollsum)
  SUM <- cbind(rollsum[1:(rmn - width)], rollsum[-(1:width)])
  
  Dmu <- (SUM[,2] - SUM[,1]) / width
  avar <- var(Dmu)
  bool <- abs(Dmu) > (ccB * sqrt(avar))
  return(which(c(rep(F, width), bool, rep(F, width-1))))
}



#' Perform jump detection
#'
#' Perform jump detection using a sequential average algorithm
#'
#' @param x Integer. A vector of integers.
#' 
#' @param width Integer. Size of window
#' 
#' @param p Numeric. Threshold level for jump detection.
#' 
#' @param idx.outlier Integer. A vector of series indexes specifying outliers to
#' control for before applying algorithm.
#' 
#' @return Index vector of detected outliers
#' 
#' @seealso Algorithm details on Hammer at ~/Documents/aF_products/ts_anomaly
#' 
#' @importFrom zoo rollapply
#'
#' @importFrom stats var qt
#'
seq_average <- function(x, width=14, p=0.01, idx.outlier=NULL)
{
  nx <- length(x)
  get_idx_block <- function(ii, reverse=FALSE)
  {
    if (reverse) {
      jj <- max(ii - width + 1L, 1L)
      return(jj:ii)
    } else {
      jj <- min(ii + width - 1L, nx)
      return(ii:jj)
    }
  }
  if (!is.null(idx.outlier)) {
    weighted_mean <- function(ii)
    {
      jj <- min(ii + width, nx)
      mean(x[seq(ii + 1, jj)])
    }
    x[idx.outlier] <- vapply(idx.outlier, weighted_mean, numeric(1L))
  }
  rollsum <- zoo::rollapply(x, width=width, sum)
  rollvar <- zoo::rollapply(x, width=width, var)
  rollmean <- rollsum / width
  
  # critical t-value
  tval <- qt((1 - p/2), 2 * (width - 1))
  # critical difference
  D <- tval * sqrt((2 * rollvar) / width)
  # mean of first block
  
  # Test all points agains their preceeding block's average
  chgs <- x[-(1:width)] - rollmean[-(nx-width)]
  checkidx <- which(abs(chgs) > D[-length(D)]) + width
  # Starting at each of those points, run RSI 
  
  RSI <- function(ii)
  {
    sign_bool <- function(x) as.logical((sign(x) + 1) / 2)
    upshift <- sign_bool(chgs[ii-width])
    xx <- x[get_idx_block(ii)]
    DD <- D[ii-width]
    xstar <- if (upshift) {
      R2 <- rollmean[ii-width] + DD
      xx - R2
    } else {
      R2 <- rollmean[ii-width] - DD
      R2 - xx
    }
    rsi <- cumsum(xstar) / (width * sqrt(rollvar[ii-width]))
    all(rsi > 0)
  }
  rsi <- vapply(checkidx, RSI, logical(1L))
  if (all(!rsi))
    return(NULL)
  else
    return(checkidx[rsi])
}



#' Perform jump detection and mean-shift detection
#'
#' @param x Integer. A vector of integers.
#' 
#' @param width Integer. Size of window
#' 
#' @param ccB Numeric. Threshold level for mean-shift detection.
#' 
#' @param p Numeric. Threshold level for jump detection.
#' 
#' @param idx.outliers Integer. A vector of series indexes specifying outliers to
#' control for before applying algorithm.
#' 
#' @return A list containing both outlier indicies from \code{\link{seq_average}}
#' and \code{\link{block_average}}.
#' 
#' @seealso Algorithm details on Hammer at ~/Documents/aF_products/ts_anomaly
#' 
#' @importFrom zoo rollapply
#' 
detect_meanshift <- function(x, width=7, ccB=2, p=0.01, idx.outliers=NULL)
{
  stopifnot(width > 0)
  if (width < 1) {
    width <- round(width * length(x))
  } else {
    width <- round(width)
  }
  rollsum <- zoo::rollapply(x, width=width, sum)
  idxBLKAVG <- block_average(x, width, ccB, idx.outliers)
  idxSEQAVG <- seq_average(x, width, p, idx.outliers)
  return(list(idxBLKAVG=idxBLKAVG, idxSEQAVG=idxSEQAVG))
}



#' Performs outlier detection
#'
#' @param x Integer. A vector of integers.
#' 
#' @param STL Logical. Should the series be seasonally detrended before isolation
#' forest is run.
#' 
#' @param sigma Numeric. The significance level for outlier detection.
#' 
#' @param ccF Numeric. The significance level for outlier detection. See Details.
#'  
#' @return A list containing both outlier indicies from \code{\link{iso_forest}}
#' and \code{\link{count_dist_outlier}}.
#'
detect_outliers <- function(x, STL=FALSE, sigma=0.01, ccF=2.5)
{
  idxDIST <- count_dist_outlier(x, sigma=sigma)
  idxISOF <- iso_forest(x, STL, ccF)
  return(list(idxDIST=idxDIST, idxISOF=idxISOF))
}



#' Performs variance shift detection
#'
#' @param x Integer. A vector of integers.
#' 
#' @param width Integer. Value need to control for outliers passed to function by
#' \code{idx.outliers}.
#' 
#' @param var.diff Logical. Should the series be differenced before algorithm is
#' run.
#' 
#' @param idx.outliers Integer. A vector of series indexes specifying outliers to
#' control for before applying algorithm.
#'  
#' @return A list containing index of variance shifts.
#' 
#' @importFrom changepoint cpt.var
#'
detect_varshift <- function(x, width=7, var.diff=FALSE, idx.outliers=NULL)
{
  nx <- length(x)
  if (!is.null(idx.outliers)) {
    weighted_mean <- function(ii)
    {
      uw <- floor(width / 2)
      lw <- uw + as.integer(width %% 2 == 1)
      jj <- max(ii - lw, 1)
      kk <- min(ii + uw, nx)
      mean(x[setdiff(seq(jj, kk), ii)])
    }
    x[idx.outliers] <- vapply(idx.outliers, weighted_mean, numeric(1L))
  }
  idxVAR <- if (var.diff) {
    changepoint::cpt.var(diff(x), method="PELT")@cpts + 1
  } else {
    changepoint::cpt.var(x, method="PELT")@cpts
  }
  return(list(idxVAR=idxVAR))
}



#' Translates anomalies in R lists to json
#'
#' @param lst List. Output from \code{\link{segment_anomaly}} or
#' \code{\link{anomaly_detection}}
#' 
#' @return A list containing index of variance shifts.
#' 
#' @importFrom jsonlite toJSON
#'
merge_anomaly <- function(lst)
{
  dates <- seq(as.Date(lst$start_date), as.Date(lst$end_date), by="day")
  outliers <- Reduce(intersect, lst[c("idxDIST", "idxISOF")])
  meanshift <- Reduce(intersect, lst[c("idxBLKAVG", "idxSEQAVG")])
  varshift <- lst[["idxVAR"]]
  uidx <- sort(unique(c(outliers, meanshift, varshift)))
  anomaly_type <- function(j)
  {
    out <- c()
    if (j %in% outliers)
      out <- c(out, "outlier")
    if (j %in% meanshift)
      out <- c(out, "meanshift")
    if (j %in% varshift)
      out <- c(out, "varshift")
    out
  }
  outjson <- lapply(uidx, anomaly_type)
  names(outjson) <- dates[uidx]
  return(outjson)
}



#' Perform Anomaly Detection
#'
#' @param x Integer. A vector of integers.
#' 
#' @param start_date Character. Date of first observation in the series \code{x}.
#' 
#' @param end_date Character. Date of last observation in the series \code{x}.
#' 
#' @param ctr.outlier Logical. Should we pass detected outliers to meanshift and
#' jump detection algorithms
#' 
#' @param var.diff Logical. Should the series be differenced before algorithm is
#' run.
#' 
#' @param type Character. Select which anomaly test to run. Defaults to all.
#' 
#' @param STL Logical. Should the series be seasonally detrended before isolation
#' forest is run.
#' 
#' @param sigma Numeric. The significance level for outlier detection.
#' 
#' @param ccF Numeric. The significance level for outlier detection. See Details.
#' 
#' @param width Integer. Size of window
#' 
#' @param ccB Numeric. Threshold level for mean-shift detection.
#' 
#' @param p Numeric. Threshold level for jump detection.
#'
#' @param consensus Logical. If \code{TRUE}, report anomalies only when all
#' algorithms agree.
#'  
#' @param json Logical. If \code{TRUE}, yield output in json.
#'
#' @return A list the \code{x}, \code{start_date}, \code{end_date}, and the output
#' from \code{\link{detect_outliers}}, \code{\link{detect_meanshift}}, and
#' \code{\link{detect_varshift}}.
#' 
#' @export
#'
anomaly_detection <- function(x, start_date=NULL, end_date=NULL,
                              type=c("outlier", "meanshift", "varshift"),
                              ctr.outlier=TRUE, var.diff=FALSE,
                              STL=FALSE, sigma=0.01, ccF=2.5,
                              width=7, ccB=2, p=0.01,
                              consensus=TRUE, json=TRUE)
{
  type <- type[type %in% c("outlier", "meanshift", "varshift")]
  outLST <- muLST <- varLST <- list()
  if (any(c("outlier" %in% type, ctr.outlier))) {
    outLST <- detect_outliers(x, STL=STL, sigma=sigma, ccF=ccF)
  }
  idx <- if (ctr.outlier)
    intersect(outLST$idxBLKAVG, outLST$dxSEQAVG)
  if ("meanshift" %in% type) {
    muLST <- detect_meanshift(x, width=width, ccB=ccB, p=p,
                              idx.outliers=idx)
  }
  if ("varshift" %in% type) {
    varLST <- detect_varshift(x, width=width, var.diff=var.diff,
                              idx.outliers=idx)
  }
  out <- c(list(series=x, start_date=start_date, end_date=end_date),
           outLST, muLST, varLST)
  if (consensus)
    out <- merge_anomaly(out)
  if (json)
    return(jsonlite::toJSON(out, auto_unbox=TRUE))
  return(structure(out, class="otdt"))
}



#' Perform Anomaly Detection with Series Segmentation
#'
#' @param x Integer. A vector of integers.
#' 
#' @param start_date Character. Date of first observation in the series \code{x}.
#' 
#' @param end_date Character. Date of last observation in the series \code{x}.
#' 
#' @param type Character. Select which anomaly test to run. Defaults to all.
#'
#' @param ctr.outlier Logical. Should we pass detected outliers to meanshift and
#' jump detection algorithms
#' 
#' @param alpha Numeric. Significance level for segmentation process
#' 
#' @param var.diff Logical. Should the series be differenced before algorithm is
#' run.
#' 
#' @param STL Logical. Should the series be seasonally detrended before isolation
#' forest is run.
#' 
#' @param sigma Numeric. The significance level for outlier detection.
#' 
#' @param ccF Numeric. The significance level for outlier detection. See Details.
#' 
#' @param width Integer. Size of window
#' 
#' @param ccB Numeric. Threshold level for mean-shift detection.
#' 
#' @param p Numeric. Threshold level for jump detection.
#'
#' @param consensus Logical. If \code{TRUE}, report anomalies only when all
#' algorithms agree.
#'
#' @param json Logical. If \code{TRUE}, yield output in json.
#'  
#' @return A list the \code{x}, \code{start_date}, \code{end_date}, and the output
#' from \code{\link{detect_outliers}}, \code{\link{detect_meanshift}}, and
#' \code{\link{detect_varshift}}.
#' 
#' @export
#'
segment_anomaly <- function(x, start_date=NULL, end_date=NULL,
                            type=c("outlier", "meanshift", "varshift"),
                            ctr.outlier=TRUE, alpha=0.01, var.diff=FALSE,
                            STL=FALSE, sigma=0.01, ccF=2.5,
                            width=7, ccB=2, p=0.01,
                            consensus=TRUE, json=TRUE)
{
  type <- type[type %in% c("outlier", "meanshift", "varshift")]
  stopifnot(!is.null(start_date), !is.null(end_date))
  dates <- seq(as.Date(start_date), as.Date(end_date), by="day")
  re_index <- function(lst, sq)
  {
    lst <- lst[grepl("^idx", names(lst))]
    Map(function(x, y) x + min(y) - 1L, lst, list(sq))
  }
  idx.seg <- segment(x, alpha=alpha) + 1L # C++ source code indexes start at 0
  if (any(c("outlier" %in% type, ctr.outlier))) {
    padidx <- regroup_segments(segment_index(idx.seg, overlap=0), mindays=20)
    lstx <- lapply(padidx, function(i) x[i])
    outLST1 <- Map(detect_outliers, x=lstx, STL=STL, sigma=sigma, ccF=ccF)
    outLST2 <- Map(re_index, outLST1, padidx)
    idx.out <- list(idxDIST=Reduce(union, lapply(outLST2, `[[`, "idxDIST")),
                    idxISOF=Reduce(union, lapply(outLST2, `[[`, "idxISOF")))
  }
  if ("meanshift" %in% type) {
    idx.ctr <- if (ctr.outlier) {
      Reduce(intersect, idx.out)
    }
    idx.mu <- detect_meanshift(x, width=width, ccB=ccB, p=p, idx.outliers=idx.ctr)
  }
  if ("varshift" %in% type) {
    idx.ctr <- if (ctr.outlier) {
      unlist(lapply(outLST2, function(lst) Reduce(intersect, lst)))
    }
    varshft <- detect_varshift(x, width=width, var.diff=var.diff, idx.outliers=idx.ctr)
  }
  out <- c(list(series=x, start_date=start_date, end_date=end_date),
           idx.out, idx.mu, varshft, segments=list(idx.seg))
  if (consensus)
    out <- merge_anomaly(out)
  if (json)
    return(jsonlite::toJSON(out, auto_unbox=TRUE))
  return(structure(out, class="otdt"))
}


#' Plot detected anomalies
#' 
#' @param obj S3 object created by call to \code{\link{anomaly_detection}} or
#' \code{\link{segment_anomaly}}.
#'
#' @param ... other paramaters to pass to \code{plot}
#'
#' @importFrom graphics abline par plot
#'
#' @importFrom grDevices rgb
#'
#' @export
#'  
oplot <- function(obj, ...)
{
  x <- obj$series
  idxISOF <- obj$idxISOF
  idxDIST <- obj$idxDIST
  idxBLKAVG <- obj$idxBLKAVG
  idxSEQAVG <- obj$idxSEQAVG
  idx <- intersect(idxBLKAVG, idxSEQAVG)
  idxVAR <- obj$idxVAR
  dates <- if (!is.null(obj$start_date) && !is.null(obj$end_date)) {
    seq(as.Date(obj$start_date), as.Date(obj$end_date), by="day")
  } else {
    seq_len(length(x))
  }
  on.exit(par(mfrow=c(1,1)))
  par(mfrow=c(3,1))
  col1 <- rgb(255, 99, 71, 150, NULL, 255)  # tomato
  col2 <- rgb(174, 54, 204, 30, NULL, 255)  # purple
  col3 <- rgb(57, 89, 71, 150, NULL, 255)   # green
  col4 <- rgb(0, 52, 127, 150, NULL, 255)   # blue
  col5 <- rgb(242, 159, 5, 30, NULL, 255)   # orange
  col6 <- rgb(140, 63, 91, 150, NULL, 255)  # fusia
  plot2 <- function(...)
  {
    plot(dates, x, type="l", ...)
    if ("segments" %in% names(obj)) {
      abline(v=dates[obj$segments], lwd=2, lty=2)
    }
  }
  # Graph 1
  plot2()
  if (length(idxISOF) > 0) {
    abline(v=dates[setdiff(idxISOF, idxDIST)], lwd=3, col=col1)
    abline(v=dates[setdiff(idxDIST, idxISOF)], lwd=3, col=col3)
    abline(v=dates[intersect(idxDIST, idxISOF)], lwd=3, col=col4)
  }
  plot2()
  if (length(idxBLKAVG) + length(idxSEQAVG) > 0) {
    abline(v=dates[setdiff(idxBLKAVG, idx)], lwd=3, col=col5)
    abline(v=dates[setdiff(idxSEQAVG, idx)], lwd=3, col=col2)
    abline(v=dates[idx], lwd=3, col=col3)
  }
  # Graph 3
  plot2()
  if (length(idxVAR) > 0)
    abline(v=dates[idxVAR], lwd=3, col=col4)
}


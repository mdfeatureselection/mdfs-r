#' Compute p-values from information gains and return MDFS
#'
#' @param IG max conditional information gains
#' @param dimensions number of dimensions
#' @param divisions number of divisions
#' @param discretizations number of discretizations (matters only for dimensions=discretizations=1)
#' @param contrast.mask boolean mask on \code{IG} specifying which variables are contrast variables (or \code{NULL} if none, otherwise at least 3 variables must be marked)
#' @param response.divisions number of response divisions (i.e. categories-1)
#' @param df vector of degrees of freedom for each variable (optional)
#' @param ig.in.bits \code{TRUE} if input is in binary log (as opposed to natural log)
#' @param ig.doubled \code{TRUE} if input is doubled (to follow the chi-squared distribution)
#' @param param.median \code{TRUE} if FILL ME WITH TEXT
#' @param n.outliers FILL ME WITH TEXT
#' @param max.outliers FILL ME WITH TEXT
#' @param level FILL ME WITH TEXT
#' @return A \code{\link{data.frame}} with class set to \code{MDFS}. Can be coerced back to \code{data.frame} using \code{\link{as.data.frame}}.
#'
#'  The following columns are present:
#'  \itemize{
#'    \item \code{IG} -- information gains (input copy)
#'    \item \code{chi.squared.p.value} -- chi-squared p-values
#'    \item \code{p.value} -- theoretical p-values
#'  }
#'
#'  Additionally the following \code{\link{attributes}} are set:
#'  \itemize{
#'   \item \code{run.params} -- run parameters
#'   \item \code{dist.param} -- distribution parameter
#'   \item \code{fdr.outliers} -- FILL ME WITH TEXT
#'   \item \code{n.outliers} -- FILL ME WITH TEXT
#'  }
#' @examples
#' ComputePValue(madelon$IG.2D, dimensions = 1, divisions = 1, discretizations=1)
#' @importFrom stats pchisq ks.test
#' @export
ComputePValue <- function(
    IG,
    dimensions,
    divisions,
    discretizations, #this matters only for dimensions=discretizations=1
    contrast.mask=NULL,
    response.divisions = 1,
    df = NULL,
    ig.in.bits = TRUE,
    ig.doubled = FALSE,
    param.median = TRUE,
    n.outliers = NULL,
    max.outliers = NULL,
    level = 0.05
) {
  #check the reasonability of input
  
  if (length(IG)<4 || sum(is.na(as.numeric(IG)))>0 || sum(!is.numeric(IG))>0) {
    stop('IG has to be a numeric vector of length above 3, without N/A values')
  }
  
  if (as.integer(dimensions) != dimensions || dimensions < 1) {
    stop('Dimensions has to be a positive integer')
  }
  
  if (as.integer(divisions) != divisions || divisions < 1) {
    stop('Divisions has to be a positive integer')
  }
  
  if (as.integer(discretizations) != discretizations || discretizations < 1) {
    stop('Discretizations has to be a positive integer')
  }
  
  if (as.integer(response.divisions) != response.divisions || response.divisions < 1) {
    stop('Response.divisions has to be a positive integer')
  }
  
  if (!is.null(df)) {
    if (sum(is.na(as.integer(df)))>0 || sum(as.integer(df)!=df)>0 || sum(df<1)>0) stop('Df has to be a vector of positive integers')
    if (length(df) != 1 && length(df) != length(IG)) stop('Df has to have the same length as IG or 1')
  }
  
  if (dimensions>1 || discretizations>1) {
    if (length(contrast.mask)==length(IG)) contrast.mask<-as.logical(contrast.mask)
    if (is.null(contrast.mask) || length(IG[contrast.mask])<4) {
      stop('Contrast.mask has to specify a vector of length above 3')
    }
    
    if (as.numeric(level) != level || level<=0 || level>=1) {
      stop('Level has to be a number between 0 and 1')
    } 
    
    if (!is.null(n.outliers)) {
      if (as.integer(n.outliers) != n.outliers || 
          n.outliers<0 || 
          n.outliers>length(IG[contrast.mask])-4) {
        stop('N.outliers must be a positive integer smaller than the number of contrast variables')
      }
    } else 
      if (is.null(max.outliers)) {
        max.outliers <- ceiling(length(IG[contrast.mask])/4)
      } else {
        if (as.integer(max.outliers) != max.outliers || 
            max.outliers<0 || 
            max.outliers>length(IG[contrast.mask])-4) {
          stop('Max.outliers must be a positive integer smaller than the number of contrast variables')
        }
      }
  } 
  
  IG.original <- IG
  
  #bits to nats and the factor 2
  if (ig.in.bits) IG <- log(2)*IG
  if (!ig.doubled) IG <- 2*IG
  
  #IG must be positive
  IG[IG<1e-20] <- 1e-20
  
  #degrees of freedom - if not specified by the user
  if (is.null(df)) df<-response.divisions*divisions*(divisions+1)^(dimensions-1)
  
  #chi-squared p-value
  chisq <- pchisq(IG,df,lower.tail=FALSE)
  
  #for dimensions=discretizations=1, p-value=chisq
  if (dimensions==1 && discretizations==1) {
    gamma <- 1
    n.important <- length(IG)
    fdr.outliers <- NULL
    n.contrast<- 0
    p.value <- chisq
  } else {
    #otherwise use contrast variables to estimate the distribution parameter 
    chisq.log <- -pchisq(IG,df,log.p=TRUE)
    
    chisq.log.contrast <- chisq.log[contrast.mask]
    chisq.log.contrast <- chisq.log.contrast[order(chisq.log.contrast)]
    n.contrast <- length(chisq.log.contrast)
    min.important<-n.contrast-max.outliers
    
    #compute cumulative means of chisq.contrast
    #  means <- cumsum(chisq.log.contrast)/(1+(1:n.contrast))
    #  medians are more robust to outliers  
    medians <- (c(0,chisq.log.contrast)[1+floor((1:n.contrast)/2)]+chisq.log.contrast[ceiling((1:n.contrast)/2)])/2
    
    #compute fdr that a least IG is an outlier.
    # fdr.outliers <- p.adjust(exp(-chisq.log.contrast/means),'fdr') 
    # p.adjust written explicitely  
    # fdr.outliers <- cummin((exp(-chisq.log.contrast/means)*n.contrast/(n.contrast:1))[(min.important+1):n.contrast])
    # medians are more robust, lack of the factor log(2) gives a linear growth
    fdr.outliers <- cummin((exp(-chisq.log.contrast/medians)*n.contrast/(n.contrast:1))[(min.important+1):n.contrast])
    
    #number of contrast scores not being outliers
    if (!is.null(n.outliers)) {
      n.important <- n.contrast-n.outliers
    } else {
      n.important<-min.important+max(0, which(fdr.outliers>level))
      
      if (n.important==min.important) {
        warning("Border value reached for number of outliers")
      } 
    }  
    
    #distribution parameter and eventual p-value
    if (param.median) {
      M<-median(chisq.log.contrast[1:n.important])
      gamma <- log(2)/M
      
      #p.value <- -expm1(-gamma*chisq.log)
      #this is the exact distribution
      p.value <- ifelse(chisq.log/M>=n.important/2,1,-expm1(-(lgamma(n.important/2-chisq.log/M)+lgamma(n.important)-lgamma(n.important-chisq.log/M)-lgamma(n.important/2))))
      
    } else {
      S<-sum(chisq.log.contrast[1:n.important])
      gamma <- (n.important+1)/S
      
      #p.value <- -expm1(-gamma*chisq.log)
      #Lomax distribution is more precise
      p.value <- -expm1(-(n.important+1)*log1p(chisq.log/S))
      
    }
  }
  
  #the final result
  result <- data.frame(
    IG = IG.original,
    chi.squared.p.value = chisq,
    p.value = p.value
  )
  class(result) <- 'MDFS'
  
  attr(result, 'run.params') <- list(
    contrast.mask           = contrast.mask,
    dimensions              = dimensions,
    divisions               = divisions,
    discretizations         = discretizations,
    response.divisions      = response.divisions,
    df                      = df,
    ig.in.bits              = ig.in.bits,
    ig.doubled              = ig.doubled,
    n.outliers              = n.outliers,
    max.outliers            = max.outliers,
    level                   = level)
  
  attr(result, 'dist.param') <- gamma
  attr(result, 'fdr.outliers') <- fdr.outliers
  attr(result, 'n.outliers') <- n.contrast-n.important
  
  return(result)
}


#' Find indices of relevant variables
#'
#' @param fs feature selector
#' @param ... arguments passed to methods
#' @return indices of important variables
#' @export
RelevantVariables <- function(fs, ...) {
  UseMethod('RelevantVariables')
}

#' Find indices of relevant variables from MDFS
#'
#' @details
#' In case of FDR control it is recommended to use Benjamini-Hochberg-Yekutieli p-value adjustment
#' method (\code{"BY"} in \code{\link[stats]{p.adjust}}) due to unknown dependencies between tests.
#'
#' @param fs an MDFS object
#' @param level statistical significance level
#' @param p.adjust.method method as accepted by \code{\link[stats]{p.adjust}} (\code{"BY"} is recommended for FDR, see Details)
#' @param ... ignored
#' @return indices of relevant variables
#' @importFrom stats p.adjust
#' @export
RelevantVariables.MDFS <- function(fs, level=0.05, p.adjust.method="holm", ...) {
  contrast.mask <- attr(fs, 'run.params')$contrast.mask
  p.value <- if(is.null(contrast.mask)) { fs$p.value } else { fs$p.value[!contrast.mask] }
  return(which(p.adjust(p.value, method=p.adjust.method)<level))
}

#' Plot MDFS details
#'
#' @param x an MDFS object
#' @param plots plots to plot (ig for max IG, c for chi-squared p-values, p for p-values)
#' @param ... passed on to \code{\link[graphics]{plot}}
#' @importFrom graphics plot
#' @export
plot.MDFS <- function(x, plots=c('ig', 'c', 'p'), ...) {
  ord <- order(x$IG,decreasing=TRUE)
  for (plt in plots) {
    switch(plt,
           ig = plot(x$IG[ord], xlab='index', ylab=expression('I'[max]*'(X)'), ...),
           c = plot(x$chi.squared.p.value[ord], seq(ord)/length(ord), xlab='chi-squared p-value', ylab='experimental p-value', ...),
           p = plot(x$p.value[ord], seq(ord)/length(ord), xlab='theoretical p-value', ylab='experimental p-value', ...),
           stop(paste('I don\'t know how to plot', plt)))
  }
}

#' as.data.frame S3 method implementation for MDFS
#'
#' @param x an MDFS object
#' @param ... ignored
#' @return data.frame
#' @export
as.data.frame.MDFS <- function(x, ...) {
  class(x) <- 'data.frame'
  return(x)
}
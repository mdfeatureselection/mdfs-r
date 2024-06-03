#' Run end-to-end MDFS for discretized data with same number of categories
#'
#' @details
#' In case of FDR control it is recommended to use Benjamini-Hochberg-Yekutieli p-value adjustment
#' method (\code{"BY"} in \code{\link[stats]{p.adjust}}) due to unknown dependencies between tests.
#'
#' @param data input data.frame where columns are variables and rows are observations (all numeric)
#' @param decision decision variable as a boolean vector of length equal to number of observations
#' @param n.contrast number of constrast variables (defaults to max of 1/10 of variables number and 30)
#' @param dimensions number of dimensions 
#' @param pc.xi parameter xi used to compute pseudocounts (the default is recommended not to be changed)
#' @param p.adjust.method method as accepted by \code{\link[stats]{p.adjust}} (\code{"BY"} is recommended for FDR, see Details)
#' @param level statistical significance level
#' @param seed seed for PRNG used during discretizations (\code{NULL} for random)
#' @return A \code{\link{list}} with the following fields:
#'  \itemize{
#'    \item \code{contrast.indices} -- indices of variables chosen to build contrast variables
#'    \item \code{statistic} -- vector of statistic's values (IGs) for corresponding variables
#'    \item \code{p.value} -- vector of p-values for corresponding variables
#'    \item \code{adjusted.p.value} -- vector of adjusted p-values for corresponding variables
#'    \item \code{relevant.variables} -- vector of relevant variables indices
#'  }
#' @examples
#' \donttest{
#' as.data.frame(lapply(as.data.frame(madelon$data), function(COL) 1.*(COL > median(COL)) ) )-> discr_data
#' MDFS.discrete(discr_data, madelon$decision, dimensions = 2)
#' }
#' @importFrom stats p.adjust
#' @export
MDFS.discrete<-function(data, 
                        decision, 
                        n.contrast = max(ncol(data), 30),
                        dimensions = 1,
                        pc.xi = 0.25,
                        p.adjust.method = "holm",
                        level = 0.05,
                        seed = NULL
) {                          
  data <- as.data.frame(data)
  if (!is.null(seed)) set.seed(seed)
  
  if (dimensions > 1  && (is.null(n.contrast) || n.contrast < 30)) stop('Specify the number of contrast variables!')
  
  #contrast variables
  if (n.contrast > 0) {
    contrast <- GenContrastVariables(data, n.contrast)
    contrast.indices <- contrast$indices
    contrast_data <- contrast$contrast_data
    contrast.mask <- c(rep.int(F, ncol(data)), rep.int(T, ncol(contrast_data)))
  }
  else {
    contrast.mask <- contrast.indices <- contrast_data <- NULL
  }
  
  #compute IG
  IG <- rep(0, ncol(data)+ncol(contrast_data))
  
  MIG.Result <- ComputeMaxInfoGainsDiscrete(data=data,
                                            decision=decision,
                                            contrast_data=contrast_data,
                                            dimensions = dimensions,
                                            pc.xi = pc.xi)
  
  IG <- IG+c(MIG.Result$IG, attr(MIG.Result, "contrast_igs"))
  
  #p-value
  divisions <- length(unique(data[,1]))-1
  PV <- ComputePValue(IG, dimensions, divisions, discretizations=1, response.divisions=1, contrast.mask=contrast.mask)
  
  #result
  result <- list(contrast.indices = contrast.indices,
                 statistic = PV$statistic,
                 p.value = PV$p.value,
                 adjusted.p.value = p.adjust(PV$p.value,method=p.adjust.method),
                 relevant.variables = which(p.adjust(PV$p.value,method=p.adjust.method)<level))
  
  return(result)
}
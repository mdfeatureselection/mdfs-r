#' Symmetrize matrix by choosing between one of lower/upper half elements
#'
#' @param A square matrix to be symmetrized
#' @param pairwiseSelectionFunction vectorized function for producing element of a symmetrized matrix based on lower/upper half elements 
#' @return A symmetric matrix
#' @examples
#' symmetrizeMatrix(as.matrix(madelon$data[1:30,1:30]))
#' @export
symmetrizeMatrix <- function (A, pairwiseSelectionFunction=pmin) {
  ltr_selector <- lower.tri(A)
  A[ltr_selector] -> lowerHalf
  A <- t(A)
  A[ltr_selector] -> upperHalf
  A <- t(A)
  replacement <- pairwiseSelectionFunction(lowerHalf, upperHalf)
  A[ltr_selector] <- replacement
  A <- t(A)
  A[ltr_selector] <- replacement
  A <- t(A)
  return(A)
}



#' Compute similarity matrix of discrete variables with same number of categories based on their mutual information
#'
#' @details
#' In case of FDR control it is recommended to use Benjamini-Hochberg-Yekutieli p-value adjustment
#' method (\code{"BY"} in \code{\link[stats]{p.adjust}}) due to unknown dependencies between tests.
#
#' @param data input data.frame where columns are variables and rows are observations (all numeric, discrete variables with same number of categories)
#' @param divisions number of classes of each variable (same for all variables in data)
#' @param p.adjust.method method as accepted by \code{\link[stats]{p.adjust}} (\code{"BY"} is recommended for FDR, see Details)
#' @param level statistical significance level
#' @return A \code{\link{list}} of symmetric matrices, with the following fields:
#'  \itemize{
#'    \item \code{S} -- similarity matrix of variables where \code{S[i,j]} is equal to \code{MI[i,j]} if \code{pv.adjusted[i,j]} is \code{<level} else 0
#'    \item \code{pv.adjusted} -- \code{pv.adjusted[i,j]} is adjusted p-value of mutual information \code{MI[i,j]}. 
#'    \item \code{pv} -- \code{pv[i,j]} is p-value of mutual information \code{MI[i,j]} 
#'    \item \code{MI} -- \code{MI[i,j]} is mutual information computed between variables \code{data[,i]} and \code{data[,j]}
#'  }
#' @examples
#' \donttest{
#' as.data.frame(lapply(as.data.frame(madelon$data), function(COL) 1.*(COL > median(COL)) ) )-> discr_data
#' pairwiseMIsimilarityDiscrete(discr_data, 1)
#' }
#' @importFrom stats p.adjust
#' @export
pairwiseMIsimilarityDiscrete			<- function(data, 
                                           divisions=1, 
                                           p.adjust.method="holm", 
                                           level=0.05
){
  data <- as.data.frame(data)
  
  
  
  #compute IG
  MI <- matrix(0, ncol(data), ncol(data))
  MI <- ComputeInterestingTuplesDiscrete(data=data,
                                         dimensions=2,
                                         return.matrix=TRUE,
                                         stat_mode="MI")*nrow(data)
  diag(MI) <- 0
  
  divisions <- length(unique(data[,1]))-1
  df <- divisions^2
  
  #p-value
  pv<-MI
  pv[,]<- ComputePValue(IG=as.vector(MI), 
                        dimensions=1, 
                        divisions=divisions, 
                        response.divisions=divisions,
                        discretizations=1,
                        df=df)$p.value
  diag(pv) <- 1
  
  #adjust p-value
  ltr_selector <- lower.tri(pv)
  adjustedFlatPV <- p.adjust(pv[ltr_selector], method=p.adjust.method)
  
  pv.adj <- matrix(nrow=ncol(data), ncol=ncol(data))
  pv.adj[ltr_selector] <- adjustedFlatPV
  pv.adj <- t(pv.adj)
  pv.adj[ltr_selector] <- adjustedFlatPV
  pv.adj <- t(pv.adj)
  diag(pv.adj) <- 1
  
  #similarity
  S <- MI
  S[pv.adj >= level] <-0
  
  return (list(S=S,
               pv.adjusted=pv.adj,
               pv=pv,
               MI=MI))
}


#' Compute similarity matrix of continuous data based on their mutual information after discretization
#'
#' @details
#' In case of FDR control it is recommended to use Benjamini-Hochberg-Yekutieli p-value adjustment
#' method (\code{"BY"} in \code{\link[stats]{p.adjust}}) due to unknown dependencies between tests.
#'
#' @param data input data.frame where columns are variables and rows are observations (all numeric, discrete variables with same number of categories)
#' @param divisions number of divisions
#' @param discretizations number of discretizations
#' @param n.contrast number of contrast variables
#' @param p.adjust.method method as accepted by \code{\link[stats]{p.adjust}} (\code{"BY"} is recommended for FDR, see Details)
#' @param level statistical significance level
#' @param seed random seed for randomized discretization procedure
#' @return A \code{\link{list}} of symmetric matrices, with the following fields:
#'  \itemize{
#'    \item \code{S} -- similarity matrix of variables where \code{S[i,j]} is equal to \code{MI[i,j]} if \code{pv.adjusted[i,j]} is \code{<level} else 0
#'    \item \code{pv.adjusted} -- \code{pv.adjusted[i,j]} is adjusted p-value of mutual information \code{MI[i,j]}. 
#'    \item \code{pv} -- \code{pv[i,j]} is p-value of mutual information \code{MI[i,j]} 
#'    \item \code{MI} -- \code{MI[i,j]} is mutual information computed between variables \code{data[,i]} and \code{data[,j]}
#'  }
#' @examples
#' \donttest{
#' pairwiseMIsimilarity(madelon$data, 1)
#' }
#' @importFrom stats p.adjust
#' @export
pairwiseMIsimilarity <- function(data, 
                                 divisions=1,
                                 discretizations=1,
                                 n.contrast = max(ncol(data), 30),
                                 p.adjust.method="holm", 
                                 level=0.05,
                                 seed=NULL){
  
  if (discretizations > 1 && (is.null(n.contrast) || n.contrast < 30)) stop('Specify the number of contrast variables!')
  
  data <- as.data.frame(data)
  if (!is.null(seed)) set.seed(seed)
  
  
  #transform many confounders into a single one
  
  
  #contrast variables
  if (n.contrast > 0) {
    contrast <- GenContrastVariables(data, n.contrast)
    contrast.indices <- contrast$indices
    contrast_data <- contrast$contrast_data
  }
  else {
    contrast.mask <- contrast.indices <- contrast_data <- NULL
  }
  
  #compute IG
  MI <-  matrix(0, ncol(data), ncol(data))
  if (n.contrast>0) MI.contrast <- matrix(0, n.contrast, n.contrast) else MI.contrast<-as.matrix(0)
  MI <- ComputeInterestingTuples(data=data,
                                 dimensions=2,
                                 divisions=divisions,
                                 discretizations=discretizations,
                                 return.matrix=TRUE,
                                 stat_mode="MI")*nrow(data)
  
  if (n.contrast > 0)   MI.contrast <- MI.contrast+ComputeInterestingTuples(data=contrast_data,
                                                                            dimensions=2,
                                                                            divisions=divisions,
                                                                            discretizations=discretizations,
                                                                            return.matrix=TRUE,
                                                                            stat_mode="MI")*nrow(data)
  
  
  diag(MI) <- diag(MI.contrast) <- 0
  
  ltr_selector <- lower.tri(MI)
  ltr_selector.contrast <- lower.tri(MI.contrast) 
  
  flatMI <- c(MI[ltr_selector],MI.contrast[ltr_selector.contrast])
  contrast.mask <- c(rep.int(F, sum(ltr_selector)), rep.int(T, sum(ltr_selector.contrast)))
  
  
  #p-value
  flatPV <- ComputePValue(IG=flatMI, 
                          dimensions=1, 
                          divisions=divisions, 
                          response.divisions=divisions,
                          discretizations=discretizations,
                          contrast.mask=contrast.mask)$p.value	    
  
  adjustedFlatPV <- p.adjust(flatPV[!contrast.mask], method=p.adjust.method)
  pv <- pv.adj <- matrix(nrow=ncol(data), ncol=ncol(data))
  pv[ltr_selector] <- flatPV[!contrast.mask]
  pv <- t(pv)
  pv[ltr_selector] <- flatPV[!contrast.mask]
  pv <- t(pv)
  
  pv.adj[ltr_selector] <- adjustedFlatPV
  pv.adj <- t(pv.adj)
  pv.adj[ltr_selector] <- adjustedFlatPV
  pv.adj <- t(pv.adj)
  
  diag(pv) <- diag(pv.adj) <- 1
  
  #similarity
  S <- MI
  S[pv.adj >= level] <-0
  
  return (list(S=S,
               pv.adjusted=pv.adj,
               pv=pv,
               MI=MI))
}



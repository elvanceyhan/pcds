#AuxDomination.R;
#Contains the ancillary functions for finding the domination number of a given digraph
#################################################################

#' @title Approximate domination number and approximate dominating set by the greedy algorithm
#'
#' @description Returns the (approximate) domination number and the indices (i.e., row numbers) for the corresponding
#' (approximate) minimum dominating set based on the incidence matrix \code{Inc.Mat} of a graph or a digraph
#' by using the greedy algorithm (\insertCite{chvatal:1979;textual}{pcds}).
#' Here the row number in the incidence matrix corresponds to the index of the vertex (i.e., index of the
#' data point). The function works whether loops are allowed or not (i.e., whether the first diagonal
#' is all 1 or all 0). This function may yield the actual domination number or overestimates it.
#'
#' @param Inc.Mat A square matrix consisting of 0's and 1's which represents the incidence matrix of
#' a graph or digraph.
#'
#' @return A \code{list} with two elements
#' \item{dom.num}{The cardinality of the (approximate) minimum dominating set found by the greedy algorithm.
#' i.e., (approximate) domination number of the graph or digraph whose incidence matrix \code{Inc.Mat} is given
#' as input.}
#' \item{ind.domset}{Indices of the rows in the incidence matrix \code{Inc.Mat} for the ((approximate) minimum
#' dominating set). The row numbers in the incidence matrix correspond to the indices of the vertices
#' (i.e., indices of the data points).}
#'
#' @seealso \code{\link{dom.exact}}, \code{\link{PEdom1D}}, \code{\link{PEdomTri}}, \code{\link{PEdom.nd}},
#' and \code{\link{IndNCSdomUBtri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' n<-5
#' M<-matrix(sample(c(0,1),n^2,replace=TRUE),nrow=n)
#' diag(M)<-1
#'
#' dom.greedy(M)
#'
#' @export dom.greedy
dom.greedy <- function(Inc.Mat)
{
  inc.mat<-as.matrix(Inc.Mat)
  nr<-nrow(inc.mat); nc<-ncol(inc.mat)

  if (nr!=nc)
  {stop('Inc.Mat must a square matrix')}

  if (sum(inc.mat!=0 & inc.mat!=1)!=0 )
  {stop('Inc.Mat must have entries equal to 0 or 1')}

  if (!(all(diag(inc.mat)==0) || all(diag(inc.mat)==1) ))
  {stop('Inc.Mat must have all zeroes (when loops are not allowed) or all ones (when loops are allowed) in the diagonal.')}

  cnt<-nr
  dom.ind<-c()
  vert.ind<-1:nr
  while (cnt>0)
  {
    rsums<-apply(inc.mat,1,sum)
    max.ind<-which(rsums==max(rsums))
    sel.ind<-ifelse(length(max.ind)>1,sample(max.ind,1),max.ind)
    dom.ind<-c(dom.ind,vert.ind[sel.ind])
    nghd<-which(inc.mat[sel.ind,]==1)  #indices in the nghd of the max dominating vertex
    vert.ind<-vert.ind[-nghd]
    inc.mat1<-inc.mat[-nghd,]
    if (length(inc.mat1)>0)
    {ifelse(ncol(as.matrix(inc.mat1))==1,inc.mat2<-matrix(inc.mat1,ncol=length(inc.mat1)),inc.mat2<-as.matrix(inc.mat1))
      inc.mat<-as.matrix(inc.mat2[,-nghd])
      cnt<-nrow(inc.mat)
    } else
    {cnt<-0}
  }
  list(
    approx.dom.num=length(dom.ind), #greedy domination number
    ind.approx.mds=dom.ind #indices of vertices or data points for the greedy domination number
  )
} #end of the function
#'

#################################################################

#' @title Indicator for an upper bound for the domination number by the exact algorithm
#'
#' @description Returns 1 if the domination number is less than or equal to the prespecified value \code{k} and also the indices
#' (i.e., row numbers) of a dominating set of size \code{k} based on the incidence matrix \code{Inc.Mat} of a graph or
#' a digraph. Here the row number in the incidence matrix corresponds to the index of the vertex
#' (i.e., index of the data point). The function works whether loops are allowed or not
#' (i.e., whether the first diagonal is all 1 or all 0).
#' It takes a rather long time for large number of vertices (i.e., large number of row numbers).
#'
#' @param Inc.Mat A square matrix consisting of 0's and 1's which represents the incidence matrix of
#' a graph or digraph.
#' @param k A positive integer for the upper bound (to be checked) for the domination number.
#'
#' @return A \code{list} with two elements
#' \item{domUB}{The upper bound (to be checked) for the domination number. It is prespecified as \code{k}
#' in the function arguments.}
#' \item{IndUBdom}{The indicator for the upper bound for domination number of the graph or digraph being
#' the specified value \code{k} or not. It returns 1 if the upper bound is \code{k}, and 0 otherwise based on the incidence
#' matrix \code{Inc.Mat} of the graph or digraph.}
#' \item{ind.domset}{Indices of the rows in the incidence matrix \code{Inc.Mat} that correspond to the vertices in the
#' dominating set of size \code{k} if it exists, otherwise it yields \code{NULL}.}
#'
#' @seealso \code{\link{dom.exact}} and \code{\link{dom.greedy}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' n<-10
#' M<-matrix(sample(c(0,1),n^2,replace=TRUE),nrow=n)
#' diag(M)<-1
#'
#' dom.greedy(M)
#' IndUBdom(M,2)
#'
#' for (k in 1:n)
#' print(c(k,IndUBdom(M,k)))
#' }
#'
#' @export IndUBdom
IndUBdom <- function(Inc.Mat,k)
{
  inc.mat<-as.matrix(Inc.Mat)
  nr<-nrow(inc.mat); nc<-ncol(inc.mat)

  if (nr!=nc)
  {stop('Inc.Mat must a square matrix')}

  if (sum(inc.mat!=0 & inc.mat!=1)!=0 )
  {stop('Inc.Mat must have entries equal to 0 or 1')}

  if (!(all(diag(inc.mat)==0) || all(diag(inc.mat)==1) ))
  {stop('Inc.Mat must have all zeroes (when loops are not allowed) or all ones (when loops are allowed) in the diagonal.')}

  if (all(diag(inc.mat)==0))
  {diag(inc.mat)<-1}

  xc<-combinat::combn(1:nr,k); N1<-choose(nr,k);
  xc<-matrix(xc,ncol=N1)  #required when N1=1
  dom<-0; j<-1; ind.domset<-c();
  while (j<=N1 && dom==0)
  { ifelse(nrow(xc)==1,dom.check<-inc.mat[xc[,j],],dom.check<-apply(inc.mat[xc[,j],],2,sum))

    if (all(dom.check>=1))  #this is where domination is checked
    {dom<-1; ind.domset<-xc[,j]}
    j<-j+1;
  }

  list(domUB=k, #upper bound for the domination number
       IndUBdom=dom, #indicator that domination number <=k
       ind.domset=ind.domset #indices of the vertices in a dominating set of size k (if exists)
  )
} #end of the function
#'

#################################################################

#' @title Exact domination number (i.e., domination number by the exact algorithm)
#'
#' @description Returns the (exact) domination number based on the incidence matrix \code{Inc.Mat} of a graph or a digraph
#' and the indices (i.e., row numbers of \code{Inc.Mat}) for the corresponding (exact) minimum dominating set.
#' Here the row number in the incidence matrix corresponds to the index of the vertex (i.e., index of the data
#' point). The function works whether loops are allowed or not (i.e., whether the first diagonal is all 1 or
#' all 0). It takes a rather long time for large number of vertices (i.e., large number of row numbers).
#'
#' @inheritParams dom.greedy
#'
#' @return A \code{list} with two elements
#' \item{dom.num}{The cardinality of the (exact) minimum dominating set i.e., (exact) domination number of the
#' graph or digraph whose incidence matrix \code{Inc.Mat} is given as input.}
#' \item{ind.mds}{The vector of indices of the rows in the incidence matrix \code{Inc.Mat} for the (exact) minimum
#' dominating set. The row numbers in the incidence matrix correspond to the indices of the vertices
#' (i.e., indices of the data points).}
#'
#' @seealso \code{\link{dom.greedy}}, \code{\link{PEdom1D}}, \code{\link{PEdomTri}}, \code{\link{PEdom.nd}},
#' and \code{\link{IndNCSdomUBtri}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' n<-10
#' M<-matrix(sample(c(0,1),n^2,replace=TRUE),nrow=n)
#' diag(M)<-1
#'
#' dom.greedy(M)
#' IndUBdom(M,2)
#' dom.exact(M)
#' }
#'
#' @export dom.exact
dom.exact <- function(Inc.Mat)
{
  inc.mat<-as.matrix(Inc.Mat)
  nr<-nrow(inc.mat); nc<-ncol(inc.mat)

  if (nr!=nc)
  {stop('Inc.Mat must a square matrix')}

  if (sum(inc.mat!=0 & inc.mat!=1)!=0 )
  {stop('Inc.Mat must have entries equal to 0 or 1')}

  if (!(all(diag(inc.mat)==0) || all(diag(inc.mat)==1) ))
  {stop('Inc.Mat must have all zeroes (when loops are not allowed) or all ones (when loops are allowed) in the diagonal.')}

  if (all(diag(inc.mat)==0))
  {diag(inc.mat)<-1}

  crit<-0; i<-1; ind.mds<-c();
  while (i<=nr && crit==0)
  { DN<-IndUBdom(inc.mat,i)
  if (DN$IndUBdom==1)  #this is where domination is checked
  {dom.num<-i; ind.mds<-DN$ind.domset;
  crit<-1}
  i<-i+1;
  }

  list(dom.num=dom.num, #domination number
       ind.mds=ind.mds #indices of the vertices in a minimum dominating set
  )
} #end of the function
#'

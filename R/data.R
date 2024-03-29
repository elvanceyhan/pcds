#' @title Tree Species in a Swamp Forest
#'
#' @description Locations and species classification of trees in a plot in the Savannah River, SC, USA.
#' Locations are given in meters, rounded to the nearest 0.1 decimal.
#' The data come from a one-hectare (200-by-50m) plot in the Savannah River Site.
#' The 734 mapped stems included 156 Carolina ashes (Fraxinus caroliniana),
#' 215 water tupelos (Nyssa aquatica), 205 swamp tupelos (Nyssa sylvatica), 98 bald cypresses (Taxodium distichum)
#' and 60 stems from 8 additional three species (labeled as Others (OT)).
#' The plots were set up by Bill Good and their spatial patterns described in (\insertCite{good:1982;textual}{pcds}),
#' the plots have been maintained and resampled by Rebecca Sharitz and her colleagues of the Savannah River
#' Ecology Laboratory. The data and some of its description are borrowed from the swamp data entry in the \code{dixon}
#' package in the CRAN repository.
#'
#' See also (\insertCite{good:1982,jones:1994,dixon:NNCTEco2002;textual}{pcds}).
#'
#' @details
#' Text describing the variable (i.e., column) names in the data set.
#' * x,y:  x and y (i.e., Cartesian) coordinates of the trees
#' * live: a categorical variable that indicates the tree is alive (labeled as 1) or dead (labeled as 0)
#' * sp: species label of the trees:
#'   \itemize{
#'   \item{FX: }{Carolina ash (Fraxinus caroliniana)}
#'   \item{NS: }{Swamp tupelo (Nyssa sylvatica)}
#'   \item{NX: }{Water tupelo (Nyssa aquatica)}
#'   \item{TD: }{Bald cypress (Taxodium distichum)}
#'   \item{OT: }{Other species}
#'    }
#' @md
#'
#' @references
#' \insertAllCited{}
#'
#' @source \href{https://pdixon.stat.iastate.edu/datasets/goodplot1.txt}{Prof. Philip Dixon's website}
#'
#' @examples
#' data(swamptrees)
#' plot(swamptrees$x,swamptrees$y, col=as.numeric(swamptrees$sp),pch=19,
#' xlab='',ylab='',main='Swamp Trees')
#'
#' @docType data
#' @keywords datasets
#' @name swamptrees
#' @usage data(swamptrees)
#' @format A data frame with 734 rows and 4 variables
NULL


#' pcds: A package for Proximity Catch Digraphs and Their Applications
#'
#' pcds is a package for generation, computation and visualization of proximity catch digraphs
#' and tests based on them.
#'
#' The pcds package contains the functions for generating patterns of segregation, association, CSR
#' (complete spatial randomness) and Uniform data in one, two and three dimensional cases,
#' for testing these patterns based on two invariants of various families of the proximity catch digraphs (PCDs),
#' (see (\insertCite{ceyhan:Phd-thesis;textual}{pcds}).
#'
#' The graph invariants used in testing spatial point data are the  domination number
#' (\insertCite{ceyhan:dom-num-NPE-Spat2011;textual}{pcds})
#' and arc density (\insertCite{ceyhan:arc-density-PE,ceyhan:arc-density-CS;textual}{pcds})
#' of for two-dimensional data  for visualization of PCDs for one,
#' two and three dimensional data. The PCD families considered are Arc-Slice PCDs,
#' Proportional-Edge PCDs and Central Similarity PCDs.
#'
#' The package also contains visualization tools for these digraphs for 1D-3D vertices. The AS-PCD
#' related tools are provided for 1D and 2D data; PE-PCD related tools are provided for 1D-3D data,
#' and CS-PCD tools are provided for 1D and 2D data.
#'
#' @section The pcds functions:
#' The pcds functions can be grouped as Auxiliary Functions, AS-PCD Functions, PE-PCD Functions,
#' and CS-PCD Functions.
#'
#' @section Auxiliary Functions:
#' Contains the auxiliary functions used in PCD calculations, such as equation of lines for two points,
#' distances between lines and points, generation of points from uniform, segregation and association patterns,
#' checking points inside the triangle etc.
#' In all these functions points are vectors, and data sets are either matrices or data frames.
#'
#' @section Arc-Slice PCD Functions:
#' Contains the functions used in AS-PCD calculations, such as generation of data in a given a triangle and
#' estimation of gamma, arc density, etc.
#'
#' @section Proportional-Edge PCD Functions:
#' Contains the functions used in PE-PCD calculations, such as generation of data in a given interval, triangle and
#' tetrahedron and estimation of gamma, arc density, etc.
#'
#' @section Central-Similarity PCD Functions:
#' Contains the functions used in CS-PCD calculations, such as generation of data in a given interval and triangle and
#' estimation of gamma, arc density, etc.
#'
#' @references
#' \insertAllCited{}
#'
#' @docType package
#' @name pcds
NULL

#' pcds: A package for Proximity Catch Digraphs and Their Applications
#'
#' \code{pcds} is a package for construction and
#' visualization of proximity catch digraphs (PCDs)
#' and computation of two graph invariants of the PCDs and
#' testing spatial patterns using these invariants.
#' The PCD families considered are Arc-Slice (AS) PCDs,
#' Proportional-Edge (PE) PCDs
#' and Central Similarity (CS) PCDs.
#'
#' The graph invariants used in testing spatial point data are
#' the domination number
#' (\insertCite{ceyhan:dom-num-NPE-Spat2011;textual}{pcds})
#' and arc density (\insertCite{ceyhan:arc-density-PE,ceyhan:arc-density-CS;textual}{pcds})
#' of for two-dimensional data.
#'
#' The \code{pcds} package also contains the functions
#' for generating patterns of segregation, association, CSR
#' (complete spatial randomness) and Uniform data
#' in one, two and three dimensional cases,
#' for testing these patterns
#' based on two invariants of various families of the proximity catch digraphs (PCDs),
#' (see (\insertCite{ceyhan:Phd-thesis;textual}{pcds}).
#'
#' Moreover, the package has visualization tools for these digraphs for 1D-3D vertices.
#' The AS-PCD related tools are provided for 1D and 2D data;
#' PE-PCD related tools are provided for 1D-3D data,
#' and CS-PCD tools are provided for 1D and 2D data.
#'
#' @section The \code{pcds} functions:
#' The \code{pcds} functions can be grouped as Auxiliary Functions,
#' AS-PCD Functions,
#' PE-PCD Functions,
#' and CS-PCD Functions.
#'
#' @section Auxiliary Functions:
#' Contains the auxiliary (or utility) functions for constructing and
#' visualizing Delaunay tessellations in 1D and 2D settings,
#' computing the domination number,
#' constructing the geometrical tools,
#'  such as equation of lines for two points,
#' distances between lines and points, checking points inside the triangle etc.,
#' finding the (local) extrema (restricted to Delaunay cells
#' or vertex or edge regions in them).
#'
#' @section Arc-Slice PCD Functions:
#' Contains the functions used in AS-PCD construction,
#' estimation of domination number,
#' arc density, etc in the 2D setting.
#'
#' @section Proportional-Edge PCD Functions:
#' Contains the functions used in PE-PCD construction,
#' estimation of domination number,
#' arc density, etc in the 1D-3D settings.
#'
#' @section Central-Similarity PCD Functions:
#' Contains the functions used in CS-PCD construction,
#' estimation of domination number,
#' arc density, etc in the 1D and 2D setting.
#'
#' @section Point Generation Functions:
#' Contains functions for generation of points from uniform (or CSR),
#' segregation and association patterns.
#'
#' In all these functions points are vectors,
#' and data sets are either matrices or data frames.
#'
#' @references
#' \insertAllCited{}
#'
#' @docType package
#' @name pcds
NULL

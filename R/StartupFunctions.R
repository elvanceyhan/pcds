## Startup functions ------------------------------------

#' .onAttach start message
#'
#' @param libname defunct
#' @param pkgname defunct
#'
#' @return invisible()
.onAttach <- function(libname, pkgname) {
  start_message <- c( "     pcds: functions for proximity catch digraphs,  \n\n",
                      "their visualisation and application in spatial data analysis\n\n",
                      "      by Dr. Elvan Ceyhan <elvanceyhan@gmail.com>\n\n"
  )
  packageStartupMessage(start_message)
  invisible()
}
#'

################################################################
#'
#' .onLoad getOption package settings
#'
#' @param libname defunct
#' @param pkgname defunct
#'
#' @return invisible()
#'
#' @examples
#' getOption("pcds.name")
.onLoad <- function(libname, pkgname) {
  op <- options()
  op.pcds <- list(
    #pcds.path = "~/R-dev",
    pcds.install.args  = "",
    pcds.name          = "Elvan Ceyhan",
    pcds.desc.author   = "Elvan Ceyhan <elvanceyhan@gmail.com> [aut, cre]",
    pcds.desc.license  = "GPL-2",
    pcds.desc.suggests = NULL,
    pcds.desc          = list()
  )
  toset <- !(names(op.pcds) %in% names(op))
  if (any(toset)) options(op.pcds[toset])

  invisible()
}

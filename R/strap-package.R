# devtools::build_win(version = "R-devel")

#' Measuring Morphological Diversity and Evolutionary Tempo
#' 
#' A series of functions for stratigraphic analysis of phylogenetic trees.
#' 
#' @name strap-package
#'
#' @aliases strap
#'
#' @docType package
#'
#' @author Mark A. Bell <mark.bell521@@gmail.com>>
#'
#' @references
#'
#' Bell, M. A. and Lloyd, G. T., 2015. strap: an R package for plotting phylogenies against stratigraphy and assessing their stratigraphic congruence. \emph{Palaeontology}, \bold{58}, 379-389.
#'
#' @keywords congruence,geology,phylogeny,stratigraphy,timescale
#'
#' @examples
#' 
#' # Calculate stratigraphic fit measures treating ages as ranges:
#' fit.to.strat.1 <- StratPhyloCongruence(trees = Dipnoi$tree,
#'   ages = Dipnoi$ages, rlen = 0, method = "basic", samp.perm=5,
#'   rand.perm = 5, hard = TRUE, randomly.sample.ages = FALSE,
#'   fix.topology = FALSE, fix.outgroup = TRUE,
#'   outgroup.taxon = "Psarolepis_romeri")
#'
#' # Show just the output for the input tree(s)
#' fit.to.strat.1$input.tree.results
#'
#' @exportPattern "^[[:alpha:]]+"
#' @import ape
#' @import geoscale

#' @importFrom graphics par points rect segments text
#' @importFrom grDevices rgb
#' @importFrom pbapply pblapply
#' @importFrom stats median pnorm runif sd
#' @importFrom utils data
NULL

#' Phylogeny and age data for the Asaphidae
#'
#' Phylogeny (162 most parsimonious trees) and age data for Asaphidae genera (Trilobita, Asaphida) taken from Bell and Braddy (2012).
#'
#'
#' @name Asaphidae
#' @docType data
#' @format A list containing 162 trees (\code{$tree}) and a matrix of first and last appearances (\code{$ages}).
#' @references Bell, M. A. and Braddy, S. J., 2012, Cope’s rule in the Ordovician trilobite Family Asaphidae (Order Asaphida): patterns across multiple most parsimonious trees: \emph{Historical Biology}, \bold{24}, 223-230.
#' @keywords datasets
NULL






#' Phylogeny and age data for dipnoans (lungfish)
#'
#' Phylogeny (first most parsimonious tree) and age data for lungfish (Osteichthyes, Sarcopterygii, Dipnoi) taken from Lloyd et al. (2012).
#'
#'
#' @name Dipnoi
#' @docType data
#' @format A list containing a tree (\code{$tree}) and a matrix of first and last appearances (\code{$ages}).
#' @references Lloyd, G. T., Wang, S. C. and Brusatte, S. L., 2012. Identifying heterogeneity in rates of morphological evolution: discrete character change in the evolution of lungfish (Sarcopterygii; Dipnoi). \emph{Evolution}, \bold{66}, 330-348.
#' @keywords datasets
NULL





#' British regional stages for the Ordovician
#' 
#' The stratigraphic ranges for the British stages of the Ordovician.
#' 
#' 
#' @name UKzones
#' @docType data
#' @format A matrix containing the start and end ages for the British stage subdivisions of the Ordovician.
#' @references Webby, B.D., Paris, F., Droser, M.L., and Percival, I.G. (editors), 2004. \emph{The Great Ordovician Biodiversity Event}. Columbia University Press, New York, 496 pp.
#' @keywords datasets
NULL




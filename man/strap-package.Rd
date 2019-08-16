\name{strap-package}
\alias{strap-package}
\alias{strap}
\docType{package}
\title{Stratigraphic Tree Analysis for Palaeontology}
\description{A series of functions for stratigraphic analysis of phylogenetic trees.}
\details{
\tabular{ll}{
Package: \tab strap\cr
Type: \tab Package\cr
Version: \tab 1.4-1\cr
Date: \tab 2014-12-03\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}

}

\author{Mark A. Bell <mark.bell521@gmail.com> and Graeme T. Lloyd <graemetlloyd@gmail.com>}

\references{
		
Wills, M. A., 1999, Congruence between phylogeny and stratigraphy: randomization tests and the gap excess ratio. Systematic Biology, 48(3), 559-580.
    
}

\keyword{phylogeny}

\examples{
# Calculate stratigraphic fit measures treating ages as ranges:
fit.to.strat.1 <- StratPhyloCongruence(trees=Dipnoi$tree, ages=Dipnoi$ages, rlen=0,
	method="basic", samp.perm=5, rand.perm=5, hard=TRUE,
	randomly.sample.ages=FALSE, fix.topology=FALSE, fix.outgroup=TRUE)
	
fit.to.strat.1$input.tree.results # Show just the output for the input tree(s)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DateNodeHedman.R
\name{DateNodeHedman}
\alias{DateNodeHedman}
\title{Calculates branch lengths for a topology}
\usage{
DateNodeHedman(tnodes, t0, resolution)
}
\arguments{
\item{tnodes}{The sequence of outgroup ages to the target node.}

\item{t0}{The arbitrary lower stratigraphic bound.}

\item{resolution}{The number of steps to take between the FAD and the lower stratigraphic bound.}
}
\value{
A phylo object with branch lengths scaled to time and the root age stored as $root.time.
}
\description{
Calculates branch lengths for a topology given a tree and age data for the tips.
}
\details{
The basic method (Norell 1992; Smith 1994) of dating a phylogenetic tree of fossil occurrences in palaeontology has been to make each internal node the age of its oldest descendant. In practical terms this means at least half or the branches in a fully bifurcating tree will have a duration of zero million years, as a hypothetical ancestor and its immediate descendant will have the same age, creaing a major problem for a variety of rate-based approaches where bracnh duration is the denominator.

Early solutions to this problem relied on adding some arbitrary value to each branch in order to enforce non-zero durations. However, more recently Ruta et al. (2006) argued for an approach that first dated the tree using the basic approach then, working from tip-to-root, whenever a zero duration branch was encountered it was assigned a share of the time available from the first directly ancestral branch of positive length. The size of this share is decided by some measure of evolutionary change along that branch. Ruta et al. (2006) used patristic dissimilarity (Wagner 1997), but conceivably any measure could be used. This approach was modified slightly by Brusatte et al. (2008), who preferred equal sharing. This has a couple of benefits over Ruta et al. (2006). Firstly, it avoids zero-length branches entirely - these could still happen with the Ruta et al. 2006 approach, as if no change occurs along a branch it gets zero share of any time. Secondly, it opens up the dating approach to trees without meaningful branch lengths, such as supertrees.

An undiscussed problem with the Ruta et al. (2006), and by extension the Brusatte et al. (2008) approach, concerns the inevitable zero-length branch at the base of the tree that has no preceding ancestral branch with which to share time. Here the obvious practical solution to this problem is implemented - to allow the user to pick a root length that the lowest branch(es) of the tree can share time with (Lloyd et al. 2012). Although selection of this value is potentially arbitrary, in most cases it will only effect a very small number of branches (potentially only a single branch). A recommended method for choosing root length is to use the difference between the oldest taxon in the tree and the age of the first outgroup to the tree that is older (ensuring a positive value).

Note that all three methods implemented here are effectively minimal approaches, in that they assume as little missing or unsampled history as possible. This is because they have their roots in maximum parsimony as an optimality criterion. Consequently the user should be aware that this function will likely return trees with relatively very short internal branch lengths, which may be a source of bias in subsequent analyses.

These approaches (with the exception of the Ruta method) are also implemented, along with others, in the timePaleoPhy function of the paleotree package.
}
\examples{

# Time-scale the lungfish tree using the "equal" method and a root length of 1 Ma:
time.tree <- DatePhylo(Dipnoi$tree, Dipnoi$ages, method = "equal", rlen = 1)

# Plot the tree with new branch lengths:
plot(time.tree, cex = 0.5)

}
\references{
Brusatte, S. L., Benton, M. J., Ruta, M. and Lloyd, G. T., 2008. Superiority, competition, and opportunism in the evolutionary radiation of dinosaurs. Science, 321, 1485-1488.
Lloyd, G. T., Wang, S. C. and Brusatte, S. L., 2012. Identifying heterogeneity in rates of morphological evolution: discrete character change in the evolution of lungfish (Sarcopterygii; Dipnoi). Evolution, 66, 330-348.
Norell, M. A., 1992. Taxic origin and temporal diversity: the effect of phylogeny. In: Extinction and Phylogeny, Novacek, M. J. and Wheeler, Q. D. (eds.). Columbia University Press, New York, p89-118.
Ruta, M., Wagner, P. J. and Coates, M. I., 2006. Evolutionary patterns in early tetrapods. I. Rapid initial diversification followed by decrease in rates of character change. Proceedings of the Royal Society B, 273, 2107-2111.
Smith, A. B., 1994. Systematics and the Fossil Record. Blackwell Scientific, London, 223pp.
Wagner, P. J., 1997. Patterns of morphologic diversification among the Rostroconchia. Paleobiology, 23, 115-15
}
\seealso{
cal3 in paleotree package
}
\author{
Matt Friedman \email{mfriedm@umich.edu} and Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}
\keyword{phylogeny}

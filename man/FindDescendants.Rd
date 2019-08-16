\name{FindDescendants}
\alias{FindDescendants}
\title{Finds the tip numbers descending from a specific node in a phylo object.}
\description{Finds the tip numbers descending from a specific node in a phylo object.}
\usage{FindDescendants(n, tree)}
\arguments{
  \item{n}{The node number.}
  \item{tree}{Tree as a phylo object.}
}


\examples{
# Find descendants of the root node in the lungfish tree:
FindDescendants(87, Dipnoi$tree)
}

\author{Graeme T. Lloyd \email{graemetlloyd@gmail.com}}


\keyword{phylogeny}

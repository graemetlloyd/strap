#' Finds the tip numbers descending from a specific node in a phylo object
#'
#' @description
#'
#' Finds the tip numbers descending from a specific node in a phylo object.
#'
#' @param n The node number.
#' @param tree Tree as a phylo object.
#'
#' @details
#'
#' A simple way to get the tips descending from a given node in a phylogenetic tree.
#'
#' @return A vector of the descendant tip numbers.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Find descendants of the root node in the lungfish tree:
#' FindDescendants(n = 87, tree = Dipnoi$tree)
#'
#' @export FindDescendants
FindDescendants <- function(n, tree) {

  # Vector to store ancestors (nodes):
  ancestors <- vector(mode = "numeric")
    
  # Vector to store descendnats (tips):
  descendants <- vector(mode = "numeric")
    
  # Keep going as long as there are nodes without recorded descendants (tips): 
  while (length(n) > 0) {
        
    # For each ancestral node find its descendants:
    for (i in 1:length(n)) ancestors <- c(tree$edge[grep(TRUE, tree$edge[, 1] == n[i]), 2], ancestors)
    
    # Reset n:
    n <- vector(mode = "numeric")
    
    # For each descendant:
    for (i in length(ancestors):1) {
      
      # If it is a tip store as a descendant:
      if (ancestors[i] <= Ntip(tree)) descendants <- c(descendants, ancestors[i])

      # If it is a node store as an ancestor:
      if (ancestors[i] > Ntip(tree)) n <- c(n, ancestors[i])
            
      # Remove ancestral node that has already has its descendants examined:
      ancestors <- ancestors[-i]

    }

  }
    
  # Output descendant tip numbers:
  return(descendants)

}

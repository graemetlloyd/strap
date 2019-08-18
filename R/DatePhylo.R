#' Calculates branch lengths for a topology
#'
#' @description
#'
#' Calculates branch lengths for a topology given a tree and age data for the tips.
#'
#' @param tree Tree as a phylo object.
#' @param ages A two-column matrix of taxa (rows) against First and Last Appearance Datums (FADs and LADs). Note that rownames should be the taxon names exactly as they appear in tree$tip.label and colnames should be "FAD" and "LAD". All ages should be in time before present.
#' @param rlen Root length. This must be greater than zero if using a method other than basic.
#' @param method The dating method used. Either basic (Norell 1992; Smith 1994), ruta (Ruta et al. 2006; requires input tree to have branch lengths) or equal (Brusatte et al. 2008).
#' @param add.terminal An optional to add the range of a taxon (FAD minus LAD) to terminal branch lengths.
#'
#' @details
#'
#' The basic method (Norell 1992; Smith 1994) of dating a phylogenetic tree of fossil occurrences in palaeontology has been to make each internal node the age of its oldest descendant. In practical terms this means at least half or the branches in a fully bifurcating tree will have a duration of zero million years, as a hypothetical ancestor and its immediate descendant will have the same age, creaing a major problem for a variety of rate-based approaches that use branch durations as a divisor.
#'
#' Early solutions to this problem relied on adding some arbitrary value to each branch in order to enforce non-zero durations. However, more recently Ruta et al. (2006) argued for an approach that first dated the tree using the basic approach then, working from tip-to-root, whenever a zero duration branch was encountered it was assigned a share of the time available from the first directly ancestral branch of positive length. The size of this share is decided by some measure of evolutionary change along that branch. Ruta et al. (2006) used patristic dissimilarity (Wagner 1997), but conceivably any measure could be used. This approach was modified slightly by Brusatte et al. (2008), who preferred equal sharing. This has a couple of benefits over Ruta et al. (2006). Firstly, it avoids zero-length branches entirely - these could still happen with the Ruta et al. 2006 approach, as if no change occurs along a branch it gets zero share of any time. Secondly, it opens up the dating approach to trees without meaningful branch lengths, such as supertrees.
#'
#' An undiscussed problem with the Ruta et al. (2006), and by extension the Brusatte et al. (2008) approach, concerns the inevitable zero-length branch at the base of the tree that has no preceding ancestral branch with which to share time. Here the obvious practical solution to this problem is implemented - to allow the user to pick a root length that the lowest branch(es) of the tree can share time with (Lloyd et al. 2012). Although selection of this value is potentially arbitrary, in most cases it will only effect a very small number of branches (potentially only a single branch). A recommended method for choosing root length is to use the difference between the oldest taxon in the tree and the age of the first outgroup to the tree that is older (ensuring a positive value).
#'
#' Note that all three methods implemented here are effectively minimal approaches, in that they assume as little missing or unsampled history as possible. This is because they have their roots in maximum parsimony as an optimality criterion. Consequently the user should be aware that this function will likely return trees with relatively very short internal branch lengths, which may be a source of bias in subsequent analyses.
#'
#' These approaches (with the exception of the Ruta method) are also implemented, along with others, in the timePaleoPhy function of the paleotree package.
#'
#' @return A phylo object with branch lengths scaled to time and the root age stored as $root.time.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Brusatte, S. L., Benton, M. J., Ruta, M. and Lloyd, G. T., 2008. Superiority, competition, and opportunism in the evolutionary radiation of dinosaurs. Science, 321, 1485-1488.
#' Lloyd, G. T., Wang, S. C. and Brusatte, S. L., 2012. Identifying heterogeneity in rates of morphological evolution: discrete character change in the evolution of lungfish (Sarcopterygii; Dipnoi). Evolution, 66, 330-348.
#' Norell, M. A., 1992. Taxic origin and temporal diversity: the effect of phylogeny. In: Extinction and Phylogeny, Novacek, M. J. and Wheeler, Q. D. (eds.). Columbia University Press, New York, p89-118.
#' Ruta, M., Wagner, P. J. and Coates, M. I., 2006. Evolutionary patterns in early tetrapods. I. Rapid initial diversification followed by decrease in rates of character change. Proceedings of the Royal Society B, 273, 2107-2111.
#' Smith, A. B., 1994. Systematics and the Fossil Record. Blackwell Scientific, London, 223pp.
#' Wagner, P. J., 1997. Patterns of morphologic diversification among the Rostroconchia. Paleobiology, 23, 115-15
#'
#' @seealso timePaleoPhy in paleotree package
#'
#' @keywords phylogeny
#'
#' @examples
#'
#' # Time-scale the lungfish tree using the "equal" method and a root length of 1 Ma:
#' time.tree <- DatePhylo(Dipnoi$tree, Dipnoi$ages, method = "equal", rlen = 1)
#'
#' # Plot the tree with new branch lengths:
#' plot(time.tree, cex = 0.5)
#'
#' @export DatePhylo
DatePhylo <- function(tree, ages, rlen = 0, method = "basic", add.terminal = FALSE) {
  
  # Stop if using Ruta method but not supplying a tree with branch lengths:
  if(is.null(tree$edge.length) && method == "ruta") stop("Tree has no branch lengths (required for Ruta method).")

  # Stop if taxon names do not match between tree and ages matrix:
  if(length(c(setdiff(tree$tip.label, rownames(ages)), setdiff(rownames(ages), tree$tip.label))) > 0) stop("Taxon names for tree and ages do not match.")

  # Stop if the ages matrix does not have columns named "FAD" and "LAD":
  if(length(match(c("FAD", "LAD"), colnames(ages))) < 2) stop("Ages matrix must have FAD and LAD (First and Last Appearance Datum) columns.")

  # Stop if any FAD is younger than its LAD:
  if(any(ages[, "FAD"] < ages[, "LAD"])) stop("FADs must all be at least as old as LADs.")

  # Stop if root length is negative:
  if(rlen < 0) stop("Root length cannot be negative.")

  # Stop if problem with root length:
  if(rlen == 0 && method != "basic") stop("If not using the basic method then rlen (root length) must be a positive value.")

  # Stop if method not available:
  if(method != "basic" && method != "ruta" && method != "equal") stop("Method must be one of basic, equal, or ruta.")

  # If method is not Ruta set all branch lengths equal to 1:
  if(is.null(tree$edge.length) || method != "ruta") tree$edge.length <- rep(1, length(tree$edge[, 1]))

  # Get node numbers to start:
  nodes <- c(ape::Ntip(tree) + 1):(ape::Nnode(tree) + ape::Ntip(tree))
  
  # Get starting node ages:
  node.ages <- unlist(lapply(as.list(nodes), function(x) max(ages[tree$tip.label[FindDescendants(n = x, tree = tree)], "FAD"])))

  # Get all ages (tips and nodes):
  all.ages <- as.vector(c(ages[tree$tip.label, "FAD"], node.ages))
  
  # Set root age to maximum age plus root length:
  all.ages[ape::Ntip(tree) + 1] <- all.ages[ape::Ntip(tree) + 1] + rlen
  
  # Create time-scaled tree:
  time.tree <- tree
  
  # Add branch lengths as time:
  time.tree$edge.length <- abs(apply(matrix(all.ages[tree$edge], ncol = 2), 1, diff))

  # Only continue if non basic dating option chosen:
  if(method != "basic") {
    
    # Keep going until there are no zero-length branches:
    while(min(time.tree$edge.length[grep(TRUE, tree$edge.length > 0)]) == 0) {
      
      # Record top zero-length branch encountered that is not also a zero change branch if using the Ruta method:
      share.branches <- intersect(grep(TRUE, time.tree$edge.length == 0), grep(TRUE, tree$edge.length > 0))[order(dist.nodes(tree)[(Ntip(tree) + 1), ][tree$edge[intersect(grep(TRUE, time.tree$edge.length == 0), grep(TRUE, tree$edge.length > 0)), 2]], decreasing = TRUE)[1]]
      
      # Keep going until there is a positive length branch:
      while(max(time.tree$edge.length[share.branches]) == 0) {
        
        # Find branches ancestral to those in memory:
        share.branches <- unique(c(share.branches, match(time.tree$edge[share.branches, 1], time.tree$edge[, 2])))

      }

      # Get total branch time:
      branch.time <- sum(time.tree$edge.length[share.branches])
      
      # Get number of branches to share:
      n.branches.to.share <- length(share.branches)
      
      # Get novel node ages (based on equal method):
      new.node.ages <- seq(from = range(all.ages[time.tree$edge[share.branches]])[1], to = range(all.ages[time.tree$edge[share.branches]])[2], by = sum(time.tree$edge.length[share.branches]) / n.branches.to.share)

      # Case if dating method is Ruta:
      if(method == "ruta") {
        
        # Get branch proportion based on branch lengths from input tree:
        branch.proportions <- ((tree$edge.length[share.branches] / sum(tree$edge.length[share.branches])))
        
        # Update novel node ages:
        new.node.ages <- c(new.node.ages[1], cumsum(branch.proportions[1:(length(branch.proportions) - 1)] * diff(range(new.node.ages))) + min(range(new.node.ages)), new.node.ages[length(new.node.ages)])

      }
      
      # Update node ages:
      all.ages[unique(as.vector(time.tree$edge[share.branches, 2:1]))] <- new.node.ages

      # Update branch lengths as time:
      time.tree$edge.length <- abs(apply(matrix(all.ages[tree$edge], ncol = 2), 1, diff))

    }

  }
  
  # Add ranges of taxa to terminal branch lengths if requested:
  if(add.terminal) time.tree$edge.length[match(1:Ntip(time.tree), time.tree$edge[, 2])] <- time.tree$edge.length[match(1:Ntip(time.tree), time.tree$edge[, 2])] + ages[time.tree$tip.label, "FAD"] - ages[time.tree$tip.label, "LAD"]

  # Get root age:
  root.age <- max(all.ages)

  # Store root age:
  time.tree$root.time <- root.age
  
  # Output tree:
  return(time.tree)

}

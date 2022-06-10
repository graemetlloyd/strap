#' Calculates branch lengths for a topology
#'
#' @description
#'
#' Calculates branch lengths for a topology given a tree and age data for the tips.
#'
#' @param tnodes The sequence of outgroup ages to the target node.
#' @param t0 The arbitrary lower stratigraphic bound.
#' @param resolution The number of steps to take between the FAD and the lower stratigraphic bound.
#'
#' @details
#'
#' The basic method (Norell 1992; Smith 1994) of dating a phylogenetic tree of fossil occurrences in palaeontology has been to make each internal node the age of its oldest descendant. In practical terms this means at least half or the branches in a fully bifurcating tree will have a duration of zero million years, as a hypothetical ancestor and its immediate descendant will have the same age, creaing a major problem for a variety of rate-based approaches where bracnh duration is the denominator.
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
#' @author Matt Friedman \email{mfriedm@@umich.edu} and Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
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
#' @seealso cal3 in paleotree package
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
#' @export DateNodeHedman

# Hedman (2010) method for estimating confidence given ages of ougroups (tnodes is the sequence of ages, from oldest to youngest; t0 is the arbitrary lower stratigraphic bound; resolution is the number of steps to take between the FAD and the lower stratigraphic bound):
DateNodeHedman <- function(tnodes, t0, resolution) { # Outputs estimate with two-tailed 95% CIs
  
  # Check t0 is older than any other node age:
  if(!all(t0 > tnodes)) stop("t0 must be older than any tnode value.")
  
  # Store requested resolution:
  requested.resolution <- resolution
  
  # Function returns p.d.f. for node ages	(t0 is an arbitrary oldest age to consider, tnodes are the ages of oldest known fossil stemming from each node, and tsteps is a vector of arbitrary time steps on which the p.d.f. is calculated):
  nodeage <- function(tsteps, tnodes, t0) {
    
    # Get number of outgroups (tnodes):
    nn  <- length(tnodes)
    
    # Get number of time steps (resolution):
    nt  <- length(tsteps)
    
    # Initialize array:
    pnodes  <- matrix(0, nn, nt)
    
    # First get pdf for node 1, at discrete values:
    ii  <- which(tsteps > t0 & tsteps < tnodes[1])
    
    # Assume uniform distribution for oldest node:
    pnodes[1, ii] <- 1.0 / length(ii)
    
    # Cycle through remaining nodes:
    for (i in 2:nn) {
      
      # Cycle through series of time steps:
      for (j in 1:nt) {
        
        # Initialize vector:
        p21  <- rep(0, nt)
        
        # Get p.d.f. for ith node ay discrete values:
        ii <- which(tsteps >= tsteps[j] & tsteps < tnodes[i])
        
        
        p21[ii]  <- 1.0 / length(ii) * pnodes[(i - 1), j]
        
        # Conditional probability of this age, given previous node age, times probability of previous node age, added to cumulative sum:
        pnodes[i, ii]  <- pnodes[i, ii] + p21[ii]
        
      }
      
    }
    
    # Get just the p.d.f. vector for the youngest node:
    out <- pnodes[nn, ]
    
    # Return output:
    return(out)
    
  }
  
  # Calculate c.d.f. from p.d.f., find median and 95% credibility limits:
  HedmanAges <- function(tnodes, t0, resolution) {
    
    # Store input resolution for later reuse:
    old.resolution <- resolution
    
    # Get uniformly spaced p-values (including 95% limits and median) for calculating CIs:
    CIs <- sort(unique(c(0.025, 0.975, 0.5, c(1:old.resolution) * (1 / old.resolution))))
    
    # Get enw resolution (may be longer by adding limits and/or median):
    resolution <- length(CIs)
    
    # Get oldest possible age to consider:
    first <- t0
    
    # Get youngest possible age to consider:
    last <- min(tnodes)
    
    # Make tnodes negative for Hedman function:
    tnodes  <- -tnodes
    
    # Get uniformly spaced time steps for Hedman function:
    tsteps <- seq(-first, -last, length=resolution)
    
    # Make t0 negative for Hedman function:
    t0 <- -t0
    
    # Run Hedman function to get p.d.f.:
    vector <- nodeage(tsteps, tnodes, t0)
    
    # Convert from p-value scale to age scale:
    integral.vector <- vector * ((abs(first - last)) / resolution)
    
    # Get sum of probabilities in order to re-scale as p.d.f.:
    probability.sum <- sum(integral.vector)
    
    # Get re-scaled p-values for CIs:
    ps.CIs <- probability.sum * CIs
    
    # Cretae empty vector to store dates:
    date.distribution <- vector(mode="numeric")
    
    # Set initial value:
    value <- 0
    
    # For each re-scaled p-value:
    for(i in length(ps.CIs):1) {
      
      # Update value:
      value <- value + integral.vector[i]
      
      # If in age window:
      if(length(which(value >= ps.CIs)) > 0) {
        
        # Store dates:
        date.distribution <- c(date.distribution, rep(tsteps[i], length(which(value >= ps.CIs))))
        
        # Update re-scaled p-values:
        ps.CIs <- ps.CIs[-grep(TRUE, value >= ps.CIs)]
        
      }
      
    }
    
    # Add t0 at end if length is short:
    while(length(date.distribution) < resolution) date.distribution <- c(date.distribution, t0)
    
    # Find median of distribution:
    Best.guess <- date.distribution[CIs == 0.5]
    
    # Get upper 95% CI:
    Best.guess.lower <- date.distribution[CIs == 0.975]
    
    # Get lower 95% CI:
    Best.guess.upper <- date.distribution[CIs == 0.025]
    
    # Combine results and reverse sign to give palaeo ages:
    results <- list(-Best.guess, -Best.guess.lower, -Best.guess.upper, -date.distribution[match(c(c(1:old.resolution) * (1 / old.resolution)), CIs)])
    
    # Name subvariables:
    names(results) <- c("Best.guess", "Best.guess.lower", "Best.guess.upper", "Age.distribution")
    
    # Return result
    return(results)
    
  }
  
  # Get Hedman ages:
  out <- HedmanAges(tnodes, t0, resolution)
  
  # If Hedman ages represent less than five unique values (leading to downstream problem of flat distributions):
  while(length(unique(out$Age.distribution[round(seq(1, length(out$Age.distribution), length.out=requested.resolution))])) < 5) {
    
    # Double the resolution size:
    resolution <- resolution * 2
    
    # Get Hedman ages:
    out <- HedmanAges(tnodes, t0, resolution)
    
  }
  
  # Update output:
  out$Age.distribution <- out$Age.distribution[round(seq(1, length(out$Age.distribution), length.out=requested.resolution))]
  
  # Return output:
  return(out)
  
}

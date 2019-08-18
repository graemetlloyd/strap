#' Calculates fit to stratigraphy metrics for a set of tree(s)
#'
#' @description
#'
#' Calculates SCI, RCI, MSM*, and GER for a number of topologies.
#'
#' @param trees Input tree(s) as either a phylo or multiphylo object.
#' @param ages A two-column matrix of taxa (rows) against First and Last Appearance Datums (FADs and LADs) to be passed to DatePhylo. Note that rownames should be the taxon names exactly as they appear in tree$tip.label and colnames should be "FAD" and "LAD". All ages should be in time before present.
#' @param rlen Root length, to be passed to DatePhylo.
#' @param method Tree dating method, to be passed to DatePhylo.
#' @param samp.perm Number of sampled trees to be produced by resolving polytomies and/or drawing random dates for the tips for the input trees.
#' @param rand.perm Number of random trees to be produced in calculating probabilities for the input trees, and (if used) the sampled trees.
#' @param hard Whether to treat polytomies as hard or soft. If FALSE polytomies are resolved randomly.
#' @param randomly.sample.ages Whether to treat FAD and LAD as a range (randomly.sample.ages = FALSE) or an uncertainty (randomly.sample.ages = TRUE). If the latter then two ages are randomly sampled from the range and these are used as the FAD and LAD.
#' @param fix.topology Whether to allow tree shape to be random (fix.topology = FALSE) or to reflect the tree shape of the input tree(s) (fix.topology = TRUE).
#' @param fix.outgroup Whether to force the randomly generated trees to share the outgroup of the first input tree (fix.outgroup = TRUE) or not (fix.outgroup = FALSE)
#' @param outgroup.taxon The outgroup taxon name (only applicable if fix.outgroup = TRUE).
#'
#' @details
#'
#' Cladograms of fossil taxa make explicit predictions about the successive appearance of taxa in the fossil record that can be compared with their observed stratigraphic ranges. Several methods have been developed to quantify this "fit to stratigraphy" of phylogenetic hypotheses, and these can be assessed in both a statistic (measuring the apparent strength of this congruence) and an associated significance test (p-value) based on generating random topologies for the same taxon set.
#'
#' This function produces both values for all four main metrics: the Stratigraphic Consistency Index (SCI; Huelsenbeck 1994), the Relative Consistency Index (RCI; Benton and Storrs 1994), the Manhattan Stratigraphic Measure (MSM*; Siddall 1998; Pol and Norell 2001), and the Gap Excess Ratio (GER; Wills 1999).
#'
#' \bold{SCI - Stratigraphic Consistency Index}
#'
#' The SCI works by assessing the "consistency" of nodes. A node is considered stratigraphically consistent if its oldest descendant is the same age or younger than the oldest descendant of the preceding node. The SCI is thus given simply as:
#'
#' SCI = C/N
#'
#' Where C is the sum of all the consistent nodes and N is the total number of nodes - 1. (As there is no node preceding the root there is no basis on which to estimate its consistency.) This value can range from zero (maximally inconsistent) to one (maximally consistent). However, a potential criticism of the SCI is that a high value may be returned when in fact a single inconsistent node may represent a very large amount of missing history (measured in unsampled units or millions of years), whereas a low SCI may represent relatively few unsampled units or millions of years.
#'
#' \bold{RCI - Relative Completeness Index}
#'
#' The RCI was the first method to explicitly account for the absolute amount of missing data implied by the tree. This figure is usually expressed as the Minimum Implied Gap (MIG), a term also used by both the MSM and GER (see below), and corresponds to the sum of the branch lengths excluding the duration of the terminals (the observed ranges of the taxa). The RCI expresses the MIG as a proportion of the sum of the observed ranges (Simple Range Length; SRL) of the taxa converted to a percentage:
#'
#' RCI = (1 − (MIG/SRL)) ∗ 100percent
#'
#' Importantly this value is not confined to a 0 to 100 percent scale, and can have both negative values and values greater than 100 percent, which can make it difficult to interpret.
#'
#' \bold{MSM - Manhattan Stratigraphic Measure}
#'
#' The MSM was the first method to account for both the absolute MIG and range on a confined zero to one scale. It is expressed as:
#'
#' MSM = Lm/L0
#'
#' Where L0 is the length of the tree expressed by optimising times of first appearance on to the tree as a Sankoff character and taking the total length. Lm represents the same process, but for the optimal possible tree given the same set of first appearances. However, Pol and Norell (2001) noted a critical flaw in this approach, specifically that the Sankoff optimisation is reversible, meaning that nodes in the topology are allowed to be younger than their descendants, leading in some cases to a poor fit to stratigraphy being perceived as a good fit. Instead they suggest modifying the character step matrix to make the cost of reversals effectively infinite and hence impossible. Thus the values for L0 and Lm are modified accordingly. This approach they termed MSM* and is the implementation of MSM used here. This statistic can be expressed as:
#'
#' MSM∗ = Gmin/MIG
#'
#' Where Gmin represents the MIG for the tree with the optimal fit to stratigraphy. In effect this is a completely unbalanced tree where the youngest pair of taxa are the most deeply nested and successive outgroups represent the next oldest taxon. Theoretically MSM* ranges from one (the best fit as the observed tree is the maximally consistent tree) to zero (the least optimal tree). However, in effect no tree can have a value of zero as its MIG would have to be equal to infinity.
#'
#' \bold{GER - Gap Excess Ratio}
#'
#' The GER represents a method that accounts for MIG, ranges from zero to one, and the best and worst fits to stratigraphy are both practically realisable. It can be expressed as:
#'
#' GER = 1 − ((MIG − Gmin)/(Gmax − Gmin))
#'
#' Where Gmax represents the MIG of the tree with the worst possible fit to stratigraphy. This is in effect any topology where the oldest taxon is the most deeply nested such that every clade in the tree contains it and hence must be minimally that old.
#'
#' \bold{P-values}
#'
#' In isolation all four methods suffer from an inability to reject the null hypothesis that an apparent good fit to stratigraphy may be generated by chance alone. In practice this can be tested by generating a set of random topologies, calculating the fit to stratigraphy measure, and then either fitting a normal distribution to the resulting values (to get an estimated p-value) or assessing the relative position of the MIG of the observed (and sampled) tree(s) to get an absolute p-value. (Note that the assumption of normality may not always hold and the former approach should be used at the user’s discretion. However, it should be noted that for the SCI, MSM*, and GER p-values are calculated after first transforming the data by taking the arcsine of the square root of each value.) The reason for having two sets of p-values is that if the observed trees fall completely outside the range of the random topologies they will be given an extreme p-value (0 or 1) that may be misleading. In such cases the estimated value may be more accurate.
#'
#' P-values should be interpreted as the probability of the null: that the observed tree(s) have an equal or worse fit to stratigraphy than the sample of random trees. Thus if the p-values are very small the user can reject the null hypothesis in favour of the alternative: that the observed tree(s) have a better fit to stratigraphy than expected by chance alone.
#'
#' \bold{Modifications of the GER}
#'
#' More recently Wills et al. (2008) introduced two new versions of the GER that take advantage of the distribution of MIGs from the set of randomly generated topologies. The first of these (GERt) uses the extreme values of the random topologies as modified versions of Gmax and Gmin, termed Gtmax and Gtmin respectively. GERt is thus expressed as:
#'
#' GERt = 1 − ((MIG − Gtmin)/(Gtmax − Gtmin))
#'
#' In practice the MIG of the observed tree(s) may fall outside of these ranges so here a correction factor is employed so that any value below zero is corrected to zero, and any value above one is corrected to one. An additional stipulation for GERt is that the overall tree topology is fixed and only the taxa themselves are shuffled. This is to give a more realistic set of random topologies as there are known biases towards unbalanced trees in many palaeontological data sets. Here this is implemented by selecting the fix.topology=TRUE option. However, here GERt can also be calculated when fix.topology=FALSE.
#'
#' A second modification of GER is to use the position of the observed tree(s) in the sample of randomly generated topologies, such that:
#'
#' GER∗ = 1 − (F ractionof distribution <= MIG)
#'
#' Thus if the MIG of the observed tree(s) is less than any randomly generated topology GER* will be one (maximally optimal fit) and if it worse than any of the randomly generated topologies it will be zero (maximally suboptimal fit).
#'
#' Note: it is recommended that you use a large number of random topologies in order to get reliable values for GERt and GER* using rand.perm=N. Wills et al. (2008) used 50000, but the user should note that for many real world examples such values will take many hours to run.
#'
#' Wills et al. (2008) also introduced the notion of referring to intervals sampled rather than absolute time by recasting MIG as MIGu: the sum of ghost ranges for intervals of unit length. Although not directly implemented here this can be done manually by converting the time values (in Ma) used to simple unit counts such that FADs and LADs of taxa are given as numbered time bins (the youngest being 1 and the oldest N, where there are N time bins).
#'
#' \bold{Polytomies and age uncertainties}
#'
#' Alongside the input trees the user can also create an additional set of sampled trees based on the input trees. This option is automatically implemented when choosing either hard = FALSE or randomly.sample.ages = TRUE, and the total number of permutations to perform dictated by samp.perm = N. This process works by first sampling from the set of input tree(s) and then randomly resolving any polytomies (if hard = FALSE) to ensure all sampled trees are fully dichotomous. (At present the function does not allow the various options laid out in Boyd et al. 2011, but the user can achieve this effect by modifying the input trees themselves.) Then if randomly.sample.ages = TRUE the FAD and LAD are treated as bounds of a uniform distribution which is sampled at random. This allows the user to get results for a set of trees that account for uncertainty in dating (as outlined in Pol and Norell 2006). (Note that two dates are picked for each taxon to avoid the problem of having an SRL of zero that would cause a divide by zero error for the RCI metric.)
#'
#' All fit to stratigraphy measures calculated for the input trees are then repeated for the sampled trees. However, if both hard = TRUE and randomly.sample.ages = FALSE (the defaults) no set of sampled trees will be created.
#'
#' In all cases when using the function users will see progress bars that indicate the general progress through the various sets of trees (input, sampled, and randomly generated). This serves as a useful indicator of the time it will take for the function to finish. Here default values for samp.perm and rand.perm are both set at 1000, but the user may wish to lower these (to decrease calculation time) or increase them (to enhance accuracy).
#'
#' \bold{Additional options}
#'
#' Note that because this function uses DatePhylo the user has the option of using different tree dating algorithms than the basic method (equivalent to the basic method in the paleotree package) employed in all the published studies cited above (and the default option here). The dating method used will apply to all trees generated, including the input, sampled, and randomly generated topologies. In all cases the time-scaled trees are returned with the function output.
#'
#' A final option (fix.outgroup = TRUE) allows the user to always use the same outgroup taxon (supplied as outgroup.taxon) for all randomly generated topologies. Because the outgroup will often be the oldest taxon and its position in the input topologies is not allowed to vary letting it do so in the random topologies may lead to inferring a better fit to stratigraphy for the observed tree(s) than is fair. Fixing the outgroup thus ameliorates this potential bias and is the default option here.
#'
#' @return
#'
#' \item{input.tree.results}{A matrix with a row for each input tree and columns indicating the values for SCI, RCI, GER and MSM* and their estimated probabilities assuming a normal distribution (est.p.SCI, est.p.RCI, est.p.GER, and est.p.MSM*) as well as GERt, GER*, MIG, and p.Wills (their probability as position within the MIGs of the random topologies).}
#' \item{samp.permutation.results}{If used, a matrix with a row for each sampled tree (up to samp.perm) and columns indicating the values for SCI, RCI, GER and MSM* and their estimated probabilities assuming a normal distribution (est.p.SCI, est.p.RCI, est.p.GER, and est.p.MSM*) as well as GERt, GER*, MIG, and p.Wills (their probability as position within the MIGs of random topologies).}
#' \item{rand.permutations}{A matrix with a row for each randomly generated tree (up to rand.perm) and columns indicating the values for SCI, RCI, GER, MSM*, and MIG.}
#' \item{input.trees}{The input tree(s) as a phylo or multiphylo object, with branches scaled to time according to the input values passed to \link{DatePhylo}.}
#' \item{samp.trees}{The sampled tree(s) as a phylo or multiphylo object, with branches scaled to time according to the input values passed to \link{DatePhylo}.}
#' \item{rand.trees}{The randomly generated tree(s) as a phylo or multiphylo object, with branches scaled to time according to the input values passed to \link{DatePhylo}.}
#'
#' @author Mark A. Bell \email{mark.bell521@@gmail.com} and Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Bell, M. A. and Lloyd, G. T., 2015. \code{strap}: an R package for plotting phylogenies against stratigraphy and assessing their stratigraphic congruence. \emph{Palaeontology}, \bold{58}, 379-389.
#'
#' Benton, M. J. and Storrs, G. W., 1994. Testing the quality of the fossil record: palaeontological knowledge is improving. \emph{Geology}, \bold{22}, 111-114.
#'
#' Boyd, C. A., Cleland, T. P., Marrero, N. L. and Clarke, J. A., 2011. Exploring the effects of phylogenetic uncertainty and consensus trees on stratigraphic consistency scores: a new program and a standardized method. \emph{Cladistics}, \bold{27}, 52-60.
#'
#' Huelsenbeck, J. P., 1994. Comparing the stratigraphic record to estimates of phylogeny. \emph{Paleobiology}, \bold{20}, 470-483.
#'
#' Pol, D. and Norell, M. A., 2001. Comments on the Manhattan Stratigraphic Measure. \emph{Cladistics}, \bold{17}, 285-289.
#'
#' Pol, D. and Norell, M. A., 2006. Uncertainty in the age of fossils and the stratigraphic fit to phylogenies. \emph{Systematic Biology}, \bold{55}, 512-521.
#'
#' Siddall, M. E., 1998. Stratigraphic fit to phylogenies: a proposed solution. \emph{Cladistics}, \bold{14}, 201-208.
#'
#' Wills, M. A., 1999. Congruence between phylogeny and stratigraphy: randomization tests and the Gap Excess Ratio. \emph{Systematic Biology}, \bold{48}, 559-580.
#'
#' Wills, M. A., Barrett, P. M. and Heathcote, J. F., 2008. The modified Gap Excess Ratio (GER*) and the stratigraphic congruence of dinosaur phylogenies. \emph{Systematic Biology}, \bold{57}, 891-904.
#'
#' @keywords congruence,phylogeny,stratigraphy
#'
#' @examples
#'
#' # Calculate stratigraphic fit measures treating ages as ranges
#' # (permutation numbers used are lower than recommended for standard use):
#' fit.to.strat.1 <- StratPhyloCongruence(trees = Dipnoi$tree, ages = Dipnoi$ages, rlen = 0,
#'   method = "basic", samp.perm = 5, rand.perm = 5, hard = TRUE,
#'   randomly.sample.ages = FALSE, fix.topology = TRUE, fix.outgroup = TRUE,
#'   outgroup.taxon = "Psarolepis_romeri")
#'
#' # View all output:
#' fit.to.strat.1
#'
#' # Show output options:
#' names(fit.to.strat.1)
#'
#' # Show just the output for the input tree(s):
#' fit.to.strat.1$input.tree.results
#'
#' # Calculate stratigraphic fit measures treating ages as uncertainties
#' # (permutation numbers used are lower than recommended for standard use):
#' fit.to.strat.2 <- StratPhyloCongruence(trees = Dipnoi$tree, ages = Dipnoi$ages, rlen = 0,
#'   method = "basic", samp.perm = 10, rand.perm = 10, hard = TRUE,
#'   randomly.sample.ages = TRUE, fix.topology = TRUE, fix.outgroup = TRUE,
#'   outgroup.taxon = "Psarolepis_romeri")
#'
#' @export StratPhyloCongruence
StratPhyloCongruence <- function(trees, ages, rlen = 0, method = "basic", samp.perm = 1000, rand.perm = 1000, hard = TRUE, randomly.sample.ages = FALSE, fix.topology = TRUE, fix.outgroup = TRUE, outgroup.taxon = NULL) {
  
  # SRL only matters for taxa with ages so add option to not do double sampling from sample permutation of ages
  # Instead of randomly shuffling tips randomly sample trees from bifurcations of consensus? Could be better at finding similar trees
  # Check for repeating trees and remove?
  # Also what about when tip count is super low?
  # More checks of data format
  # CHECK TREES HAVE SAME NUMBER OF TIPS AND SAME NAMES
  # ADD USER INPUT OPTIONS FOR RANDOM TREES SAMPLE?
  # MAKE SCI OPTIONAL AS SLOWEST PART OF FUNCTION
  # MAKE SRL 1 IF EVER GETS TO ZERO?
  
  # Subfunction to calculate the SCI:
  StratigraphicConsistencyIndex <- function(tree, ages) {
    
    # Build vector of internal node numbers excluding the root:
    NodeNumbers <- (Ntip(tree) + 2):(Ntip(tree) + tree$Nnode)
    
    # Count consistent nodes:
    ConsistentNodeCount <- sum(unlist(lapply(NodeNumbers, function(x) {CurrentNodeDescendants <- tree$tip.label[FindDescendants(x, tree)]; PreviousNodeDescendants <- setdiff(tree$tip.label[FindDescendants(tree$edge[tree$edge[, 2] == x, 1], tree)], CurrentNodeDescendants); max(ages[CurrentNodeDescendants, "FAD"]) >= max(ages[PreviousNodeDescendants, "FAD"])})))
    
    # Return SCI:
    return(ConsistentNodeCount / (tree$Nnode - 1))
    
  }
  
  # If a single tree convert to list to ensure a standard format for input tree(s):
  if(class(trees) == "phylo") trees <- list(trees)
  
  # Strip names off trees:
  names(trees) <- NULL
  
  # Get total number of input trees:
  Ntrees <- length(trees)
  
  # Check that if using Ruta time-scaling method all trees include branch lengthsa nd stop and warn user if not:
  if(!length(unlist(lapply(trees, function(x) grep("edge.length", names(x))))) == length(trees) && method == "ruta") stop("If using \"ruta\" method then all trees must include branch lengths.")

  # Make sure permutation numbers are not negative:
  if(samp.perm <= 0 || rand.perm <= 0) stop("Number of permutations must be positive.")
  
  # Check that if fixing outgroup that an outgroup taxon is specified and stop and warn user if not:
  if(fix.outgroup && is.null(outgroup.taxon)) stop("If using fix.outgroup = TRUE must specify an outgroup.taxon.")
  
  # If outgroup taxon is being used check it is found in all trees and stop and warn user if not:
  if(!is.null(outgroup.taxon)) if(!all(unlist(lapply(trees, function(x) length(intersect(outgroup.taxon, x$tip.label)))) == 1)) stop("outgroup.taxon must appear exactly once in all trees.")
  
  # Set null variables that may not be used later but are still included in the output:
  samp.permutations <- samp.trees <- NULL
  
  # If sampling tree set SamplingTrees to TRUE else set to FALSE:
  SamplingTrees <- ifelse(hard == FALSE || randomly.sample.ages == TRUE, TRUE, FALSE)
  
  # Create matrix to store output with extra column for GERt:
  input.permutations <- matrix(ncol = 15, nrow = length(trees), dimnames = list(paste("tree_", c(1:length(trees)), sep = ""), c("SRL", "MIG", "GMax", "GMin", "SCI", "RCI", "GER", "MSM*", "est.p.SCI", "est.p.RCI", "est.p.GER", "est.p.MSM*", "GER*", "GERt", "p.Wills")))
  
  # Create matrix to store random permutation results:
  rand.permutations <- matrix(nrow = rand.perm, ncol = 8, dimnames = list(c(), c("SRL", "MIG", "GMax", "GMin", "SCI", "RCI", "GER", "MSM*")))
  
  # If sampling trees create matrix to store permutation results with extra column for GERt::
  if(SamplingTrees) samp.permutations <- matrix(nrow = samp.perm, ncol = 15, dimnames = list(c(), c("SRL", "MIG", "GMax", "GMin", "SCI", "RCI", "GER", "MSM*", "est.p.SCI", "est.p.RCI", "est.p.GER", "est.p.MSM*", "GER*", "GERt", "p.Wills")))
  
  # Calculate sum of stratigraphic ranges of taxa (SRL):
  SRL <- sum(ages[, "FAD"] - ages[, "LAD"])
  
  # If sum of stratigraphic ranges is zero (all FADs are equal to their LADs):
  if(SRL == 0) {
  
    # If user has requested random sampling of ages method:
    if(randomly.sample.ages == TRUE) {
    
      # Warn user:
      cat("Can not perform random sampling of ages method as sum of stratigraphic ranges is zero.\nFor RCI sum of stratigraphic ranges is set to one.\n\n")
      
      # Set randomly sampled ages as false:
      randomly.sample.ages <- FALSE
      
    # If randomly sampled ages is set to false:
    } else {

      # Warn user:
      cat("Sum of stratigraphic ranges is set to one (to avoid dividing by zero).\n\n")
    
    }
    
    # Set SRL as one:
    SRL <- 1
    
  }
  
  # If user wants to fix the topology as in GERt:
  if(fix.topology) {
    
    # Create random trees by resampling input trees:
    rand.trees <- trees[sample(1:length(trees), rand.perm, replace = TRUE)]
    
    # If not using hard = TRUE then make sure trees are all fully bifurcating:
    if(!hard) rand.trees <- lapply(rand.trees, multi2di)
    
    # If fixing the outgroup taxon:
    if(fix.outgroup) {
      
      # Shuffle all tip labels except the outgroup taxon:
      rand.trees <- lapply(rand.trees, function(x) {x$tip.label[-(x$tip.label == outgroup.taxon)] <- sample(x$tip.label[-(x$tip.label == outgroup.taxon)]); x})
      
    # If not fixing the outgroup taxon:
    } else {
      
      # Shuffle all tip labels including the outgroup taxon:
      lapply(rand.trees, function(x) {x$tip.label <- sample(x$tip.label); x})
    
    }
    
    # If using Ruta method will need branch lengths adding to random trees:
    if(method == "ruta") {
      
      # Build vector of edge lengths from input trees to sample from for random trees:
      EdgeLengthsToSampleFrom <- unlist(lapply(trees, function(x) x$edge.length))
      
      # Assign edge lengths to random trees by randomly sampling from sampled tree edge lengths:
      rand.trees <- lapply(rand.trees, function(x) {x$edge.length <- sample(EdgeLengthsToSampleFrom, nrow(x$edge), replace = TRUE); x})
      
    }
    
  # If user does not want to fix the topology:
  } else {
  
    # Make random trees for permutations:
    rand.trees <- rmtree(rand.perm, Ntip(trees[[1]]))
    
    # Add real tip names as randomly shuffled values:
    rand.trees <- lapply(rand.trees, function(x) {x$tip.label <- sample(trees[[1]]$tip.label); x})
    
    # If fixing outgroup then reroot random trees on outgroup:
    if(fix.outgroup) rand.trees <- lapply(rand.trees, function(x) root(x, outgroup.taxon))
  
  }
  
  # If using the sample permutation loop:
  if(SamplingTrees) {
    
    # If polytomies are to be considered hard:
    if(hard == TRUE) {
      
      # Store randomly sampled tree:
      samp.trees <- trees[sample(1:Ntrees, size = samp.perm, replace = TRUE)]
      
    # If polytomies are to be considered soft:
    } else {
      
      # Store randomly sampled randomly bifurcated tree:
      samp.trees <- lapply(trees[sample(1:Ntrees, size = samp.perm, replace = TRUE)], function(x) ape::multi2di(x))
      
    }
    
  # If not permuting ages of polytomie resolutions:
  } else {
    
    # Reset samp.perm to zero (for correct counter values later):
    samp.perm <- 0
    
  }
  
  # If randomly sampling ages:
  if(randomly.sample.ages) {
    
    # Build individual ages matrix for each tree by randomly sampling twice between FAD and LAD:
    rand.trees <- lapply(rand.trees, function(x) {x$ages <- t(apply(ages, 1, function(y) sort(stats::runif(2, max = y[1], min = y[2]), decreasing  = TRUE))); colnames(x$ages) <- c("FAD", "LAD"); x})
    
    # Build individual ages matrix for each tree by randomly sampling once between FAD and LAD:
    #rand.trees <- lapply(rand.trees, function(x) {x$ages <- t(apply(ages, 1, function(y) rep(stats::runif(1, max = y[1], min = y[2]), 2))); colnames(x$ages) <- c("FAD", "LAD"); x})
    
  # If not randomly sampling ages:
  } else {
    
    # Build individual ages matrix for each tree by copying across ages matrix:
    rand.trees <- lapply(rand.trees, function(x) {x$ages <- ages; x})
    
  }
  
  # If sampling trees:
  if(SamplingTrees) {
    
    # If randomly sampling ages:
    if(randomly.sample.ages) {
      
      # Build individual ages matrix for each tree by randomly sampling twice between FAD and LAD:
      samp.trees <- lapply(samp.trees, function(x) {x$ages <- t(apply(ages, 1, function(y) sort(stats::runif(2, max = y[1], min = y[2]), decreasing  = TRUE))); colnames(x$ages) <- c("FAD", "LAD"); x})
      
      # Build individual ages matrix for each tree by randomly sampling once between FAD and LAD:
      #samp.trees <- lapply(samp.trees, function(x) {x$ages <- t(apply(ages, 1, function(y) rep(stats::runif(1, max = y[1], min = y[2]), 2))); colnames(x$ages) <- c("FAD", "LAD"); x})
      
    # If not randomly sampling ages:
    } else {
      
      # Build individual ages matrix for each tree by copying across ages matrix:
      samp.trees <- lapply(samp.trees, function(x) {x$ages <- ages; x})
      
    }
    
  }
  
  # Report to user current task being performed:
  cat("Time-scaling input trees...\nUNDERLAY")
  
  # Date every input tree using the date phylo options given:
  trees <- pbapply::pblapply(trees, function(x) DatePhylo(tree = x, ages = ages, rlen = rlen, method = method, add.terminal = FALSE))
  
  # Report to user current task being performed:
  cat("Time-scaling randomly generated trees...\nUNDERLAY")
  
  # Date every random tree using the date phylo options given:
  rand.trees <- pbapply::pblapply(rand.trees, function(x) DatePhylo(tree = x, ages = x$ages, rlen = rlen, method = method, add.terminal = FALSE))
  
  # If sampling trees:
  if(SamplingTrees) {
    
    # Report to user current task being performed:
    cat("Time-scaling randomly sampled trees...\nUNDERLAY")
    
    # Date every sampled tree using the date phylo options given:
    samp.trees <- pbapply::pblapply(samp.trees, function(x) DatePhylo(tree = x, ages = x$ages, rlen = rlen, method = method, add.terminal = FALSE))
    
  }
  
  # Report to user current task being performed:
  cat("Calculating Stratigraphic Consistency Index for input trees...\nUNDERLAY")
  
  # Calculate Stratigraphic Consistency Index for every input tree:
  input.permutations[, "SCI"] <- unlist(pbapply::pblapply(trees, function(x) StratigraphicConsistencyIndex(x, ages)))
  
  # Report to user current task being performed:
  cat("Calculating Stratigraphic Consistency Index for randomly generated trees...\nUNDERLAY")
  
  # Calculate Stratigraphic Consistency Index for every random tree:
  rand.permutations[, "SCI"] <- unlist(pbapply::pblapply(rand.trees, function(x) StratigraphicConsistencyIndex(x, x$ages)))
  
  # If sampling trees:
  if(SamplingTrees) {
    
    # Report to user current task being performed:
    cat("Calculating Stratigraphic Consistency Index for sampled trees...\nUNDERLAY")
    
    # Calculate Stratigraphic Consistency Index for every sampled tree:
    samp.permutations[, "SCI"] <- unlist(pbapply::pblapply(samp.trees, function(x) StratigraphicConsistencyIndex(x, x$ages)))
    
  }

  # Calculate Simple Range Length for every input tree:
  input.permutations[, "SRL"] <- unlist(lapply(trees, function(x) abs(sum(apply(ages, 1, diff)))))
  
  # Calculate Simple Range Length for every random tree:
  rand.permutations[, "SRL"] <- unlist(lapply(rand.trees, function(x) abs(sum(apply(x$ages, 1, diff)))))
  
  # If sampling trees calculate Simple Range Length for every sampled tree:
  if(SamplingTrees) samp.permutations[, "SRL"] <- unlist(lapply(samp.trees, function(x) abs(sum(apply(x$ages, 1, diff)))))
  
  # Calculate Minimum Implied Gap for every input tree:
  input.permutations[, "MIG"] <- unlist(lapply(trees, function(x) sum(x$edge.length)))
  
  # Calculate Minimum Implied Gap for every random tree:
  rand.permutations[, "MIG"] <- unlist(lapply(rand.trees, function(x) sum(x$edge.length)))
  
  # If sampling trees then calculate Minimum Implied Gap for every input tree:
  if(SamplingTrees) samp.permutations[, "MIG"] <- unlist(lapply(samp.trees, function(x) sum(x$edge.length)))
  
  # Calculate GMax for every input tree:
  input.permutations[, "GMax"] <- unlist(lapply(trees, function(x) sum(x$root.time - ages[, "FAD"])))
  
  # Calculate GMax for every random tree:
  rand.permutations[, "GMax"] <- unlist(lapply(rand.trees, function(x) sum(x$root.time - x$ages[, "FAD"])))
  
  # If sampling trees calculate GMax for every sampled tree:
  if(SamplingTrees) samp.permutations[, "GMax"] <- unlist(lapply(samp.trees, function(x) sum(x$root.time - x$ages[, "FAD"])))
  
  # Calculate GMin for every input tree:
  input.permutations[, "GMin"] <- unlist(lapply(trees, function(x) (x$root.time - max(ages[, "FAD"])) + diff(range(ages[, "FAD"]))))
  
  # Calculate GMin for every random tree:
  rand.permutations[, "GMin"] <- unlist(lapply(rand.trees, function(x) (x$root.time - max(x$ages[, "FAD"])) + diff(range(x$ages[, "FAD"]))))
  
  # If sampling trees calculate GMin for every sampled tree:
  if(SamplingTrees) samp.permutations[, "GMin"] <- unlist(lapply(samp.trees, function(x) (x$root.time - max(x$ages[, "FAD"])) + diff(range(x$ages[, "FAD"]))))
  
  # Calculate Relative Completeness Index for every input tree:
  input.permutations[, "RCI"] <- (1 - (input.permutations[, "MIG"] / input.permutations[, "SRL"])) * 100
  
  # Calculate Relative Completeness Index for every random tree:
  rand.permutations[, "RCI"] <- (1 - (rand.permutations[, "MIG"] / rand.permutations[, "SRL"])) * 100
  
  # If sampling trees calculate Relative Completeness Index for every sampled tree:
  if(SamplingTrees) samp.permutations[, "RCI"] <- (1 - (samp.permutations[, "MIG"] / samp.permutations[, "SRL"])) * 100
  
  # Calculate Gap Excess Ratio for every input tree:
  input.permutations[, "GER"] <- 1 - ((input.permutations[, "MIG"] - input.permutations[, "GMin"]) / (input.permutations[, "GMax"] - input.permutations[, "GMin"]))
  
  # Calculate Gap Excess Ratio for every random tree:
  rand.permutations[, "GER"] <- 1 - ((rand.permutations[, "MIG"] - rand.permutations[, "GMin"]) / (rand.permutations[, "GMax"] - rand.permutations[, "GMin"]))
  
  # If sampling trees calculate Gap Excess Ratio for every sampled tree:
  if(SamplingTrees) samp.permutations[, "GER"] <- 1 - ((samp.permutations[, "MIG"] - samp.permutations[, "GMin"]) / (samp.permutations[, "GMax"] - samp.permutations[, "GMin"]))
  
  # Calculate Manhattan Stratigraphic Measure* for every input tree:
  input.permutations[, "MSM*"] <- input.permutations[, "GMin"] / input.permutations[, "MIG"]
  
  # Calculate Manhattan Stratigraphic Measure* for every random tree:
  rand.permutations[, "MSM*"] <- rand.permutations[, "GMin"] / rand.permutations[, "MIG"]
  
  # If sampling trees calculate Manhattan Stratigraphic Measure* for every sampled tree:
  if(SamplingTrees) samp.permutations[, "MSM*"] <- samp.permutations[, "GMin"] / samp.permutations[, "MIG"]
  
  # Calculate Gap Excess Ratio* for every input tree:
  input.permutations[, "GER*"] <- unlist(lapply(as.list(input.permutations[, "MIG"]), function(x) sum(x <= rand.permutations[, "MIG"]) / rand.perm))
  
  # If sampling trees calculate Gap Excess Ratio* for every sampled tree:
  if(SamplingTrees) samp.permutations[, "GER*"] <- unlist(lapply(as.list(samp.permutations[, "MIG"]), function(x) sum(x <= rand.permutations[, "MIG"]) / rand.perm))
  
  # Calculate Gap Excess Ratiot for every input tree:
  input.permutations[, "GERt"] <- unlist(lapply(as.list(input.permutations[, "MIG"]), function(x) {y <- 1 - ((x - min(rand.permutations[, "MIG"])) / (max(rand.permutations[, "MIG"]) - min(rand.permutations[, "MIG"]))); if(y > 1) y <- 1; if(y < 0) y <- 0; y}))
  
  # If sampling trees calculate Gap Excess Ratiot for every sampled tree:
  if(SamplingTrees) samp.permutations[, "GERt"] <- unlist(lapply(as.list(samp.permutations[, "MIG"]), function(x) {y <- 1 - ((x - min(rand.permutations[, "MIG"])) / (max(rand.permutations[, "MIG"]) - min(rand.permutations[, "MIG"]))); if(y > 1) y <- 1; if(y < 0) y <- 0; y}))
  
  # Estimate p-value for SCI:
  input.permutations[, "est.p.SCI"] <- stats::pnorm(asin(sqrt(input.permutations[, "SCI"])), mean(asin(sqrt(rand.permutations[, "SCI"]))), stats::sd(asin(sqrt(rand.permutations[, "SCI"]))), lower.tail = FALSE)
  
  # Estimate p-value for RCI:
  input.permutations[, "est.p.RCI"] <- stats::pnorm(input.permutations[, "RCI"], mean(rand.permutations[, "RCI"]), stats::sd(rand.permutations[, "RCI"]), lower.tail = FALSE)
  
  # Estimate p-value for GER:
  input.permutations[, "est.p.GER"] <- stats::pnorm(asin(sqrt(input.permutations[, "GER"])), mean(asin(sqrt(rand.permutations[, "GER"]))), stats::sd(asin(sqrt(rand.permutations[, "GER"]))), lower.tail = FALSE)
  
  # Estimate p-value for MSM*:
  input.permutations[, "est.p.MSM*"] <- stats::pnorm(asin(sqrt(input.permutations[, "MSM*"])), mean(asin(sqrt(rand.permutations[, "MSM*"]))), stats::sd(asin(sqrt(rand.permutations[, "MSM*"]))), lower.tail = FALSE)
  
  # If sampling trees:
  if(SamplingTrees) {
    
    # Estimate p-value for SCI:
    samp.permutations[, "est.p.SCI"] <- stats::pnorm(asin(sqrt(samp.permutations[, "SCI"])), mean(asin(sqrt(rand.permutations[, "SCI"]))), stats::sd(asin(sqrt(rand.permutations[, "SCI"]))), lower.tail = FALSE)
    
    # Estimate p-value for RCI:
    samp.permutations[, "est.p.RCI"] <- stats::pnorm(samp.permutations[, "RCI"], mean(rand.permutations[, "RCI"]), stats::sd(rand.permutations[, "RCI"]), lower.tail = FALSE)
    
    # Estimate p-value for GER:
    samp.permutations[, "est.p.GER"] <- stats::pnorm(asin(sqrt(samp.permutations[, "GER"])), mean(asin(sqrt(rand.permutations[, "GER"]))), stats::sd(asin(sqrt(rand.permutations[, "GER"]))), lower.tail = FALSE)
    
    # Estimate p-value for MSM*:
    samp.permutations[, "est.p.MSM*"] <- stats::pnorm(asin(sqrt(samp.permutations[, "MSM*"])), mean(asin(sqrt(rand.permutations[, "MSM*"]))), stats::sd(asin(sqrt(rand.permutations[, "MSM*"]))), lower.tail = FALSE)
    
  }
  
  # Calculate Wills p-value for MIG:
  input.permutations[, "p.Wills"] <- unlist(lapply(as.list(input.permutations[, "MIG"]), function(x) 1 - (sum(x < rand.permutations[, "MIG"]) / rand.perm)))
  
  # If sampling trees calculate Wills p-value for MIG:
  if(SamplingTrees) samp.permutations[, "p.Wills"] <- unlist(lapply(as.list(samp.permutations[, "MIG"]), function(x) 1 - (sum(x < rand.permutations[, "MIG"]) / rand.perm)))
  
  # If sampled trees exist set class of each tree variable as multiPhylo:
  if(samp.perm > 0) class(rand.trees) <- class(samp.trees) <- class(trees) <- "multiPhylo"
  
  # If sampled trees do not exist set class of each tree variable as multiPhylo:
  if(samp.perm == 0) class(rand.trees) <- class(trees) <- "multiPhylo"

  # Compile output as list:
  output <- list(input.permutations, samp.permutations, rand.permutations, trees, samp.trees, rand.trees)

  # Add names to output:
  names(output) <- c("input.tree.results", "samp.permutation.results", "rand.permutations", "input.trees", "samp.trees", "rand.trees")

  # Return output:
  return(invisible(output))

}

# Test data run:
#x <- StratPhyloCongruence(trees = lapply(Asaphidae$tree, function(x) {x$edge.length <- runif(nrow(x$edge)); x}), ages = Asaphidae$ages, rlen = 0, method = "basic", samp.perm = 100, rand.perm = 100, hard = TRUE, randomly.sample.ages = TRUE, fix.topology = TRUE, fix.outgroup = TRUE, outgroup.taxon = "Dikelocephalus")

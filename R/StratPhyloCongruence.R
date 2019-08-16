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
#' @param randomly.sample.ages Whether to treat FAD and LAD as a range (randomly.sample.ages=FALSE) or an uncertainty (randomly.sample.ages=TRUE). If the latter then two ages are randomly sampled from the range and these are used as the FAD and LAD.
#' @param fix.topology Whether to allow tree shape to be random (fix.topology=FALSE) or to reflect the tree shape of the input tree(s) (fix.topology=TRUE).
#' @param fix.outgroup Whether to force the randomly generated trees to share the outgroup of the first input tree (fix.outgroup=TRUE) or not (fix.outgroup=FALSE)
#'
#' @details
#'
#' Cladograms of fossil taxa make explicit predictions about the successive appearance of taxa in the fossil record that can be compared with their observed stratigraphic ranges. Several methods have been developed to quantify this "fit to stratigraphy" of phylogenetic hypotheses, and these can be assessed in both a statistic (measuring the apparent strength of this congruence) and an associated significance test (p-value) based on generating random topologies for the same taxon set.
#'
#' This function produces both values for all four main metrics: the Stratigraphic Consistency Index (SCI; Huelsenbeck 1994), the Relative Consistency Index (RCI; Benton and Storrs 1994), the Manhattan Stratigraphic Measure (MSM*; Siddall 1998; Pol and Norell 2001), and the Gap Excess Ratio (GER; Wills 1999).
#'
#' SCI - Stratigraphic Consistency Index
#'
#' The SCI works by assessing the "consistency" of nodes. A node is considered stratigraphically consistent if its oldest descendant is the same age or younger than the oldest descendant of the preceding node. The SCI is thus given simply as:
#'
#' SCI = C/N
#'
#' Where C is the sum of all the consistent nodes and N is the total number of nodes - 1. (As there is no node preceding the root there is no basis on which to estimate its consistency.) This value can range from zero (maximally inconsistent) to one (maximally consistent). However, a potential criticism of the SCI is that a high value may be returned when in fact a single inconsistent node may represent a very large amount of missing history (measured in unsampled units or millions of years), whereas a low SCI may represent relatively few unsampled units or millions of years.
#'
#' RCI - Relative Completeness Index
#'
#' The RCI was the first method to explicitly account for the absolute amount of missing data implied by the tree. This figure is usually expressed as the Minimum Implied Gap (MIG), a term also used by both the MSM and GER (see below), and corresponds to the sum of the branch lengths excluding the duration of the terminals (the observed ranges of the taxa). The RCI expresses the MIG as a proportion of the sum of the observed ranges (Simple Range Length; SRL) of the taxa converted to a percentage:
#'
#' RCI = (1 − (MIG/SRL)) ∗ 100percent
#'
#' Importantly this value is not confined to a 0 to 100 percent scale, and can have both negative values and values greater than 100 percent, which can make it difficult to interpret.
#'
#' MSM - Manhattan Stratigraphic Measure
#'
#' The MSM was the first method to account for both the absolute MIG and range on a confined zero to one scale. It is expressed as:
#'
#' MSM = Lm/L0
#'
#' Where L0 is the length of the tree expressed by optimising times of first appearance on to the tree as a Sankoff character and taking the total length. Lm represents the same process, but for the optimal possible tree given the same set of first appearances. However, Pol and Norell (2001) noted a critical flaw in this approach, specifically that the Sankoff optimisation is reversible, meaning that nodes in the topology are allowed to be younger than their descendants, leading in some cases to a poor fit to stratigraphy being perceived as a good fit. Instead they suggest modifying the character step matrix to make the cost of reversals effectively infinite and hence impossible. Thus the values for L0 and Lm are modified accordingly. This approach they termed MSM* and is the implementation of MSM used here. This statistic can be expressed as:
#'
#' MSM∗ = Gmin/M IG
#'
#'Where Gmin represents the MIG for the tree with the optimal fit to stratigraphy. In effect this is a completely unbalanced tree where the youngest pair of taxa are the most deeply nested and successive outgroups represent the next oldest taxon. Theoretically MSM* ranges from one (the best fit as the observed tree is the maximally consistent tree) to zero (the least optimal tree). However, in effect no tree can have a value of zero as its MIG would have to be equal to infinity.
#'
#' GER - Gap Excess Ratio
#'
#' The GER represents a method that accounts for MIG, ranges from zero to one, and the best and worst fits to stratigraphy are both practically realisable. It can be expressed as:
#'
#' GER = 1 − ((MIG − Gmin)/(Gmax − Gmin))
#'
#' Where Gmax represents the MIG of the tree with the worst possible fit to stratigraphy. This is in effect any topology where the oldest taxon is the most deeply nested such that every clade in the tree contains it and hence must be minimally that old.
#'
#' P-values
#'
#' In isolation all four methods suffer from an inability to reject the null hypothesis that an apparent good fit to stratigraphy may be generated by chance alone. In practice this can be tested by generating a set of random topologies, calculating the fit to stratigraphy measure, and then either fitting a normal distribution to the resulting values (to get an estimated p-value) or assessing the relative position of the MIG of the observed (and sampled) tree(s) to get an absolute p-value. (Note that the assumption of normality may not always hold and the former approach should be used at the user’s discretion. However, it should be noted that for the SCI, MSM*, and GER p-values are calculated after first transforming the data by taking the arcsine of the square root of each value.) The reason for having two sets of p-values is that if the observed trees fall completely outside the range of the random topologies they will be given an extreme p-value (0 or 1) that may be misleading. In such cases the estimated value may be more accurate.
#'
#' P-values should be interpreted as the probability of the null: that the observed tree(s) have an equal or worse fit to stratigraphy than the sample of random trees. Thus if the p-values are very small the user can reject the null hypothesis in favour of the alternative: that the observed tree(s) have a better fit to stratigraphy than expected by chance alone.
#'
#' Modifications of the GER
#'
#' More recently Wills et al. (2008) introduced two new versions of the GER that take advantage of the distribution of MIGs from the set of randomly generated topologies. The first of these (GERt)uses the extreme values of the random topologies as modified versions of Gmax and Gmin, termed Gtmax and Gtmin respectively. GERt is thus expressed as:
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
#' Polytomies and age uncertainties
#'
#' Alongside the input trees the user can also create an additional set of sampled trees based on the input trees. This option is automatically implemented when choosing either hard=FALSE or randomly.sample.ages=TRUE, and the total number of permutations to perform dictated by samp.perm=N. This process works by first sampling from the set of input tree(s) and then randomly resolving any polytomies (if hard=FALSE) to ensure all sampled trees are fully dichotomous. (At present the function does not allow the various options laid out in Boyd et al. 2011, but the user can achieve this effect by modifying the input trees themselves.) Then if randomly.sample.ages=TRUE the FAD and LAD are treated as bounds of a uniform distribution which is sampled at random. This allows the user to get results for a set of trees that account for uncertainty in dating (as outlined in Pol and Norell 2006). (Note that two dates are picked for each taxon to avoid the problem of having an SRL of zero that would cause a divide by zero error for the RCI metric.)
#'
#' All fit to stratigraphy measures calculated for the input trees are then repeated for the sampled trees. However, if both hard=TRUE and randomly.sample.ages=FALSE (the defaults) no set of sampled trees will be created.
#'
#' In all cases when using the function users will see a progress bar that indicates the general progress through the combined set of trees (input, sampled, and randomly generated). This serves as a useful indicator of the time it will take for the function to finish. Here default values for samp.perm and rand.perm are both set at 1000, but the user may wish to lower these (to decrease calculation time) or increase them (to enhance accuracy).
#'
#' Additional options
#'
#' Note that because this function uses DatePhylo the user has the option of using different tree dating algorithms than the basic method (equivalent to the basic method in the paleotree package) employed in all the published studies cited above (and the default option here). The dating method used will apply to all trees generated, including the input, sampled, and randomly generated topologies. In all cases the time-scaled trees are returned with the function output.
#'
#' A final option (fix.outgroup=TRUE) allows the user to always use the same outgroup taxon (based on the first input tree) for all randomly generated topologies. Because the outgroup will often be the oldest taxon and its position in the input topologies is not allowed to vary letting it do so in the random topologies may lead to inferring a better fit to stratigraphy for the observed tree(s) than is fair. Fixing the outgroup thus ameliorates this potential bias and is the default option here.
#'
#' @return
#'
#' \item{Output}{Explanantion.}
#' \item{input.tree.results}{A matrix with a row for each input tree and columns indicating the values for SCI, RCI, GER and MSM* and their estimated probabilities assuming a normal distribution (est.p.SCI, est.p.RCI, est.p.GER, and est.p.MSM*) as well as GERt, GER*, MIG, and p.Wills (their probability as position within the MIGs of the random topologies).}
#' \item{samp.permutation.results}{If used, a matrix with a row for each sampled tree (up to samp.perm) and columns indicating the values for SCI, RCI, GER and MSM* and their estimated probabilities assuming a normal distribution (est.p.SCI, est.p.RCI, est.p.GER, and est.p.MSM*) as well as GERt, GER*, MIG, and p.Wills (their probability as position within the MIGs of random topologies).}
#' \item{rand.permutations}{A matrix with a row for each randomly generated tree (up to rand.perm) and columns indicating the values for SCI, RCI, GER, MSM*, and MIG.}
#' \item{input.trees}{The input tree(s) as a phylo or multiphylo object, with branches scaled to time according to the input values passed to DatePhylo.}
#' \item{samp.trees}{The sampled tree(s) as a phylo or multiphylo object, with branches scaled to time according to the input values passed to DatePhylo.}
#' \item{rand.trees}{The randomly generated tree(s) as a phylo or multiphylo object, with branches scaled to time according to the input values passed to DatePhylo.}
#'
#' @author Mark A. Bell \email{mark.bell521@@gmail.com} and Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Bell, M. A. and Lloyd, G. T., 2015. strap: an R package for plotting phylogenies against stratigraphy and assessing their stratigraphic congruence. Palaeontology, 58, 379-389.
#' Benton, M. J. and Storrs, G. W., 1994. Testing the quality of the fossil record: palaeontological knowledge is improving. Geology, 22, 111-114.
#' Boyd, C. A., Cleland, T. P., Marrero, N. L. and Clarke, J. A., 2011. Exploring the effects of phylogenetic uncertainty and consensus trees on stratigraphic consistency scores: a new program and a standardized method. Cladistics, 27, 52-60.
#' Huelsenbeck, J. P., 1994. Comparing the stratigraphic record to estimates of phylogeny. Paleobiology, 20, 470-483.
#' Pol, D. and Norell, M. A., 2001. Comments on the Manhattan Stratigraphic Measure. Cladistics, 17, 285-289.
#' Pol, D. and Norell, M. A., 2006. Uncertainty in the age of fossils and the stratigraphic fit to phylogenies. Systematic Biology, 55, 512-521.
#' Siddall, M. E., 1998. Stratigraphic fit to phylogenies: a proposed solution. Cladistics, 14, 201-208.
#' Wills, M. A., 1999. Congruence between phylogeny and stratigraphy: randomization tests and the Gap Excess Ratio. Systematic Biology, 48, 559-580.
#' Wills, M. A., Barrett, P. M. and Heathcote, J. F., 2008. The modified Gap Excess Ratio (GER*) and the stratigraphic congruence of dinosaur phylogenies. Systematic Biology, 57, 891-904.
#'
#' @keywords congruence,phylogeny,stratigraphy
#'
#' @examples
#'
#' # Calculate stratigraphic fit measures treating ages as ranges
#' # (permutation numbers used are lower than recommended for standard use):
#' fit.to.strat.1 <- StratPhyloCongruence(trees = Dipnoi$tree, ages = Dipnoi$ages, rlen = 0,
#'   method = "basic", samp.perm = 5, rand.perm = 5, hard = TRUE,
#'   randomly.sample.ages = FALSE, fix.topology = TRUE, fix.outgroup = TRUE)
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
#'   randomly.sample.ages = TRUE, fix.topology = TRUE, fix.outgroup = TRUE)
#'
#' @export StratPhyloCongruence
StratPhyloCongruence <- function(trees, ages, rlen = 0, method = "basic", samp.perm = 1000, rand.perm = 1000, hard = TRUE, randomly.sample.ages = FALSE, fix.topology = TRUE, fix.outgroup = TRUE) {

  # Make sure permutation numbers are not negative:
  if(samp.perm <= 0 || rand.perm <= 0) stop("Number of permutations must be positive.")

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
  
  # Create separate ages matrices for sampled and random permutations:
  samp.ages <- rand.ages <- ages

  # Set as null variables that may not be used later but are still included in the output:
  samp.permutations <- samp.trees <- NULL
  
  # Set root age:
  root.age <- max(ages[, "FAD"]) + rlen
  
  # If a single tree:
  if(class(trees) == "phylo") {
    
    # Convert to list to ensure a standard format for input tree(s):
    trees <- list(trees)
    
    # Set class as multiPhylo:
    class(trees) <- "multiPhylo"

  }

  # Get outgroup taxon for later re-rooting:
  if(fix.outgroup) {

    # Find maximum node height (indicating rooting taxon):
    max.node.height <- max(node.height(trees[[1]])[1:Ntip(trees[[1]])])

    # Find rows with maximum node height:
    max.node.height.rows <- grep(TRUE, node.height(trees[[1]])[1:Ntip(trees[[1]])] == max.node.height)
    
    # If there is not a single outgroup taxon for rooting:
    if(length(max.node.height.rows) > 1) {

      # Reset fix.outgroup to FALSE:
      fix.outgroup <-  FALSE

      # Warn user:
      print("Tree doesn't have single outgroup, fix.outgroup set to FALSE.")

    } else {

      # Set outgroup taxon:
      outgroup.taxon <- trees[[1]]$tip.label[max.node.height.rows]

    }
  
  }

  # Get total number of input trees:
  Ntrees <- length(trees)
  
  # If trees have edge lengths and are using Ruta dating method:
  if(length(grep("edge.length", names(unlist(trees)))) > 0 && method == "ruta") {
  
    # Collect all edge lengths in a single vector:
    edge.lengths <- as.numeric(unlist(trees)[grep("edge.length", names(unlist(trees)))])
  
  # If trees do not have edge lengths or if not using Ruta method:
  } else {
    
    # Set edge lengths vector to NULL:
    edge.lengths <- NULL

  }
  
  # If user wants to fix the topology as in GERt:
  if(fix.topology) {
    
    # Create random trees by resampling input trees:
    rand.trees <- trees[sample(1:length(trees), rand.perm, replace = TRUE)]
    
    # Reset class as multiPhylo:
    class(rand.trees) <- "multiPhylo"
    
    # For each random tree:
    for(i in 1:length(rand.trees)) {
      
      # Make sure it is fully bifurcating:
      rand.trees[[i]] <- multi2di(rand.trees[[i]])
      
      # If fixing the outgroup taxon (recommended):
      if(fix.outgroup) {

        # Randomly sample taxon names:
        sampled.names <- sample(rand.trees[[i]]$tip.label)

        # Ensure outgroup taxon is listed first:
        sampled.names <- c(outgroup.taxon, sampled.names[grep(FALSE, sampled.names == outgroup.taxon)])

        # Overwrite tip names with sampeld names:
        rand.trees[[i]]$tip.label <- sampled.names

      # If not fixing the outgroup taxon:
      } else {
  
        # Shuffle all tip labels:
        rand.trees[[i]]$tip.label <- sample(rand.trees[[i]]$tip.label)

      }
  
    }
    
  # If user does not want to fix the topology:
  } else {
  
    # Make random trees for permutations:
    rand.trees <- rmtree(rand.perm, Ntip(trees[[1]]))
  
  }

  # If edge.lengths is not NULL:
  if(!is.null(edge.lengths)) {
    
    # For each random tree:
    for(i in 1:rand.perm) {
      
      # Draw edge lengths from a uniform distribution between the limits of observed edge lengths:
      rand.trees[[i]]$edge.length <- runif(length(rand.trees[[i]]$edge.length), min = min(edge.lengths), max = max(edge.lengths))
    
    }

  }

  # For each tree in trees:
  for(i in 1:length(trees)) {
  
    # Date tree using settings specified in function:
    trees[[i]] <- DatePhylo(tree=trees[[i]], ages = ages, rlen = rlen, method = method, add.terminal = FALSE)
  
  }
  
  # If using the sample permutation loop:
  if(hard == FALSE || randomly.sample.ages == TRUE) {

    # Randomly sample numbers from 1 to Ntrees:
    tree.nos <- sample(c(1:Ntrees), samp.perm, replace = TRUE)
    
    # Create sampled trees list:
    samp.trees <- rmtree(samp.perm, Ntip(trees[[1]]))
    
    # For each sample permutation:
    for(i in 1:samp.perm) {
    
      # If polytomies are to be considered hard:
      if(hard == TRUE) {
        
        # Store randomly sampled tree:
        samp.trees[[i]] <- trees[[tree.nos[i]]]
      
      # If polytomies are to be considered soft:
      } else {
        
        # Store randomly sampled randomly bifurcated tree:
        samp.trees[[i]] <- multi2di(trees[[tree.nos[i]]])

      }
      
    }
    
  # If not permuting ages of polytomie resolutions:
  } else {
  
    # Reset samp.perm to zero:
    samp.perm <- 0

  }

  # Additional steps if randomly sampling ages:
  if(randomly.sample.ages) {
    
    # Create matrix to store randomly sampled ages ages:
    randomly.sampled.ages <- matrix(nrow = nrow(ages), ncol = max(c(samp.perm, rand.perm)) * 2)
    
    # For each taxon:
    for(i in 1:nrow(ages)) {
    
      # Draw random ages from a uniform distribution between FAD and LAD:
      randomly.sampled.ages[i, ] <- runif(samp.perm * 2, min = ages[i, "LAD"], ages[i, "FAD"])

    }
    
    # Take first half of matrix:
    firsthalf <- as.vector(randomly.sampled.ages[, 1:(length(randomly.sampled.ages[1, ]) / 2)])
    
    # Take second half of matrix:
    secondhalf <- as.vector(randomly.sampled.ages[, ((length(randomly.sampled.ages[1, ]) / 2) + 1):length(randomly.sampled.ages[1, ])])
    
    # Sort halves relative to each other to get FADs and LADs:
    randomly.sampled.ages <- apply(cbind(firsthalf, secondhalf), 1, sort)
    
    # Define randomly sampled ages FADs:
    randomly.sampled.ages.FAD <- matrix(randomly.sampled.ages[2, ], nrow = nrow(ages))
    
    # Define randomly sampled ages LADs:
    randomly.sampled.ages.LAD <- matrix(randomly.sampled.ages[1, ], nrow = nrow(ages))

    # Add rownames:
    rownames(randomly.sampled.ages.FAD) <- rownames(randomly.sampled.ages.LAD) <- rownames(ages)
    
    # Set randomly sampled ages root lengths:
    randomly.sampled.ages.rlen <- rlen + max(ages) - apply(randomly.sampled.ages.FAD, 2, max)
    
    # Calculate Gmax (the worst possible fit to stratigraphy):
    randomly.sampled.ages.Gmax <- apply(root.age - randomly.sampled.ages.FAD, 2, sum)
    
    # Calculate Gmin (the best possible fit to stratigraphy - if no ADRs):
    randomly.sampled.ages.Gmin <- (apply(apply(apply(randomly.sampled.ages.FAD, 2, sort), 2, diff), 2, sum)) + (2 * (root.age - apply(randomly.sampled.ages.FAD, 2, max)))
    
    # Calculate sum of stratigraphic ranges of taxa (SRL):
    randomly.sampled.ages.SRL <- apply(randomly.sampled.ages.FAD - randomly.sampled.ages.LAD, 2, sum)
    
  }
  
  # Create matrix to store output with extra column for GERt:
  input.permutations <- matrix(ncol = 10, nrow = length(trees))
  
  # Annotate matrix:
  colnames(input.permutations) <- c("SCI", "RCI", "GER", "MSM*", "est.p.SCI", "est.p.RCI", "est.p.GER", "est.p.MSM*", "GER*", "GERt")
  rownames(input.permutations) <- paste("tree_", c(1:length(trees)),sep = "")

  # Calculate Gmax (the worst possible fit to stratigraphy):
  Gmax <- sum(root.age - ages[, "FAD"])
  
  # Calculate Gmin (the best possible fit to stratigraphy - if no ADRs):
  Gmin <- diff(range(ages[, "FAD"])) + (2 * (root.age - max(ages[, "FAD"])))

  # Start progress bar:
  pb <- utils::txtProgressBar(min = 0, max = (Ntrees + samp.perm + rand.perm), style = 3)

  # Text bar counter:
  counter <- 0
  
  # Vector to store random MIGs for calculating GER* and GERt:
  rand.MIG <- samp.MIG <- obs.MIG <- vector(mode = "numeric")

  # Cycle through each tree in the list:
  for (i in 1:length(trees)) {
  
    # Isolate ith tree:
    tree <- trees[[i]]
    
    # MIG calculation:
    obs.MIG[i] <- MIG <- sum(tree$edge.length)
    
    # Create node matrix for SCI calculation:
    newnodes <- matrix(ncol = 2, nrow = Nnode(tree) - 1)
    
    # Add rownames to matrix:
    rownames(newnodes) <- c((Ntip(tree) + 2):(Ntip(tree) + Nnode(tree)))
    
    # Add column names to matrix:
    colnames(newnodes) <- c("anc", "dec")
    
    # For each node (excluding the root):
    for (j in 1:nrow(newnodes)) {
      
      # Get descendants for current node:
      dec <- tree$tip.label[FindDescendants(as.numeric(rownames(newnodes)[j]), tree)]
      
      # Get descendants for immediate ancestral node:
      anc <- tree$tip.label[FindDescendants(tree$edge[match(as.numeric(rownames(newnodes)[j]), tree$edge[, 2]), 1], tree)]
      
      # Remove descendnats from ancestral node already found in current node:
      anc <- anc[match(anc, dec, nomatch = 0) == 0]
      
      # Store oldest age for current node descendants:
      newnodes[j, 2] <- max(ages[dec, "FAD"])
      
      # Store oldest age for ancestral node descendants:
      newnodes[j, 1] <- max(ages[anc, "FAD"])
      
    }
    
    # SCI calculation:
    input.permutations[i, "SCI"] <- sum(newnodes[, "anc"] >= newnodes[, "dec"]) / nrow(newnodes)
    
    # RCI calculation:
    input.permutations[i, "RCI"] <- (1 - (MIG / SRL)) * 100
    
    # GER calculation:
    input.permutations[i, "GER"] <- 1 - ((MIG - Gmin) / (Gmax - Gmin))
    
    # MSM* calculation:
    input.permutations[i, "MSM*"] <- Gmin / MIG
    
    # Update counter:
    counter <- counter + 1
    
    # Update progress bar:
    utils::setTxtProgressBar(pb, counter)
    
  }
  
  # If polytomies are not hard (i.e. if they are soft):
  if(hard == FALSE || randomly.sample.ages == TRUE) {
    
    # Create matrix to store permutation results with extra column for GERt:
    samp.permutations <- matrix(nrow = samp.perm, ncol = 10)

    # Add column names:
    colnames(samp.permutations) <- c("SCI", "RCI", "GER", "MSM*", "est.p.SCI", "est.p.RCI", "est.p.GER", "est.p.MSM*", "GER*", "GERt")

    # For each sample permutation:
    for(i in 1:samp.perm) {
      
      # Isolate ith sample tree:
      samp.tree <- samp.trees[[i]]
      
      # If using randomly sampled ages:
      if(randomly.sample.ages) {
        
        # Overwrite rand.ages with ith randomly.sampled.ages
        samp.ages[, "FAD"] <- randomly.sampled.ages.FAD[rownames(rand.ages), i]
        
        # Overwrite rand.ages with ith randomly.sampled.ages
        samp.ages[, "LAD"] <- randomly.sampled.ages.LAD[rownames(rand.ages), i]
        
        # Overwrite rlen with ith randomly sampled ages rlen:
        rlen <- randomly.sampled.ages.rlen[i]
        
        # Overwrite Gmin with randomly sampled ages Gmin:
        Gmin <- randomly.sampled.ages.Gmin[i]
        
        # Overwrite Gmax with randomly sampled ages Gmax:
        Gmax <- randomly.sampled.ages.Gmax[i]
        
        # Overwrite SRL with randomly sampled ages SRL:
        SRL <- randomly.sampled.ages.SRL[i]
        
      }
      
      # Date tree using settings specified in function:
      samp.tree <- DatePhylo(tree = samp.tree, ages = samp.ages, rlen = rlen, method = method, add.terminal = FALSE)
      
      # MIG calculation:
      samp.MIG[i] <- MIG <- sum(samp.tree$edge.length)
      
      # Create node matrix for SCI calculation:
      newnodes <- matrix(ncol = 2, nrow = Nnode(samp.tree) - 1)
      
      # Add rownames to matrix:
      rownames(newnodes) <- c((Ntip(samp.tree) + 2):(Ntip(samp.tree) + Nnode(samp.tree)))
      
      # Add column names to matrix:
      colnames(newnodes) <- c("anc", "dec")
      
      # For each node (excluding the root):
      for (j in 1:nrow(newnodes)) {
        
        # Get descendants for current node:
        dec <- samp.tree$tip.label[FindDescendants(as.numeric(rownames(newnodes)[j]), samp.tree)]
        
        # Get descendants for immediate ancestral node:
        anc <- samp.tree$tip.label[FindDescendants(samp.tree$edge[match(as.numeric(rownames(newnodes)[j]), samp.tree$edge[, 2]), 1], samp.tree)]
        
        # Remove descendnats from ancestral node already found in current node:
        anc <- anc[match(anc, dec, nomatch = 0) == 0]
        
        # Store oldest age for current node descendants:
        newnodes[j, 2] <- max(samp.ages[dec, "FAD"])
        
        # Store oldest age for ancestral node descendants:
        newnodes[j, 1] <- max(samp.ages[anc, "FAD"])
        
      }
      
      # SCI calculation:
      samp.permutations[i, "SCI"] <- sum(newnodes[, "anc"] >= newnodes[, "dec"]) / nrow(newnodes)
      
      # RCI calculation:
      samp.permutations[i, "RCI"] <- (1 - (MIG / SRL)) * 100
      
      # GER calculation:
      samp.permutations[i, "GER"] <- 1 - ((MIG - Gmin) / (Gmax - Gmin))
      
      # MSM* calculation:
      samp.permutations[i, "MSM*"] <- Gmin / MIG

      # Save sampled tree:
      samp.trees[[i]] <- samp.tree
      
      # Update counter:
      counter <- counter + 1
      
      # Update progress bar:
      utils::setTxtProgressBar(pb, counter)
      
    }

    # Add MIG to the sample permutaion results:
    samp.permutations <- cbind(samp.permutations, samp.MIG)

    # Add column heading for MIG:
    colnames(samp.permutations)[length(colnames(samp.permutations))] <- "MIG"

  }

  # Create matrix to store permutation results:
  rand.permutations <- matrix(nrow = rand.perm, ncol = 4)
  
  # Add column names:
  colnames(rand.permutations) <- c("SCI", "RCI", "GER", "MSM*")
  
  # For each permutation:
  for (i in 1:rand.perm) {
    
    # Isolate ith random tree:
    rand.tree <- rand.trees[[i]]

    # If fixing the outgroup (recommended):
    if(fix.outgroup == TRUE && fix.topology == FALSE) {

      # Overwrite random names with real taxon names
      rand.tree$tip.label <- sample(trees[[1]]$tip.label)

      # Re-root tree on outgroup taxon:
      rand.tree <- root(rand.tree, outgroup.taxon, resolve.root = TRUE)
    
    }

    # If not fixing the outgroup:
    if(fix.outgroup == FALSE) {

      # Overwrite random names with real taxon names
      rand.tree$tip.label <- sample(trees[[1]]$tip.label)

    }

    # If using randomly sampled ages:
    if(randomly.sample.ages) {
      
      # Overwrite rand.ages with ith randomly.sampled.ages
      rand.ages[, "FAD"] <- randomly.sampled.ages.FAD[rownames(rand.ages), i]
      
      # Overwrite rand.ages with ith randomly.sampled.ages
      rand.ages[, "LAD"] <- randomly.sampled.ages.LAD[rownames(rand.ages), i]

      # Overwrite rlen with ith randomly sampled ages rlen:
      rlen <- randomly.sampled.ages.rlen[i]
      
      # Overwrite Gmin with randomly sampled ages Gmin:
      Gmin <- randomly.sampled.ages.Gmin[i]

      # Overwrite Gmax with randomly sampled ages Gmax:
      Gmax <- randomly.sampled.ages.Gmax[i]
      
      # Overwrite SRL with randomly sampled ages SRL:
      SRL <- randomly.sampled.ages.SRL[i]
      
    }
    
    # Date tree using settings specified in function:
    rand.tree <- DatePhylo(tree = rand.tree, ages = rand.ages, rlen = rlen, method = method, add.terminal = FALSE)
    
    # MIG calculation:
    rand.MIG[i] <- MIG <- sum(rand.tree$edge.length)
    
    # Create node matrix for SCI calculation:
    newnodes <- matrix(ncol = 2, nrow = Nnode(rand.tree) - 1)
    
    # Add rownames to matrix:
    rownames(newnodes) <- c((Ntip(rand.tree) + 2):(Ntip(rand.tree) + Nnode(rand.tree)))
    
    # Add column names to matrix:
    colnames(newnodes) <- c("anc", "dec")
    
    # For each node (excluding the root):
    for (j in 1:nrow(newnodes)) {
      
      # Get descendants for current node:
      dec <- rand.tree$tip.label[FindDescendants(as.numeric(rownames(newnodes)[j]), rand.tree)]
      
      # Get descendants for immediate ancestral node:
      anc <- rand.tree$tip.label[FindDescendants(rand.tree$edge[match(as.numeric(rownames(newnodes)[j]), rand.tree$edge[, 2]), 1], rand.tree)]
      
      # Remove descendnats from ancestral node already found in current node:
      anc <- anc[match(anc, dec, nomatch = 0) == 0]
      
      # Store oldest age for current node descendants:
      newnodes[j, 2] <- max(rand.ages[dec, "FAD"])
      
      # Store oldest age for ancestral node descendants:
      newnodes[j, 1] <- max(rand.ages[anc, "FAD"])
      
    }
    
    # SCI calculation:
    rand.permutations[i, "SCI"] <- sum(newnodes[, "anc"] >= newnodes[, "dec"]) / nrow(newnodes)
    
    # RCI calculation:
    rand.permutations[i, "RCI"] <- (1 - (MIG / SRL)) * 100
    
    # GER calculation:
    rand.permutations[i, "GER"] <- 1 - ((MIG - Gmin) / (Gmax - Gmin))
    
    # MSM* calculation:
    rand.permutations[i, "MSM*"] <- Gmin / MIG
    
    # Save randomisation tree:
    rand.trees[[i]] <- rand.tree
    
    # Update counter:
    counter <- counter + 1
    
    # Update progress bar:
    setTxtProgressBar(pb, counter)
  
  }
  
  # Close progress bar:
  close(pb)
  
  # Calculate input tree p.values for SCI:
  input.permutations[, "est.p.SCI"] <- stats::pnorm(asin(sqrt(input.permutations[, "SCI"])), mean(asin(sqrt(rand.permutations[, "SCI"]))), stats::sd(asin(sqrt(rand.permutations[, "SCI"]))), lower.tail = FALSE)

  # Calculate input tree p.values for RCI:
  input.permutations[, "est.p.RCI"] <- stats::pnorm(input.permutations[, "RCI"], mean(rand.permutations[, "RCI"]), stats::sd(rand.permutations[, "RCI"]), lower.tail = FALSE)
  
  # Calculate input tree p.values for GER:
  input.permutations[, "est.p.GER"] <- stats::pnorm(asin(sqrt(input.permutations[, "GER"])), mean(asin(sqrt(rand.permutations[, "GER"]))), stats::sd(asin(sqrt(rand.permutations[, "GER"]))), lower.tail = FALSE)
  
  # Calculate input tree p.values for MSM*:
  input.permutations[, "est.p.MSM*"] <- stats::pnorm(asin(sqrt(input.permutations[, "MSM*"])), mean(asin(sqrt(rand.permutations[, "MSM*"]))), stats::sd(asin(sqrt(rand.permutations[, "MSM*"]))), lower.tail = FALSE)
  
  # If sample permutations were performed:
  if(samp.perm > 0) {

    # Calculate sample permutation p.values for SCI:
    samp.permutations[, "est.p.SCI"] <- stats::pnorm(asin(sqrt(samp.permutations[, "SCI"])), mean(asin(sqrt(rand.permutations[, "SCI"]))), stats::sd(asin(sqrt(rand.permutations[, "SCI"]))), lower.tail = FALSE)
  
    # Calculate sample permutation p.values for RCI:
    samp.permutations[, "est.p.RCI"] <- stats::pnorm(samp.permutations[, "RCI"], mean(rand.permutations[, "RCI"]), stats::sd(rand.permutations[, "RCI"]), lower.tail = FALSE)
  
    # Calculate sample permutation p.values for GER:
    samp.permutations[, "est.p.GER"] <- stats::pnorm(asin(sqrt(samp.permutations[, "GER"])), mean(asin(sqrt(rand.permutations[, "GER"]))), stats::sd(asin(sqrt(rand.permutations[, "GER"]))), lower.tail = FALSE)
  
    # Calculate sample permutation p.values for MSM*:
    samp.permutations[, "est.p.MSM*"] <- stats::pnorm(asin(sqrt(samp.permutations[, "MSM*"])), mean(asin(sqrt(rand.permutations[, "MSM*"]))), stats::sd(asin(sqrt(rand.permutations[, "MSM*"]))), lower.tail = FALSE)
  
    # For each sampled tree:
    for(i in 1:samp.perm) {

      # Calculate GER* for ith sampled tree:
      samp.permutations[i, "GER*"] <- length(grep(TRUE, samp.MIG[i] <= rand.MIG)) / rand.perm

      # Calculate GERt:
      samp.permutations[i, "GERt"] <- 1 - ((samp.MIG[i] - min(rand.MIG)) / (max(rand.MIG) - min(rand.MIG)))
  
      # If GERt is negative re-scale to 0:
      if(samp.permutations[i, "GERt"] < 0) samp.permutations[i, "GERt"] <- 0
  
      # If GERt is greater than 1 re-scale to 1:
      if(samp.permutations[i, "GERt"] > 1) samp.permutations[i, "GERt"] <- 1

    }

  }

  # For each input tree:
  for (i in 1:Ntrees) {

    # Calculate GER* for ith input tree:
    input.permutations[i, "GER*"] <- length(grep(TRUE, obs.MIG[i] <= rand.MIG)) / rand.perm

    # Calculate GERt:
    input.permutations[i, "GERt"] <- 1 - ((obs.MIG[i] - min(rand.MIG)) / (max(rand.MIG) - min(rand.MIG)))

    # If GERt is negative re-scale to 0:
    if(input.permutations[i, "GERt"] < 0) input.permutations[i, "GERt"] <- 0

    # If GERt is greater than 1 re-scale to 1:
    if(input.permutations[i, "GERt"] > 1) input.permutations[i, "GERt"] <- 1
  
  }

  # Add MIG to input permutations:
  input.permutations <- cbind(input.permutations, obs.MIG)

  # Add MIG to random permutations:
  rand.permutations <- cbind(rand.permutations, rand.MIG)

  # Add MIG to column headings:
  colnames(input.permutations)[length(colnames(input.permutations))] <- "MIG"

  # Add MIG to column headings:
  colnames(rand.permutations)[length(colnames(rand.permutations))] <- "MIG"

  # If sample permutations were performed:
  if(samp.perm > 0) {

    # Add Wills p column:
    samp.permutations <- cbind(samp.permutations, rep(0, samp.perm))

    # Add Wills p to column headings:
    colnames(samp.permutations)[length(colnames(samp.permutations))] <- "p.Wills"

    # Add Wills p values:
    for(i in 1:samp.perm) samp.permutations[i, "p.Wills"] <- 1 - (sum(samp.permutations[i, "MIG"] < rand.permutations[, "MIG"]) / rand.perm)

  }

  # Add Wills p column:
  input.permutations <- cbind(input.permutations, rep(0, nrow(input.permutations)))

  # Add Wills p to column headings:
  colnames(input.permutations)[length(colnames(input.permutations))] <- "p.Wills"

  # Add Wills p values:
  for(i in 1:nrow(input.permutations)) input.permutations[i, "p.Wills"] <- 1 - (sum(input.permutations[i, "MIG"] < rand.permutations[, "MIG"]) / rand.perm)

  # Compile output as list:
  output <- list(input.permutations, samp.permutations, rand.permutations, trees, samp.trees, rand.trees)

  # Add names to output:
  names(output) <- c("input.tree.results", "samp.permutation.results", "rand.permutations", "input.trees", "samp.trees", "rand.trees")

  # Notify user that run is complete:
  cat(paste("Fit to stratigraphy measures for", nrow(input.permutations), "input trees have been calculated."))

  # Return output:
  return(invisible(output))

}

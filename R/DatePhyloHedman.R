#' Calculates branch lengths for a topology
#'
#' @description
#'
#' Calculates branch lengths for a topology given a tree and age data for the tips.
#'
#' @param tree A phylo object representing the tree the user wishes to time-scale.
#' @param tip.ages The ages of the tips of the tree to use for time-scaling.
#' @param outgroup.ages A vector of numeric values representing the ages of the outgroup taxa that fall immediately outside of the root node.
#' @param t0 The absolute maaximum age allowable. This must be older than anything in either tip.ages or outgroup.ages.
#' @param resolution The number of ages to sample from the posterior distributionof each inferred node age.
#' @param conservative A logical indicating whether or not to apply the conservatove approach of Lloyd et al. (2016). TRUE is the default and recommended option.
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
#' @export DatePhyloHedman
DatePhyloHedman <- function(
  tree,
  tip.ages,
  outgroup.ages,
  t0,
  resolution = 1000,
  conservative = TRUE
) {
  
  # TO DO - MORE CONDITIONALS NEEDED FOR WHEN RESOLUTION IS LOW
  
  # Function to find nodes that can be realistically dated by the Hedman technique:
  find.dateable.nodes <- function(tree, tip.ages) {
    
    # List all internal nodes except root:
    list.nodes <- (ape::Ntip(tree) + 2):(ape::Ntip(tree) + ape::Nnode(tree))
    
    # Create vector to store nodes that are dateable:
    defo.nodes <- vector(mode="numeric")
    
    # For each node in the list:
    for(i in length(list.nodes):1) {
      
      # Find descendant nodes:
      descs <- tree$edge[grep(TRUE, tree$edge[, 1] == list.nodes[i]), 2]
      
      # If all immediate descendants are internal nodes than it cannot be dated and is removed:
      if(length(which(descs > Ntip(tree))) == length(descs)) list.nodes <- list.nodes[-i]
      
      # If all immediate descendants are terminal nodes than it can definitely be dated...:
      if(length(which(descs <= Ntip(tree))) == length(descs)) {
        
        # ...and is retained...:
        defo.nodes <- c(defo.nodes, list.nodes[i])
        
        # ...but can be removed from list.nodes:
        list.nodes <- list.nodes[-i]
        
      }
      
    }
    
    # For each internal node with both a terminal and internal node descendant:
    for(i in length(list.nodes):1) {
      
      # Find node age:
      node.age <- max(tip.ages[FindDescendants(list.nodes[i], tree)])
      
      # List its descendants:
      descs.ages <- descs <- tree$edge[grep(TRUE, tree$edge[, 1] == list.nodes[i]), 2]
      
      # For each descendant:
      for(j in length(descs.ages):1) {
        
        # If an internal node date as oldest descendant:
        if(descs[j] > Ntip(tree)) descs.ages[j] <- max(tip.ages[FindDescendants(descs[j], tree)])
        
        # Remove if terminal node:
        if(descs[j] <= Ntip(tree)) descs.ages <- descs.ages[-j]
        
      }
      
      # If node age is older than any descendant internal nodes:
      if(node.age > max(descs.ages)) {
        
        # Then can be dated so add to list:
        defo.nodes <- c(defo.nodes, list.nodes[i])
        
      }
      
    }
    
    # Output all dateable nodes:
    return(defo.nodes)
    
  }
  
  # Function that retrieves dates for use as tnodes variable in Hedman function
  get.tnodes <- function(node, tree, tip.ages) {
    
    # Number of root node (point at which searching for ancestral nodes stops):
    root.node <- ape::Ntip(tree) + 1
    
    # If the input node is the root node:
    if(node == root.node) {
      
      # Error and warning:
      stop("ERROR: Input node is root node - tnodes has length one!")
      
      # If input node is not root node:
    } else {
      
      # Find first ancestor node:
      internodes <- tree$edge[match(node, tree$edge[, 2]), 1]
      
      # As long as we have not yet reached the root:
      while(internodes[length(internodes)] != root.node) {
        
        # Add next internode to set:
        internodes <- c(internodes, tree$edge[match(internodes[length(internodes)], tree$edge[, 2]), 1])
      }
      
      # Add in original node as that is counted too!:
      internodes <- c(node, internodes)
      
      # Create tnodes vector for storage:
      tnodes <- vector(mode="numeric")
      
      # Date first (original) node:
      tnodes[1] <- max(tip.ages[FindDescendants(internodes[1], tree)])
      
      # For each internode:
      for(i in 2:length(internodes)) {
        
        # Find descendants:
        descs <- tree$edge[grep(TRUE, tree$edge[, 1] == internodes[i]), 2]
        
        # Remove previously dated node:
        descs <- descs[-match(internodes[(i - 1)], descs)]
        
        # For each other descendant:
        for(j in 1:length(descs)) {
          
          # If descendant node is internal:
          if(descs[j] > Ntip(tree)) {
            
            # Get maximum possible age:
            descs[j] <- max(tip.ages[FindDescendants(descs[j], tree)])
            
            # If descendant node is terminal:
          } else {
            
            # Get maximum possible age:
            descs[j] <- tip.ages[descs[j]]
            
          }
          
        }
        
        # Store in tnodes:
        tnodes[i] <- max(descs)
        
      }
      
      # Vector for storing node ages to be deleted:
      deletes <- vector(mode="numeric")
      
      # For each potential tnode:
      for(i in 2:length(tnodes)) {
        
        # Find nodes too young for dating and store in deletes:
        if(max(tnodes[1:(i-1)]) > tnodes[i]) deletes <- c(deletes, i)
        
      }
      
      # Remove these nodes (if there are any) from tnodes:
      if(length(deletes) > 0) tnodes <- tnodes[-deletes]
      
      # Give in order from oldest to youngest:
      tnodes <- rev(tnodes)
      
      # Return tnodes:
      return(tnodes)
      
    }
    
  }
  
  # Check tree is fully bifurcating:
  if(!is.binary.tree(tree)) stop("Tree must be fully bifurcating.")
  
  # Ensure tip ages are in tree tip order:
  tip.ages <- tip.ages[tree$tip.label]
  
  # Find root node:
  root.node <- Ntip(tree) + 1
  
  # Create variables to store age estimates:
  age.estimates <- matrix(nrow = Nnode(tree), ncol = 3)
  
  # Create variables to store age distributions:
  age.distributions <- matrix(nrow = Nnode(tree), ncol = resolution)
  
  # Set column headings:
  colnames(age.estimates) <- c("Best.guess", "Best.guess.lower", "Best.guess.upper")
  
  # Set row names:
  rownames(age.estimates) <- rownames(age.distributions) <- c((Ntip(tree) + 1):(Ntip(tree) + Nnode(tree)))
  
  # Vector to store outgroup sequences:
  og.seq <- vector(mode = "numeric")
  
  # Report progress:
  cat("Identifying nodes that are date-able using the Hedman technique_")
  
  # Find dateable conservative nodes:
  if(conservative) dateable.nodes <- find.dateable.nodes(tree, tip.ages)
  
  # Find dateable non-conservative nodes:
  if(!conservative) dateable.nodes <- c((Ntip(tree) + 1):(Ntip(tree) + Nnode(tree)))
  
  # Function copied from Claddis rather than adding more dependencies:
  date_nodes <- function(time_tree) {

    # Get N tips:
    n_tips <- ape::Ntip(phy = time_tree)

    # Get N nodes:
    n_nodes <- ape::Nnode(phy = time_tree)

    # Store root node number:
    root_node <- n_tips + 1

    # If tree is a complete polytomy:
    if (time_tree$Nnode == 1) {

      # Create paths for just tips:
      paths <- as.list(x = 1:n_tips)

      # Add root to each path:
      for (i in 1:length(x = paths)) paths[[i]] <- c(paths[[i]], n_tips + 1)

    # If tree is not a complete polytomy:
    } else {

      # Create initial paths list with end nodes (terminal and internal, excluding the root):
      paths <- split(c(1:n_tips, (n_tips + 2):(n_tips + n_nodes)), f = 1:(n_tips + time_tree$Nnode - 1))

      # Strip names:
      names(paths) <- NULL

      # For each path:
      for (i in 1:length(x = paths)) {

        # Set counter as 1:
        j <- 1

        # Identify current node:
        current_node <- paths[[i]][j]

        # While current node is not the root (path has not terminated):
        while (current_node != root_node) {

          # Update current node and add to path:
          current_node <- paths[[i]][j + 1] <- time_tree$edge[match(current_node, time_tree$edge[, 2]), 1]

          # Update counter:
          j <- j + 1
        }
      }
    }

    # Create vector to store node ages:
    date_nodes <- vector(mode = "numeric", length = n_tips + time_tree$Nnode)

    # For each path:
    for (i in 1:length(x = paths)) {

      # Store path lengths from root:
      date_nodes[paths[[i]][1]] <- sum(time_tree$edge.length[match(paths[[i]][1:(length(x = paths[[i]]) - 1)], time_tree$edge[, 2])])
    }

    # Subtract path lengths from root time:
    date_nodes <- time_tree$root.time - date_nodes

    # Add node numbers:
    names(date_nodes) <- 1:(n_tips + time_tree$Nnode)

    # Return node ages:
    return(date_nodes)
  }

  # Report progress:
  cat("Done\nDating nodes using Hedman technique_")
  
  # If not using the conservative approach:
  if(!conservative) {
    
    # Create ages matrix:
    ages <- cbind(tip.ages, tip.ages)
    
    # Add rownames:
    rownames(ages) <- names(tip.ages)
    
    # Add column names:
    colnames(ages) <- c("FAD", "LAD")
    
    # Sort by taxon order in tree:
    ages <- ages[tree$tip.label, ]
    
    # Date tree using basic algorithm:
    stree <- DatePhylo(tree, ages, 0, "basic", FALSE)
    
    # Get node ages:
    sages <- date_nodes(stree)
    
    # Set first outgroup sequence (i.e. root) as this is special case:
    og.seq[1] <- paste(c(outgroup.ages, sages[root.node]), collapse = "%%")
    
    # For each non-root node:
    for(i in 2:length(dateable.nodes)) {
      
      # Identify node:
      node <- dateable.nodes[i]
      
      # Find preceding node:
      internodes <- tree$edge[match(node, tree$edge[, 2]), 1]
      
      # Keep going until we reach the root and add next internode to set:
      while(internodes[length(internodes)] != root.node) internodes <- c(internodes, tree$edge[match(internodes[length(internodes)], tree$edge[, 2]), 1])
      
      # Collate nodes and put in order:
      internodes <- sort(c(node, internodes))
      
      # Find outgroup age sequence and collaspe to single string and store:
      og.seq[i] <- paste(c(outgroup.ages, sages[internodes]), collapse="%%")
      
    }
    
  # If using the conservative approach:
  } else {
    
    # Find age estimates for each dateable node:
    for(i in 1:length(dateable.nodes)) {
      
      # Find tnodes for a specific node:
      tnodes <- get.tnodes(dateable.nodes[i], tree, tip.ages)
      
      # Add additional outgroup tnodes:
      tnodes <- c(outgroup.ages, tnodes)
      
      # Turn outgroup sequence into single string and store:
      og.seq[i] <- paste(tnodes, collapse="%%")
      
    }
    
  }
  
  # Get unique outgroup sequences:
  unq.og.seq <- unique(og.seq)
  
  # For each unique outgroup sequence:
  for(i in 1:length(unq.og.seq)) {
    
    # List nodes with this outgroup sequence:
    nodes <- dateable.nodes[which(og.seq == unq.og.seq[i])]
    
    # Define tnodes:
    tnodes <- as.numeric(strsplit(unq.og.seq[i], "%%")[[1]])
    
    # Get Hedman dates:
    hedman.out <- DateNodeHedman(tnodes = tnodes, t0 = t0, resolution = resolution)
    
    # Update age distributions:
    for(j in 1:length(nodes)) age.distributions[as.character(nodes[j]), ] <- hedman.out$Age.distribution
    
  }
  
  # Separate check to see if root is dateable (first find root descendants that are tips, if any):
  root.desc.tips <- sort(tree$tip.label[tree$edge[which(root.node == tree$edge[, 1]), 2]])
  
  # Now check that there are descendants that are tips:
  if(length(root.desc.tips) > 0) {
    
    # Get maximum root descendant age:
    root.tip.age <- max(tip.ages[root.desc.tips])
    
    # Now check that descendants are older than any other tips and get Hedman dates for root and update age distributions:
    if(root.tip.age > max(tip.ages[setdiff(tree$tip.label, root.desc.tips)])) age.distributions[as.character(root.node), ] <- DateNodeHedman(tnodes = c(outgroup.ages, root.tip.age), t0 = t0, resolution = resolution)$Age.distribution
    
  }
  
  # If not using the conservative approach:
  if(conservative == FALSE) {
    
    # Report progress:
    cat("Done\nIdentifying remaining undated nodes_")
    
    # Report progress:
    cat("Done\nDating remaining undated nodes using randomisation technique_")
    
    # If using the conservative approach:
  } else {
    
    # Report progress:
    cat("Done\nIdentifying remaining undated nodes_")
    
    # Find undated nodes by first establishing actually dated nodes:
    dated.nodes <- sort(dateable.nodes)
    
    # List all internal nodes:
    all.nodes <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))
    
    # Get undated nodes:
    undated.nodes <- setdiff(all.nodes, dated.nodes)
    
    # Find branches that define undated nodes:
    internal.branches <- tree$edge[which(tree$edge[, 2] > Ntip(tree)), ]
    
    # Is root dated? (default to "yes" before actually checking):
    root.dated <- "Yes"
    
    # If the root is undated (then it needs to be):
    if(is.na(match(root.node, dated.nodes))) {
      
      # Update root.dated:
      root.dated <- "no"
      
      # Create tnodes from outgroup ages ONLY (this will later be used to constrain the root age through the saem randomisation process as the other undated nodes):
      tnodes <- outgroup.ages
      
      # Date root node:
      hedman.out <- DateNodeHedman(tnodes = tnodes, t0 = t0, resolution = resolution)
      
      # Add to age estimates:
      age.estimates <- rbind(age.estimates, c(hedman.out$Best.guess, hedman.out$Best.guess.lower, hedman.out$Best.guess.upper))
      
      # Add to age distributions:
      age.distributions <- rbind(age.distributions, hedman.out$Age.distribution)
      
      # Add node name 0:
      rownames(age.estimates)[length(rownames(age.estimates))] <- rownames(age.distributions)[length(rownames(age.distributions))] <- "0"
      
    }
    
    # Vector for storing strings of sequences of undated node(s) bound by dated nodes:
    anc.strings.trimmed <- anc.strings <- vector(mode = "character")
    
    # Work through dated nodes to find sequences of undated nodes bracketed by them:
    for(i in 1:length(dated.nodes)) {
      
      # As long as the dated node is not the root (which has no ancestors):
      if(dated.nodes[i] != root.node) {
        
        # Get ancestors:
        ancestors <- internal.branches[which(internal.branches[, 2] == dated.nodes[i]), 1]
        
        # As long as the ancestor is neither the root or a dated node:
        if(ancestors != root.node && length(sort(match(ancestors, dated.nodes))) == 0) {
          
          # And as long as it remains so find next ancestral node:
          while(ancestors[length(ancestors)] != root.node && length(sort(match(ancestors[length(ancestors)], dated.nodes))) == 0) ancestors <- c(ancestors, internal.branches[match(ancestors[length(ancestors)], internal.branches[, 2]), 1])
          
        }
        
        # Add initial dated node:
        ancestors <- c(dated.nodes[i], ancestors)
        
        # Only need retain those with at least one undated node between dated nodes:
        if(length(ancestors) > 2) {
          
          # Add to ancestor strings:
          anc.strings <- c(anc.strings, paste(ancestors, collapse = " "))
          
          # Add to trimmed ancestor strings:
          anc.strings.trimmed <- c(anc.strings.trimmed, paste(ancestors[2:(length(ancestors) - 1)], collapse = " "))
          
        }
        
      }
      
    }
    
    # If root is undated then anc.strings with the root present must be modified to include node 0 at the end:
    if(root.dated == "no") {
      
      # Go through each anc.string:
      for(i in 1:length(anc.strings)) {
        
        # Get ancestors vector:
        ancestors <- as.numeric(strsplit(anc.strings[i], " ")[[1]])
        
        # If last node is the root:
        if(ancestors[length(ancestors)] == root.node) {
          
          # Add 0 node to end:
          ancestors <- c(ancestors, 0)
          
          # Re-store in anc.strings:
          anc.strings[i] <- paste(ancestors, collapse = " ")
          
        }
        
      }
      
    }
    
    # Find oldest node from each ancestral sequence (first set up empty vector):
    oldest.nodes <- vector(mode="numeric")
    
    # Go through each anc.string:
    for(i in 1:length(anc.strings)) {
      
      # Get ancestors vector:
      ancestors <- as.numeric(strsplit(anc.strings[i], " ")[[1]])
      
      # Store oldest node:
      oldest.nodes <- c(oldest.nodes, ancestors[length(ancestors)])
      
    }
    
    # Collapse to unique oldest nodes only:
    oldest.nodes <- sort(unique(oldest.nodes))
    
    # Clump anc.strings by shared oldest node (i.e. blocks that will be dated at the same time) - first create empty vector:
    clumped.anc.strings <- vector(mode="character")
    
    # For each oldest node find anc.strings with oldest node and collapse with %%:
    for(i in 1:length(oldest.nodes)) clumped.anc.strings <- c(clumped.anc.strings, paste(anc.strings[grep(paste(" ", oldest.nodes[i], sep = ""), anc.strings)], collapse = "%%"))
    
    # Report progress:
    cat("Done\nDating remaining undated nodes using randomisation technique_")
    
    # Can now date undated nodes using random draws from bounding Hedman dated nodes:
    for(i in 1:length(clumped.anc.strings)) {
      
      # First retrieve anc.strings with shared oldest node:
      anc.strings <- strsplit(clumped.anc.strings[i], "%%")[[1]]
      
      # Get young dated nodes (for drawing from for upper bounds):
      young.dated.nodes <- vector(mode="numeric")
      
      # For each ancestor string find youngest dated node, i.e., lower bounds::
      for(j in 1:length(anc.strings)) young.dated.nodes <- c(young.dated.nodes, as.numeric(strsplit(anc.strings[j], " ")[[1]][1]))
      
      # Get age distributions for young dated nodes:
      young.distributions <- age.distributions[as.character(young.dated.nodes), ]
      
      # Force into a single row matrix if only a vector:
      if(!is.matrix(young.distributions)) young.distributions <- t(as.matrix(young.distributions))
      
      # Set young distribution half up from young distributions:
      young.distributions.old.half <- young.distributions
      
      # Get old half as age distributions for young dated nodes:
      for(j in 1:length(apply(young.distributions, 1, median))) {
        
        # As long as the median is not the maximum:
        if(max(young.distributions[j, ]) > stats::median(young.distributions[j, ])) {
          
          # Split distribution into older half using median:
          temp.distribution <- young.distributions[j, which(young.distributions[j, ] > stats::median(young.distributions[j, ]))]
          
          # If the median is the maximum:
        } else {
          
          # Just use the maximum:
          temp.distribution <- max(young.distributions[j, ])
          
        }
        
        # ???:
        young.distributions.old.half[j, ] <- c(temp.distribution, rep(NA, length(young.distributions[j, ]) - length(temp.distribution)))
        
      }
      
      # Get old dated node (for drawing from for lower bounds):
      old.dated.node <- as.numeric(strsplit(anc.strings[1], " ")[[1]][length(strsplit(anc.strings[1], " ")[[1]])])
      
      # Get age distributions for old dated node:
      old.distribution <- age.distributions[as.character(old.dated.node), ]
      
      # As long as the median is not the minimum:
      if(min(old.distribution) < stats::median(old.distribution)) {
        
        # Split distribution using median:
        old.distribution.young.half <- old.distribution[which(old.distribution < stats::median(old.distribution))]
        
        # If the median is the minimum:
      } else {
        
        # Just use the minimum:
        old.distribution.young.half <- min(old.distribution)
        
      }
      
      # Get undated nodes (vector for storing output):
      undated.nodes <- vector(mode="numeric")
      
      # For each anc.string:
      for(j in 1:length(anc.strings)) {
        
        # Get ancestors:
        ancestors <- strsplit(anc.strings[j], " ")[[1]]
        
        # Store undated nodes:
        undated.nodes <- c(undated.nodes, as.numeric(ancestors[2:(length(ancestors) - 1)]))
        
      }
      
      # Collapse to unique nodes only:
      undated.nodes <- sort(unique(undated.nodes))
      
      # Find young constraints for each undated node (vector for storing results):
      young.constraints <- vector(mode = "numeric")
      
      # For each undated node add young constraints to list:
      for(j in 1:length(undated.nodes)) young.constraints <- c(young.constraints, paste(young.dated.nodes[grep(paste(" ", undated.nodes[j], sep = ""), anc.strings)], collapse = " "))
      
      # Find old constraints for each undated node (i.e. all nodes that are lower in tree as these will potentially be dated first and replace the current oldest age; vector for storing results):
      old.constraints <- vector(mode="numeric")
      
      # For each undated node:
      for(j in 1:length(undated.nodes)) {
        
        # Get just anc.strings where undated node is present:
        node.strings <- anc.strings[grep(paste(" ", undated.nodes[j], sep=""), anc.strings)]
        
        # Vector for storing older nodes:
        older.nodes <- vector(mode="numeric")
        
        # For each anc.string where the undated node exists ?????:
        for(k in 1:length(node.strings)) older.nodes <- c(older.nodes, as.numeric(strsplit(strsplit(node.strings[k], undated.nodes[j])[[1]][2], " ")[[1]]))
        
        # ????:
        old.constraints <- c(old.constraints, paste(unique(sort(older.nodes)), collapse=" "))
        
      }
      
      # Find undated nodes constrained by dated nodes (vector for storage):
      young.dated.nodes.constrain <- vector(mode = "character")
      
      # For each young dated node constraint find undated nodes which it contrains:
      for(j in 1:length(young.dated.nodes)) young.dated.nodes.constrain[j] <- paste(as.numeric(strsplit(anc.strings[j], " ")[[1]][2:(length(strsplit(anc.strings[j], " ")[[1]]) - 1)]), collapse = " ")
      
      # Main dating loop (repeats random draw N times, where N is defined by the variable resolution):
      for(j in 1:resolution) {
        
        # Modify age distributions used if medians of undated nodes exceed bounds of medians of dated nodes (needs to have been at least two entries already):
        if(j >= 3 && j >= floor(resolution / 2)) {
          
          # If there is more than one undated node establish present medians for undated nodes:
          if(length(undated.nodes) > 1) undated.medians <- apply(age.distributions[as.character(undated.nodes), 1:(j - 1)], 1, median)
          
          # If there is only one undated node
          if(length(undated.nodes) == 1) {
            
            # Establish present median of undated node:
            undated.medians <- stats::median(age.distributions[as.character(undated.nodes), 1:(j - 1)])
            
            # Add name for reference later:
            names(undated.medians) <- as.character(undated.nodes)
            
          }
          
          # Case if oldest median of the undated nodes exceeds the median of the bounding old dated node:
          if(max(undated.medians) >= stats::median(sort(old.distribution))) {
            
            # Update active distribution with older half only (to force undated nodes towards median ages consistent with dated nodes):
            active.old.distribution <- old.distribution.young.half
            
            # If oldest median ages do not conflict:
          } else {
            
            # Draw from full distribution:
            active.old.distribution <- old.distribution
            
          }
          
          # Set default active young distributions:
          active.young.distributions <- young.distributions
          
          # For each young (top bounding) dated node:
          for(k in 1:length(young.dated.nodes)) {
            
            # Case if youngest median of the undated nodes exceeds the median of the bounding young dated node - update active distribution with younger half only (to force undated nodes towards median ages consistent with dated nodes):
            if (min(undated.medians[strsplit(young.dated.nodes.constrain[k], " ")[[1]]]) <= stats::median(sort(age.distributions[as.character(young.dated.nodes[k]), ]))) active.young.distributions[k, ] <- young.distributions.old.half[k, ]
            
          }
          
          # Case if still in first half of resolution:
        } else {
          
          # Use uncorrected distribution for upper bound:
          active.old.distribution <- old.distribution
          
          # Use uncorrected distributions for lower bound:
          active.young.distributions <- young.distributions
          
        }
        
        # Special case of last iteration where we want to ensure we draw a date between the constraining medians:
        if(j == resolution) {
          
          # Set old distribution to young half only:
          active.old.distribution <- old.distribution.young.half
          
          # Set default active young distributions:
          active.young.distributions <- young.distributions
          
          # For each young (top bounding) dated node set all young distribution to old half only:
          for(k in 1:length(young.dated.nodes)) active.young.distributions[k, ] <- young.distributions.old.half[k, ]
          
        }
        
        # Modify young and old distributions to remove tails which will always violate node order (young node older than old node and vice versa) to speed up the random draw step below; if part of young distributions are older than the oldest part of the old distribution:
        if(max(sort(as.vector(active.young.distributions))) >= max(sort(active.old.distribution))) {
          
          # Find oldest part of old distribution
          upper.limit <- max(sort(active.old.distribution))
          
          # ????:
          active.young.distributions[which(active.young.distributions >= upper.limit)] <- NA
          
        }
        
        # Vector to store young distribution minima:
        active.young.distribution.mins <- vector(mode="numeric")
        
        # For each young node store minima:
        for(k in 1:length(active.young.distributions[, 1])) active.young.distribution.mins <- c(active.young.distribution.mins, min(sort(active.young.distributions[k, ])))
        
        # If part of old distribution is younger than the youngest part of the youngest young distribution:
        if(min(sort(active.old.distribution)) <= min(active.young.distribution.mins)) {
          
          # Find lower limit (minimum of minima):
          lower.limit <- min(active.young.distribution.mins)
          
          # Remove values less than or equal to the lower limit from the old distribution:
          active.old.distribution <- active.old.distribution[grep(FALSE, active.old.distribution <= lower.limit)]
          
        }
        
        # Draw upper and lower bounds at random from old and young nodes:
        young.random.ages <- vector(mode="numeric")
        
        # For each young node draw a random age from its distributon:
        for(k in 1:length(young.dated.nodes)) young.random.ages <- c(young.random.ages, sort(active.young.distributions[k, ])[ceiling(stats::runif(1, 0, length(sort(active.young.distributions[k, ]))))])
        
        # Repeat for old age:
        old.random.age <- sort(active.old.distribution)[ceiling(stats::runif(1, 0, length(sort(active.old.distribution))))]
        
        # Ensure old node is older than all young nodes; while old node is younger or equal in age to oldest young nodes:
        while(old.random.age <= max(young.random.ages)) {
          
          # Re-draw ages as above:
          young.random.ages <- vector(mode="numeric")
          
          # For each young node draw a random age from its distributon:
          for(k in 1:length(young.dated.nodes)) young.random.ages <- c(young.random.ages, sort(active.young.distributions[k, ])[ceiling(stats::runif(1, 0, length(sort(active.young.distributions[k, ]))))])
          
          # Repeat for old age:
          old.random.age <- sort(active.old.distribution)[ceiling(stats::runif(1, 0, length(sort(active.old.distribution))))]
          
        }
        
        # Need to add names so can extract data later:
        names(young.random.ages) <- young.dated.nodes
        
        # Find oldest (i.e. constraining) young node for each undated node (vector for storing output):
        young.constraints.node <- young.constraints.age <- vector(mode="numeric")
        
        # For each undated node:
        for(k in 1:length(undated.nodes)) {
          
          # Store young age:
          young.constraints.age <- c(young.constraints.age, max(young.random.ages[strsplit(young.constraints[k], " ")[[1]]]))
          
          # Store young node:
          young.constraints.node <- c(young.constraints.node, as.numeric(names(young.random.ages[strsplit(young.constraints[k], " ")[[1]]])[which(young.random.ages[strsplit(young.constraints[k], " ")[[1]]] == max(young.random.ages[strsplit(young.constraints[k], " ")[[1]]]))[1]]))
          
        }
        
        # Vector for storing node dates:
        node.dates <- rep(NA, length(undated.nodes))
        
        # Date undated nodes using randomisation process:
        for(k in 1:length(unique(young.constraints.node)) ) {
          
          # Find top.node:
          top.node <- unique(young.constraints.node)[k]
          
          # Find target (.e. undated) nodes:
          target.nodes <- undated.nodes[which(young.constraints.node == top.node)]
          
          # Find bounding dates (potential old constraints):
          old.constraining.nodes <- as.numeric(sort(unique(strsplit(paste(old.constraints[match(target.nodes, undated.nodes)], collapse=" "), " ")[[1]])))
          
          # Remove any target nodes present
          if(length(sort(match(target.nodes, old.constraining.nodes))) > 0) old.constraining.nodes <- old.constraining.nodes[-sort(match(target.nodes, old.constraining.nodes))]
          
          # Find youngest old node as actual constraint:
          old.constraining.node <- max(old.constraining.nodes)
          
          # If lower bound is the previously dated old node:
          if(old.constraining.node == old.dated.node) {
            
            # Use that as maximum limit
            max <- old.random.age
            
            # If lower bound is a previously undated node:
          } else {
            
            # Use that as maximum limit
            max <- node.dates[match(old.constraining.node, undated.nodes)]
            
          }
          
          # Set minimum age:
          min <- young.random.ages[as.character(top.node)]
          
          # Get node dates:
          node.dates[match(target.nodes, undated.nodes)] <- sort(stats::runif(length(target.nodes), min = min, max = max), decreasing = TRUE) # Draw node ages from bounding minima and maxima
          
        }
        
        # Store results:
        age.distributions[as.character(undated.nodes), j] <- node.dates
        
      }
      
    }
    
    # Can now remove node "0" (i.e., if root is undated by Hedman method) if present:
    if(root.dated == "no") {
      
      # Remove zero node from distributions:
      age.distributions <- age.distributions[-which(rownames(age.distributions) == "0"), ]
      
      # Remove zero node from estimates:
      age.estimates <- age.estimates[-which(rownames(age.estimates) == 0), ]
      
    }
    
  }
  
  # Report progress:
  cat("Done\nTidying up and returning results_")
  
  # Create vector of nodes estimated using Hedman method:
  Hedman.estimated <- rep(1, length(rownames(age.estimates)))
  
  # Update those not estimated using Hedman method:
  if(conservative == TRUE) Hedman.estimated[match(setdiff(rownames(age.estimates), dateable.nodes), rownames(age.estimates))] <- 0
  
  # Add to age estimates:
  age.estimates <- cbind(age.estimates, Hedman.estimated)
  
  # Sort age distributions in advance of picking confidence intervals:
  age.distributions <- t(apply(age.distributions, 1, sort))
  
  # Add median values to age estimates:
  age.estimates[, "Best.guess"] <- apply(age.distributions, 1, median)
  
  # Add upper CI to age estimates:
  age.estimates[, "Best.guess.upper"] <- age.distributions[, ceiling(0.975 * resolution)]
  
  # Add lower CI to age estimates:
  age.estimates[, "Best.guess.lower"] <- age.distributions[, max(c(1, floor(0.025 * resolution)))]
  
  # Create vectors to store node ages for each branch:
  to.ages <- from.ages <- age.estimates[match(tree$edge[, 1], rownames(age.estimates)), "Best.guess"]
  
  # Get to ages for terminal branches:
  to.ages[which(tree$edge[, 2] <= Ntip(tree))] <- tip.ages[tree$edge[which(tree$edge[, 2] <= Ntip(tree)), 2]]
  
  # Get to ages for internal branches:
  to.ages[which(tree$edge[, 2] > Ntip(tree))] <- age.estimates[match(tree$edge[which(tree$edge[, 2] > Ntip(tree)), 2], rownames(age.estimates)), "Best.guess"]
  
  # Update tree with branch lengths scaled to time using median ages:
  tree$edge.length <- from.ages-to.ages
  
  # Update tree to include root time:
  tree$root.time <- age.estimates[as.character(root.node), "Best.guess"]
  
  # Collate results:
  results <- list(age.estimates, age.distributions, tree)
  
  # Add names:
  names(results) <- c("age.estimates", "age.distributions", "tree")
  
  # Repprt progress
  cat("Done")
  
  # Return output:
  return(results)
  
}


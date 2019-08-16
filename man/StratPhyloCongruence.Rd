\name{StratPhyloCongruence}

\alias{StratPhyloCongruence}

\title{Calculates fit to stratigraphy metrics for a set of tree(s).}

\description{Calculates SCI, RCI, MSM*, and GER for a number of topologies.}

\usage{StratPhyloCongruence(trees, ages, rlen=0, method="basic", samp.perm=1000,
	rand.perm=1000, hard=TRUE, randomly.sample.ages=FALSE, fix.topology=TRUE,
	fix.outgroup=TRUE)}

\arguments{
\item{trees}{Input tree(s) as either a phylo or multiphylo object.}
\item{ages}{A two-column matrix of taxa (rows) against First and Last Appearance Datums (FADs and LADs) to be passed to \link{DatePhylo}. Note that rownames should be the taxon names exactly as they appear in \code{tree$tip.label} and colnames should be "FAD" and "LAD". All ages should be in time before present.}
\item{rlen}{Root length, to be passed to \link{DatePhylo}.}
\item{method}{Tree dating method, to be passed to \link{DatePhylo}.}
\item{samp.perm}{Number of sampled trees to be produced by resolving polytomies and/or drawing random dates for the tips for the input trees.}
\item{rand.perm}{Number of random trees to be produced in calculating probabilities for the input trees, and (if used) the sampled trees.}
\item{hard}{Whether to treat polytomies as hard or soft. If \code{FALSE} polytomies are resolved randomly.}
\item{randomly.sample.ages}{Whether to treat FAD and LAD as a range (\code{randomly.sample.ages=FALSE}) or an uncertainty (\code{randomly.sample.ages=TRUE}). If the latter then two ages are randomly sampled from the range and these are used as the FAD and LAD.}
\item{fix.topology}{Whether to allow tree shape to be random (\code{fix.topology=FALSE}) or to reflect the tree shape of the input tree(s) (\code{fix.topology=TRUE}).}
\item{fix.outgroup}{Whether to force the randomly generated trees to share the outgroup of the first input tree (\code{fix.outgroup=TRUE}) or not (\code{fix.outgroup=FALSE}).}
}

\details{
Cladograms of fossil taxa make explicit predictions about the successive appearance of taxa in the fossil record that can be compared with their observed stratigraphic ranges. Several methods have been developed to quantify this "fit to stratigraphy" of phylogenetic hypotheses, and these can be assessed in both a statistic (measuring the apparent strength of this congruence) and an associated significance test (p-value) based on generating random topologies for the same taxon set.

This function produces both values for all four main metrics: the Stratigraphic Consistency Index (SCI; Huelsenbeck 1994), the Relative Consistency Index (RCI; Benton and Storrs 1994), the Manhattan Stratigraphic Measure (MSM*; Siddall 1998; Pol and Norell 2001), and the Gap Excess Ratio (GER; Wills 1999).

\bold{SCI - Stratigraphic Consistency Index}

	The SCI works by assessing the "consistency" of nodes. A node is considered stratigraphically consistent if its oldest descendant is the same age or younger than the oldest descendant of the preceding node. The SCI is thus given simply as:

		\deqn{SCI = C / N}

	Where \eqn{C} is the sum of all the consistent nodes and \eqn{N} is the total number of nodes - 1. (As there is no node preceding the root there is no basis on which to estimate its consistency.) This value can range from zero (maximally inconsistent) to one (maximally consistent). However, a potential criticism of the SCI is that a high value may be returned when in fact a single inconsistent node may represent a very large amount of missing history (measured in unsampled units or millions of years), whereas a low SCI may represent relatively few unsampled units or millions of years.

\bold{RCI - Relative Completeness Index}

	The RCI was the first method to explicitly account for the absolute amount of missing data implied by the tree. This figure is usually expressed as the Minimum Implied Gap (\eqn{MIG}), a term also used by both the MSM and GER (see below), and corresponds to the sum of the branch lengths excluding the duration of the terminals (the observed ranges of the taxa). The RCI expresses the MIG as a proportion of the sum of the observed ranges (Simple Range Length; \eqn{SRL}) of the taxa converted to a percentage:

		\deqn{RCI = (1 - (MIG / SRL)) * 100 percent}

	Importantly this value is not confined to a 0 to 100 percent scale, and can have both negative values and values greater than 100 percent, which can make it difficult to interpret.

\bold{MSM - Manhattan Stratigraphic Measure}

	The MSM was the first method to account for both the absolute \eqn{MIG} and range on a confined zero to one scale. It is expressed as:

		\deqn{MSM = L_m / L_0}

	Where \eqn{L_0} is the length of the tree expressed by optimising times of first appearance on to the tree as a Sankoff character and taking the total length. \eqn{L_m} represents the same process, but for the optimal possible tree given the same set of first appearances. However, Pol and Norell (2001) noted a critical flaw in this approach, specifically that the Sankoff optimisation is reversible, meaning that nodes in the topology are allowed to be younger than their descendants, leading in some cases to a poor fit to stratigraphy being perceived as a good fit. Instead they suggest modifying the character step matrix to make the cost of reversals effectively infinite and hence impossible. Thus the values for \eqn{L_0} and \eqn{L_m} are modified accordingly. This approach they termed MSM* and is the implementation of MSM used here. This statistic can be expressed as:

		\deqn{MSM* = G_min / MIG}

	Where \eqn{G_min} represents the \eqn{MIG} for the tree with the optimal fit to stratigraphy. In effect this is a completely unbalanced tree where the youngest pair of taxa are the most deeply nested and successive outgroups represent the next oldest taxon. Theoretically MSM* ranges from one (the best fit as the observed tree is the maximally consistent tree) to zero (the least optimal tree). However, in effect no tree can have a value of zero as its \eqn{MIG} would have to be equal to infinity.

\bold{GER - Gap Excess Ratio}

	The GER represents a method that accounts for \eqn{MIG}, ranges from zero to one, and the best and worst fits to stratigraphy are both practically realisable. It can be expressed as:

		\deqn{GER = 1 - ((MIG - G_min) / (G_max - G_min))}

	Where \eqn{G_max} represents the \eqn{MIG} of the tree with the worst possible fit to stratigraphy. This is in effect any topology where the oldest taxon is the most deeply nested such that every clade in the tree contains it and hence must be minimally that old.

\bold{P-values}

	In isolation all four methods suffer from an inability to reject the null hypothesis that an apparent good fit to stratigraphy may be generated by chance alone. In practice this can be tested by generating a set of random topologies, calculating the fit to stratigraphy measure, and then either fitting a normal distribution to the resulting values (to get an estimated p-value) or assessing the relative position of the MIG of the observed (and sampled) tree(s) to get an absolute p-value. (Note that the assumption of normality may not always hold and the former approach should be used at the user's discretion. However, it should be noted that for the SCI, MSM*, and GER p-values are calculated after first transforming the data by taking the arcsine of the square root of each value.) The reason for having two sets of p-values is that if the observed trees fall completely outside the range of the random topologies they will be given an extreme p-value (0 or 1) that may be misleading. In such cases the estimated value may be more accurate.

	P-values should be interpreted as the probability of the null: that the observed tree(s) have an equal or worse fit to stratigraphy than the sample of random trees. Thus if the p-values are very small the user can reject the null hypothesis in favour of the alternative: that the observed tree(s) have a better fit to stratigraphy than expected by chance alone.

\bold{Modifications of the GER}

	More recently Wills et al. (2008) introduced two new versions of the GER that take advantage of the distribution of \eqn{MIG}s from the set of randomly generated topologies. The first of these (GERt) uses the extreme values of the random topologies as modified versions of \eqn{G_max} and \eqn{G_min}, termed \eqn{Gt_max} and \eqn{Gt_min} respectively. GERt is thus expressed as:

		\deqn{GERt = 1 - ((MIG - Gt_min) / (Gt_max - Gt_min))}

	In practice the \eqn{MIG} of the observed tree(s) may fall outside of these ranges so here a correction factor is employed so that any value below zero is corrected to zero, and any value above one is corrected to one. An additional stipulation for GERt is that the overall tree topology is fixed and only the taxa themselves are shuffled. This is to give a more realistic set of random topologies as there are known biases towards unbalanced trees in many palaeontological data sets. Here this is implemented by selecting the \code{fix.topology=TRUE} option. However, here GERt can also be calculated when \code{fix.topology=FALSE}.

	A second modification of GER is to use the position of the observed tree(s) in the sample of randomly generated topologies, such that:

		\deqn{GER* = 1 - (Fraction of distribution <= MIG)}

	Thus if the \eqn{MIG} of the observed tree(s) is less than any randomly generated topology GER* will be one (maximally optimal fit) and if it worse than any of the randomly generated topologies it will be zero (maximally suboptimal fit).

	Note: it is recommended that you use a large number of random topologies in order to get reliable values for GERt and GER* using \code{rand.perm=N}. Wills et al. (2008) used 50000, but the user should note that for many real world examples such values will take many hours to run.

	Wills et al. (2008) also introduced the notion of referring to intervals sampled rather than absolute time by recasting \eqn{MIG} as \eqn{MIGu}: the sum of ghost ranges for intervals of unit length. Although not directly implemented here this can be done manually by converting the time values (in Ma) used to simple unit counts such that FADs and LADs of taxa are given as numbered time bins (the youngest being 1 and the oldest N, where there are N time bins).

\bold{Polytomies and age uncertainties}

	Alongside the input trees the user can also create an additional set of sampled trees based on the input trees. This option is automatically implemented when choosing either \code{hard=FALSE} or \code{randomly.sample.ages=TRUE}, and the total number of permutations to perform dictated by \code{samp.perm=N}. This process works by first sampling from the set of input tree(s) and then randomly resolving any polytomies (if \code{hard=FALSE}) to ensure all sampled trees are fully dichotomous. (At present the function does not allow the various options laid out in Boyd et al. 2011, but the user can achieve this effect by modifying the input trees themselves.) Then if \code{randomly.sample.ages=TRUE} the FAD and LAD are treated as bounds of a uniform distribution which is sampled at random. This allows the user to get results for a set of trees that account for uncertainty in dating (as outlined in Pol and Norell 2006). (Note that two dates are picked for each taxon to avoid the problem of having an SRL of zero that would cause a divide by zero error for the RCI metric.)

	All fit to stratigraphy measures calculated for the input trees are then repeated for the sampled trees. However, if both \code{hard=TRUE} and \code{randomly.sample.ages=FALSE} (the defaults) no set of sampled trees will be created.

	In all cases when using the function users will see a progress bar that indicates the general progress through the combined set of trees (input, sampled, and randomly generated). This serves as a useful indicator of the time it will take for the function to finish. Here default values for \code{samp.perm} and \code{rand.perm} are both set at 1000, but the user may wish to lower these (to decrease calculation time) or increase them (to enhance accuracy).

\bold{Additional options}

	Note that because this function uses \code{\link{DatePhylo}} the user has the option of using different tree dating algorithms than the \code{basic} method (equivalent to the \code{basic} method in the \code{paleotree} package) employed in all the published studies cited above (and the default option here). The dating method used will apply to all trees generated, including the input, sampled, and randomly generated topologies. In all cases the time-scaled trees are returned with the function output.

	A final option (\code{fix.outgroup=TRUE}) allows the user to always use the same outgroup taxon (based on the first input tree) for all randomly generated topologies. Because the outgroup will often be the oldest taxon and its position in the input topologies is not allowed to vary letting it do so in the random topologies may lead to inferring a better fit to stratigraphy for the observed tree(s) than is fair. Fixing the outgroup thus ameliorates this potential bias and is the default option here.
}

\value{
\item{input.tree.results}{A matrix with a row for each input tree and columns indicating the values for SCI, RCI, GER and MSM* and their estimated probabilities assuming a normal distribution (est.p.SCI, est.p.RCI, est.p.GER, and est.p.MSM*) as well as GERt, GER*, MIG, and p.Wills (their probability as position within the MIGs of the random topologies).}
\item{samp.permutation.results}{If used, a matrix with a row for each sampled tree (up to \code{samp.perm}) and columns indicating the values for SCI, RCI, GER and MSM* and their estimated probabilities assuming a normal distribution (est.p.SCI, est.p.RCI, est.p.GER, and est.p.MSM*) as well as GERt, GER*, MIG, and p.Wills (their probability as position within the MIGs of random topologies).}
\item{rand.permutations}{A matrix with a row for each randomly generated tree (up to \code{rand.perm}) and columns indicating the values for SCI, RCI, GER, MSM*, and MIG.}
\item{input.trees}{The input tree(s) as a phylo or multiphylo object, with branches scaled to time according to the input values passed to \link{DatePhylo}.}
\item{samp.trees}{The sampled tree(s) as a phylo or multiphylo object, with branches scaled to time according to the input values passed to \link{DatePhylo}.}
\item{rand.trees}{The randomly generated tree(s) as a phylo or multiphylo object, with branches scaled to time according to the input values passed to \link{DatePhylo}.}
}

\references{
Benton, M. J. and Storrs, G. W., 1994. Testing the quality of the fossil record: palaeontological knowledge is improving. Geology, 22, 111-114.

Boyd, C. A., Cleland, T. P., Marrero, N. L. and Clarke, J. A., 2011. Exploring the effects of phylogenetic uncertainty and consensus trees on stratigraphic consistency scores: a new program and a standardized method. Cladistics, 27, 52-60.

Huelsenbeck, J. P., 1994. Comparing the stratigraphic record to estimates of phylogeny. Paleobiology, 20, 470-483.

Pol, D. and Norell, M. A., 2001. Comments on the Manhattan Stratigraphic Measure. Cladistics, 17, 285-289.

Pol, D. and Norell, M. A., 2006. Uncertainty in the age of fossils and the stratigraphic fit to phylogenies. Systematic Biology, 55, 512-521.

Siddall, M. E., 1998. Stratigraphic fit to phylogenies: a proposed solution. Cladistics, 14, 201-208.

Wills, M. A., 1999. Congruence between phylogeny and stratigraphy: randomization tests and the Gap Excess Ratio. Systematic Biology, 48, 559-580.

Wills, M. A., Barrett, P. M. and Heathcote, J. F., 2008. The modified Gap Excess Ratio (GER*) and the stratigraphic congruence of dinosaur phylogenies. Systematic Biology, 57, 891-904.
}

\author{Mark A. Bell \email{mark.bell521@gmail.com} and Graeme T. Lloyd \email{graemetlloyd@gmail.com}}

\keyword{congruence}

\examples{
# Calculate stratigraphic fit measures treating ages as ranges
# (permutation numbers used are lower than recommended for standard use):
fit.to.strat.1 <- StratPhyloCongruence(trees=Dipnoi$tree, ages=Dipnoi$ages, rlen=0,
	method="basic", samp.perm=5, rand.perm=5, hard=TRUE,
	randomly.sample.ages=FALSE, fix.topology=TRUE, fix.outgroup=TRUE)
fit.to.strat.1 # View all output
names(fit.to.strat.1) # Show output options
fit.to.strat.1$input.tree.results # Show just the output for the input tree(s)

# Calculate stratigraphic fit measures treating ages as uncertainties
# (permutation numbers used are lower than recommended for standard use):
#fit.to.strat.2 <- StratPhyloCongruence(trees=Dipnoi$tree, ages=Dipnoi$ages, rlen=0,
#	method="basic", samp.perm=10, rand.perm=10, hard=TRUE,
#	randomly.sample.ages=TRUE, fix.topology=TRUE, fix.outgroup=TRUE)
}

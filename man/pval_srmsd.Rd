% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pval_srmsd.R
\name{pval_srmsd}
\alias{pval_srmsd}
\title{Signed RMSD measure of null p-value uniformity}
\usage{
pval_srmsd(pvals, causal_indexes, detailed = FALSE)
}
\arguments{
\item{pvals}{The vector of association p-values to analyze.
This function assumes all p-values are provided (a mix of null and alternative tests).
\code{NA} values are allowed in input and removed.
Non-\code{NA} values outside of [0, 1] will trigger an error.}

\item{causal_indexes}{The vector of causal indexes, whose p-values will be omitted.
Values of \code{causal_indexes} as returned by \code{sim_trait} work.
This parameter is required to prevent use of this function except when the true status of every test (null vs alternative) is known.
Set to \code{NULL} if all loci are truly null (non-causal).
Otherwise, \code{causal_indexes} must have at least one causal index.}

\item{detailed}{If \code{FALSE} (default) only SRMSD is returned.
If \code{TRUE}, sorted null p-values without NAs and their expectations are returned (useful for plots).}
}
\value{
If \code{detailed} is \code{FALSE}, returns the signed RMSD between the observed p-value order statistics and their expectation under true uniformity.
If \code{detailed} is \code{TRUE}, returns data useful for plots, a named list containing:
\itemize{
\item \code{srmsd}: The signed RMSD between the observed p-value order statistics and their expectation under true uniformity.
\item \code{pvals_null}: Sorted null p-values (observed order statistics).  If any input null p-values were \code{NA}, these have been removed here (removed by \code{\link[=sort]{sort()}}).
\item \code{pvals_unif}: Expected order statistics assuming uniform distribution, same length as \code{pvals_null}.
}

If the input \code{pvals} is \code{NULL} (taken for case of singular association test, which is rare but may happen), then the returned value is \code{NA} if \code{detailed} was \code{FALSE}, or otherwise the list contains \code{NA}, \code{NULL} and \code{NULL} for the above three items.
}
\description{
Quantifies null p-value uniformity by computing the RMSD (root mean square deviation) between the sorted observed null (truly non-causal) p-values and their expected quantiles under a uniform distribution.
Meant as a more robust alternative to the "inflation factor" common in the GWAS literature, which compares median values only and uses all p-values (not just null p-values).
Our signed RMSD, to correspond with the inflation factor, includes a sign that depends on the median null p-value:
positive if this median is \verb{<= 0.5} (corresponds with test statistic inflation), negative otherwise (test statistic deflation).
Zero corresponds to uniform null p-values, which arises in expectation only if test statistics have their assumed null distribution (there is no misspecification, including inflation).
}
\examples{
# simulate truly null p-values, which should be uniform
pvals <- runif(10)
# for toy example, take the first p-value to be truly causal (will be ignored below)
causal_indexes <- 1
# calculate desired measure
pval_srmsd( pvals, causal_indexes )

}
\seealso{
\code{\link[=rmsd]{rmsd()}} for the generic root-mean-square deviation function.

\code{\link[=pval_infl]{pval_infl()}} for the more traditional inflation factor, which focuses on the median of the full distribution (combination of causal and null cases).

\code{\link[=pval_type_1_err]{pval_type_1_err()}} for classical type I error rate estimates.
}

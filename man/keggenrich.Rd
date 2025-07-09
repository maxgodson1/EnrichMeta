\name{keggenrich}
\alias{keggenrich}
\title{Perform KEGG Metabolite Enrichment Analysis}
\usage{
keggenrich(
  KEGGid,
  species,
  pathway_data,
  p.adjust.method = "BH"
)
}
\arguments{
\item{KEGGid}{Character vector of metabolite KEGG IDs (e.g., c("C00031", "C00022"))}

\item{species}{KEGG species code (e.g., "hsa" for human)}

\item{pathway_data}{Pathway data list obtained from \code{\link{getkpd}}}

\item{p.adjust.method}{Method for p-value adjustment (default="BH"). 
Options: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".}
}
\value{
A data frame with enrichment results containing the following columns:
\itemize{
  \item \code{pathwayID}: KEGG pathway identifier
  \item \code{Description}: Pathway description
  \item \code{MetaRatio}: Ratio of metabolites in input list found in pathway (format k/N)
  \item \code{BgRatio}: Background ratio (format n/M)
  \item \code{pvalue}: Raw p-value from hypergeometric test
  \item \code{EnrichmentRatio}: Enrichment ratio (k*M)/(n*N)
  \item \code{keggID}: KEGG IDs of overlapping metabolites (semicolon separated)
  \item \code{Count}: Number of overlapping metabolites
  \item \code{p.adjust}: Adjusted p-value
}
}
\description{
Performs metabolite set enrichment analysis using KEGG pathways and hypergeometric test.
Identifies pathways significantly enriched with input metabolites.
}
\details{
This function takes a list of metabolite KEGG IDs and tests for enrichment in KEGG pathways.
It requires pre-fetched pathway data from \code{\link{getkpd}}. The analysis uses:
\itemize{
  \item N: Number of input metabolites
  \item M: Total metabolites in background (all pathway metabolites)
  \item n: Metabolites in a specific pathway
  \item k: Overlap between input and pathway metabolites
}
The hypergeometric test formula: P(X >= k) ~ phyper(k-1, n, M-n, N, lower.tail=FALSE)
}
\examples{
\dontrun{
# First get pathway data
pathway_data <- getkpd(species = "hsa", cache_dir = "cache")

# Perform enrichment analysis
metabolites <- c("C00031", "C00022", "C00041", "C00033", "C00092")
results <- keggenrich(KEGGid = metabolites, 
                      species = "hsa", 
                      pathway_data = pathway_data)

# View top results
head(results[order(results$p.adjust), ])

# Visualize results
enrichdot(results)
enrichbar(results)
}
}
\seealso{
\code{\link{getkpd}} to obtain pathway data, 
         \code{\link{enrichdot}} and \code{\link{enrichbar}} for visualization
}
\author{
Your Name
}
\keyword{enrichment}
\keyword{metabolomics}

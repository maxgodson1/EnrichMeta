\name{pathwaynetwork}
\alias{pathwaynetwork}
\title{Visualize KEGG Pathway Network}
\usage{
pathwaynetwork(
  pathway_data,
  min_shared = 1,
  layout = "fr",
  point_size = 3,
  label_size = 1,
  edge_width = 0.5,
  vertex_color = "#8DA0CB",
  vertex_alpha = 0.7,
  edge_color = "#666666",
  edge_alpha = 0.5,
  main_title = "KEGG Pathway Network",
  sub_title = NULL,
  legend_position = "bottomright",
  show_legend = TRUE,
  margin = c(0, 0, 0, 0),
  legend_inset = c(0.02, 0.02)
)
}
\arguments{
\item{pathway_data}{Pathway data list obtained from \code{\link{getkpd}} or RDS file.
Must contain:
\itemize{
  \item \code{pathways}: Named vector of pathway IDs to pathway names
  \item \code{pathway2compound}: List mapping pathway IDs to compound IDs
}}

\item{min_shared}{Minimum number of shared compounds to create an edge between pathways (default=1)}

\item{layout}{Network layout algorithm (default="fr"). Options:
\itemize{
  \item "fr": Fruchterman-Reingold force-directed layout
  \item "kk": Kamada-Kawai force-directed layout
  \item "dh": Davidson-Harel simulated annealing layout
  \item "circle": Circular layout
}}

\item{point_size}{Base multiplier for vertex size calculation (default=3). 
Actual size = sqrt(compound_count) * point_size + 3}

\item{label_size}{Vertex label size multiplier (default=1)}

\item{edge_width}{Multiplier for edge width (default=0.5). 
Actual width = shared_count * edge_width}

\item{vertex_color}{Vertex fill color (default="#8DA0CB" - light purple)}

\item{vertex_alpha}{Vertex transparency (0=transparent, 1=opaque, default=0.7)}

\item{edge_color}{Edge color (default="#666666" - dark gray)}

\item{edge_alpha}{Edge transparency (0=transparent, 1=opaque, default=0.5)}

\item{main_title}{Main plot title (default="KEGG Pathway Network")}

\item{sub_title}{Subtitle text (default=NULL auto-generates based on min_shared)}

\item{legend_position}{Legend position ("bottomright", "topleft", etc., default="bottomright")}

\item{show_legend}{Show explanatory legend? (default=TRUE)}

\item{margin}{Plot margins as vector c(bottom, left, top, right) (default=c(0,0,0,0))}

\item{legend_inset}{Legend inset from plot edges c(x, y) (default=c(0.02, 0.02))}
}
\value{
Invisible igraph network object containing:
\itemize{
  \item Vertices with attributes: name, label, size
  \item Edges with attributes: weight (shared compounds), width
}
}
\description{
Creates a network visualization of metabolic pathways based on shared compounds.
Pathways are represented as nodes, and shared compounds as edges between them.
}
\details{
The network is constructed where:
\itemize{
  \item \strong{Nodes (vertices)} represent KEGG pathways
  \item \strong{Node size} is proportional to the number of compounds in the pathway
  \item \strong{Edges} connect pathways that share compounds
  \item \strong{Edge width} is proportional to the number of shared compounds
}

The function automatically calculates compound sharing between all pathway pairs.
Edges with fewer shared compounds than \code{min_shared} are filtered out.
}
\examples{
\dontrun{
# Load human pathway data
pathway_data <- getkd(species = "hsa", cache_dir = "kegg_cache")

# Basic network plot
pathwaynetwork(pathway_data)

# Customized plot with circular layout
pathwaynetwork(
  pathway_data,
  min_shared = 3,
  layout = "circle",
  vertex_color = "#66C2A5",
  edge_color = "#FC8D62",
  main_title = "Human Metabolic Pathway Network"
)

# High-resolution plot with more space for labels
pathwaynetwork(
  pathway_data,
  min_shared = 2,
  layout = "kk",
  label_size = 1.2,
  margin = c(0.1, 0.1, 0.1, 0.1),
  legend_inset = c(0.1, 0.1)
)
}
}
\seealso{
\code{\link{getkpd}} to obtain pathway data, 
         \code{\link[igraph]{plot.igraph}} for advanced network customization
}
\author{

}
\keyword{network}
\keyword{visualization}

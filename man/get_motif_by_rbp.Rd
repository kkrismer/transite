% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{get_motif_by_rbp}
\alias{get_motif_by_rbp}
\title{Retrieve motif objects by gene symbol}
\usage{
get_motif_by_rbp(rbp)
}
\arguments{
\item{rbp}{character vector of gene symbols of RNA-binding proteins}
}
\value{
A list of objects of class \code{RBPMotif}
}
\description{
Retrieves one or more motif objects identified by gene symbol.
}
\examples{
get_motif_by_rbp("ELAVL1")

get_motif_by_rbp(c("ELAVL1", "ELAVL2"))
}
\seealso{
Other motif functions: 
\code{\link{generate_iupac_by_kmers}()},
\code{\link{generate_iupac_by_matrix}()},
\code{\link{generate_kmers_from_iupac}()},
\code{\link{get_motif_by_id}()},
\code{\link{get_motifs}()},
\code{\link{get_motifs_meta_info}()},
\code{\link{get_ppm}()},
\code{\link{init_iupac_lookup_table}()},
\code{\link{set_motifs}()}
}
\concept{motif functions}

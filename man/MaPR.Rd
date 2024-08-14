\name{MaPR}
\alias{MaPR}
\title{
Master Protein abundance Regulator Prediction
}
\description{
This function predict master protein abundance regulator based on the permutation of post-transcriptional network.
}
\usage{
MaPR(Network, adjust = "BH", adjustCutoff = 0.01, rand = 1000)
}
\arguments{
  \item{Network}{
  A dataframe storing post-transcriptional regulatory network information derived from \code{\link{Network}} function.
}
  \item{adjust}{
  Correction method (default "BH") in one of \code{\link{p.adjust.methods}} for the p value of network edge.
}
  \item{adjustCutoff}{
  Correction cutoff for the definition of significance.
}
  \item{rand}{
  An integer to specify how many times the network will be permuted when predicting master protein abundance regulator.
}
}
\details{
This function is required to load in network generated from \code{\link{Network}} function and will output one file named NetworkRegulatorDegreeSignicance.txt, providing the permutation significance of a gene predicted as master protein abundance regulator.
}
\examples{
# Load in post-transcriptional regulatory network constructed by Network function
infile <- gzfile('NetWork.txt.gz','r')
Network <- read.table(infile,header = T,sep = '\t',stringsAsFactors = F)
close(infile)
# Predict master protein abundance regulator based on permutation of post-transcriptional regulatory network
MaPR(Network = Network, adjust = 'BH', adjustCutoff = 0.01, rand = 1000)
}
\keyword{MaPR}

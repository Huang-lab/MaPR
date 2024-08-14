\name{Network}
\alias{Network}
\title{
Post-transcriptional regulatory network construction
}
\description{
This function constructs post-transcriptional regulatory network by integrating sample-paired mRNA profile and protome profile.
}
\usage{
Network(mRNAProfile, proteinProfile, 
        NetworkType = "ProteinSpecific", OutputFileName = "NetWork.txt", 
        method = "spearman", adjust = "BH", adjustCutoff = 0.01, 
        mRNAProfileNARatio = NULL, proteinProfileNARatio = NULL, 
        mRNAProfileZeroRatio = NULL, proteinProfileZeroRatio = NULL, 
        ProfileNARatio = 0.3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mRNAProfile}{
a numeric matrix representing expression profile of mRNA whose rowname is gene name.
}
  \item{proteinProfile}{
a numeric matrix representing expression profile of protein whose rowname is gene name and column number must be equal to \emph{mRNAProfile}. 
}
  \item{NetworkType}{
a character (default "ProteinSpecific") to specify which type of network will be constructed.
One of "ProteinSpecific" (default), "mRNASpecific", or "BothSig", can be abbreviated.
}
  \item{OutputFileName}{
a character to specify output file name of translation network
}
  \item{method}{
a character string (default "spearman") indicating which correlation coefficient is to be computed.
One of "pearson" (default), "kendall", or "spearman", can be abbreviated.
}
  \item{adjust}{
Correction method (default "BH") in one of \code{\link{p.adjust.methods}}.
}
  \item{adjustCutoff}{
Correction cutoff for the definition of significance.
}
  \item{mRNAProfileNARatio}{
a numeric (rangs from 0 to 1,default 0.3) specifying proportion of NA value of mRNA profile, which would be removed from further analysis.
}
  \item{proteinProfileNARatio}{
a numeric (rangs from 0 to 1,default 0.3) specifying proportion of NA value of protein profile, which would be removed from further analysis.
}
  \item{mRNAProfileZeroRatio}{
a numeric (rangs from 0 to 1,default 0.3) specifying proportion of zero value of mRNA profile, which would be removed from further analysis.
}
  \item{proteinProfileZeroRatio}{
a numeric (rangs from 0 to 1,default 0.3) specifying proportion of zero value of protein files, which would be removed from further analysis.
}
  \item{ProfileNARatio}{
a numeric (rangs from 0 to 1,default 0.3) specifying proportion of NA value of genes during semi-partial correlation process, which would be removed from further analysis.
}
}
\details{
Note:All the arguments without default value must be assigned.

This function would output five files as following:
   mRNACoexpCorNetWork.txt.gz stores mRNA expression correlation matrix across all the genes
   mRNACoexpPvalFDRNetWork.txt.gz stores mRNA expression correlation significance with ajusted p value above the diagonal.
   ProteinCoexpCorNetWork.txt.gz stores protein expression correlation matrix across all the genes
   ProteinCoexpPvalFDRNetWork.txt.gz stores protein expression correlation significance with ajusted p value above the diagonal.
   NetWork.txt.gz stores all the candidate translational regulations information of four columns: Regulator, Target, p.value, estimate.
}
\examples{
##Translation network construction
Network(mRNAProfile = datatests$mRNAProfile, proteinProfile = datatests$proteinProfile,
        NetworkType = 'ProteinSpecific', OutputFileName = 'NetWork.txt.gz',
        method = 'spearman', adjust = 'BH', adjustCutoff = 0.01, 
        mRNAProfileNARatio = NULL, proteinProfileNARatio = 0.3,
        mRNAProfileZeroRatio = 0.3, proteinProfileZeroRatio = NULL,
        ProfileNARatio = 0.3)
}
\keyword{Network}

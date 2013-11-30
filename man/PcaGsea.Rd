\name{PcaGsea}
\alias{PcaGsea}
\title{PcaGsea}
\description{find loadings of a numeric matrix and do GSEA on each loadings vector}
\usage{PcaGsea(exp, geneset = NULL, n = 20, cumulative.proportion = NULL, 
    proportion.variance = NULL, ...)}
\arguments{
  \item{exp}{pxn expression matrix (genes in rows, samples in columns)}
  \item{geneset}{a list where each element contains a vector of genes, and the list
names are the names of the genesets}
  \item{n}{the number of principal components to do the calculations on (n=10 suggested).}
  \item{cumulative.proportion}{Keep enough principal components to account for this total proportion of the variance.  Over-rides n.}
  \item{proportion.variance}{Keep any principal components which account for at least this
proportion of the variance.  Over-rides n and
cumulative.proportion.}
  \item{\dots}{arguments passed on to prcomp function.  scale.=TRUE is advisable.}
}

\value{a list with four elements. 1: the loadings of the expression matrix
for the components kept, 2: the scores, 3: the unadjusted p-values
from GSEA, and 4: the prcomp object which can be used to make a
screeplot, calculate scores from new data, etc.}

\author{Princy Parsana and Levi Waldron}






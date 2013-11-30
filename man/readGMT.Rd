\name{readGMT}
\alias{readGMT}
\title{Read a gmt (gene matrix transposed) file}
\description{Reads a gmt (gene matrix transposed) file as provided by mSigDB or
geneSigDB, to an R list object.}
\usage{readGMT(fname, name.column = 1)}
\arguments{
  \item{fname}{fully specified path to a .gmt file, ie from mSigDB or geneSigDB}
  \item{name.column}{name.column = which column to use for names (typically 1 or 2).}
}

\value{a list object where each element is a string of genes, and the
elements are named by the first column of the .gmt file.}

\author{Princy Parsana and Levi Waldron}






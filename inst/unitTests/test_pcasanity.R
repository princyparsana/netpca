test_pcasanity <- function(){
    library(pensim)
    library(netpca)
    set.seed(1)
    ##Five variables with zero covariance, five with block covariance 0.1 to 0.5:
    x <- create.data(nvars=rep(100, 10),
                     cors=c(0.5, 0.4, 0.3, 0.2, 0.1, 0, 0, 0, 0, 0),
                     associations=rep(0, 10),
                     firstonly=rep(TRUE, 10),
                     nsamples=500,
                     response="binary",logisticintercept=0.5)
    exp.mat <- t( x$data[, -ncol(x$data)] )
    ## sample each variable type to create gene sets
    set.seed(1)
    genesets <- lapply(letters[1:10], function(i){
        sample(grep(i, rownames(exp.mat), val=TRUE), 25)
    })
    names(genesets) <- letters[1:10]
    ##Do PCA
    pca.out <- PcaGsea(exp.mat,genesets,10)
    ##loadings 1-5 should be strongly enriched for genesets 1-5, and nothing else:
    for (i in 1:5){
        print(i)
        i.found <- which.min(pca.out$raw.pval[, i])
        names(i.found) <- NULL
        checkEquals( i.found, i )
        checkTrue( min( pca.out$raw.pval[, i] ) < 1e-10 )
        checkTrue( min( pca.out$raw.pval[-i, i] ) > 1e-10 )
    }
    ##loadings 6-10 have little enrichment for anything:
    for (i in 6:10){
        print(i)
        checkTrue( min( pca.out$raw.pval[, i] ) > 1e-10 )
    }
}

PcaGsea <- function
###find loadings of a numeric matrix and do GSEA on each loadings vector
(exp,
###pxn expression matrix (genes in rows, samples in columns)
 geneset=NULL,
###a list where each element contains a vector of genes, and the list
###names are the names of the genesets
 n=20,
###the number of principal components to do the calculations on (n=10 suggested).
 cumulative.proportion=NULL,
###Keep enough principal components to account for this total proportion of the variance.  Over-rides n.
 proportion.variance=NULL,
###Keep any principal components which account for at least this
###proportion of the variance.  Over-rides n and
###cumulative.proportion.
 ...
### arguments passed on to prcomp function.  scale.=TRUE is advisable.
 ){
    if(nrow(exp) < ncol(exp))
        warning("exp should probably have more rows (features) than columns (samples).")
    exp <- as.matrix(exp)
    exp.pca <- prcomp(t(exp), ...)
    prop.of.variance <- summary(exp.pca)$importance["Proportion of Variance", ]
    cum.proportion <- summary(exp.pca)$importance["Cumulative Proportion", ]
    if( !is.null(proportion.variance) ){
        n <- sum(prop.of.variance > proportion.variance)
        print(paste("Using n = ", n, " components explaining at least ", 100*proportion.variance,
                    "% of total variance.", sep=""))
    }else if( !is.null(cumulative.proportion) ){
        n <- sum(cum.proportion < cumulative.proportion)
        print(paste("Using n = ", n, " components which cumulatively account for ", 100*cumulative.proportion,
                    "% of total variance.", sep=""))
    }
    exp.loadings <- exp.pca$rotation[, c(1:n)]
    exp.scores <- predict(exp.pca)[, c(1:n)]
    if(is.null(geneset)){
        gsea <- NULL
    }else{
        match.geneset <- vector("list", length=length(geneset))
        names(match.geneset) <- names(geneset)
        for(i in 1:length(geneset)){
            geneset.genes <- intersect(geneset[[i]], rownames(exp.loadings))
            match.geneset[[i]] <- geneset.genes
        }
        gsea <- matrix(NA, nrow=length(match.geneset), ncol=n)
        rownames(gsea) <- names(match.geneset)
        colnames(gsea) <- colnames(exp.loadings)
        for(k in 1:n){
            gene.loadings <- exp.loadings[, k]
            names(gene.loadings) <- rownames(exp.loadings)
            rank.loadings <- order(gene.loadings, decreasing=TRUE) 
            gene.loadings <- gene.loadings[rank.loadings]
            ranked.genes <- names(gene.loadings) #Rank the genes depending on loadings
            for(j in 1:length(match.geneset)){
                geneset.pvalue <- limma::wilcoxGST(match.geneset[[j]], gene.loadings) 
                gsea[j, k] <- geneset.pvalue #p-value
            }
        }
    }
    pc.gsea <- list(exp.loadings, exp.scores, gsea, exp.pca)
    names(pc.gsea) <- c("loadings", "scores", "raw.pval", "prcomp.obj")
    return(pc.gsea)
###a list with four elements. 1: the loadings of the expression matrix
###for the components kept, 2: the scores, 3: the unadjusted p-values
###from GSEA, and 4: the prcomp object which can be used to make a
###screeplot, calculate scores from new data, etc.
}

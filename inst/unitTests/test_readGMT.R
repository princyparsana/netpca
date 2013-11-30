testGMT <- function(){
    biocarta <- readGMT(system.file("extdata/c2.cp.biocarta.v3.1.entrez.gmt", package = "netpca"))
    checkEquals(length(biocarta), 217)
    checkEquals(class(biocarta), "list")
}

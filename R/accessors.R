get_mles <- function(fits){
    if(class(fits) == "list"){
        min.idx <- map(fits, 'nll') %>% map(.f=function(X) which(X == min(X)))
        mles <- lapply(1:length(fits), function(X) fits[[X]][min.idx[[X]],])
        names(mles) <- names(fits)

        labels <- map(mles, rownames)
        specTrees <- c(sum(labels == "ABC"), sum(labels == "BAC"),
                       sum(labels == "BCA"), sum(labels == "CBA"),
                       sum(labels == "ACB"), sum(labels == "CAB"))
        names(specTrees) <- c("ABC", "BAC", "BCA", "CBA", "ACB", "CAB")
        mles$speciesTrees <- specTrees
    } else if(class(fits) == "data.frame"){
        min.idx <- which(fits$nll == min(fits$nll))
        mles <- fits[min.idx,]

        labels <- rownames(mles)
        specTrees <- c(sum(labels == "ABC"), sum(labels == "BAC"),
                       sum(labels == "BCA"), sum(labels == "CBA"),
                       sum(labels == "ACB"), sum(labels == "CAB"))
        names(specTrees) <- c("ABC", "BAC", "BCA", "CBA", "ACB", "CAB")
        mles <- list(mles, specTrees)
    }
    mles
}

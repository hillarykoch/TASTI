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


get_triples <- function(file, branchLengths = TRUE){
    trees <- read_table(file, col_names = FALSE)
    taxa <- read.tree(text = trees$X1[1])$tip.label
    splits <- str_split(trees$X1, pattern = "")

    input <- matrix(rep(0, 5*nrow(trees)), ncol = 5)
    colnames(input) <- c(taxa, "t1", "t2")
    for(i in 1:nrow(trees)){
        test <- splits[[i]]
        if(all(test[1:2] == c("(", "("))){
            commaIdx <- grep(test, pattern = ",")
            colonIdx <- grep(test, pattern = ":")

            outgroup <- paste0(test[(commaIdx[2]+1):(colonIdx[4]-1)], collapse = "")
            t1 <- as.numeric(paste0(test[(colonIdx[1]+1):(commaIdx[1]-1)], collapse = ""))
            t2 <- as.numeric(paste0(test[(colonIdx[3]+1):(commaIdx[2]-1)], collapse = ""))
        } else{
            commaIdx <- grep(test, pattern = ",")
            colonIdx <- grep(test, pattern = ":")
            parenIdx <- grep(test, pattern = "\\(")

            outgroup <- paste0(test[(parenIdx[1]+1):(colonIdx[1]-1)], collapse = "")
            t1 <- as.numeric(paste0(test[(colonIdx[2]+1):(commaIdx[2]-1)], collapse = ""))
            t2 <- as.numeric(paste0(test[(colonIdx[1]+1):(commaIdx[1]-1)], collapse = ""))
        }

        input[i,outgroup] <- 1
        input[i,"t1"] <- t1
        input[i,"t2"] <- t2
    }

    if(branchLengths){
        input
    } else{
        colSums(input)[1:3]
    }
}


parse.nucmer.coords=function(dfx)
{
    names(dfx) = c("start1", "end1", "start2", "end2", "length1", "length2", "identity", "name1", "name2")

    dfx$is.ref1 = !grepl("\\|", dfx$name1)
    dfx$is.ref2 = !grepl("\\|", dfx$name2)

    dfx$sample1 = ifelse(dfx$is.ref1, "ref", sapply(strsplit(dfx$name1, "\\|"), function(x) { x[1] }))
    dfx$cycle1 = ifelse(dfx$is.ref1, dfx$name1, sapply(strsplit(dfx$name1, "\\|"), function(x) { x[2] }))
    dfx$sample2 = ifelse(dfx$is.ref2, "ref", sapply(strsplit(dfx$name2, "\\|"), function(x) { x[1] }))
    dfx$cycle2 = ifelse(dfx$is.ref2, dfx$name2, sapply(strsplit(dfx$name2, "\\|"), function(x) { x[2] }))
    dfx = dfx[,-match(c("name1", "name2"), names(dfx))]
    dfx$key1 = paste0(dfx$sample1, "_", dfx$cycle1)
    dfx$key2 = paste0(dfx$sample2, "_", dfx$cycle2)
    dfx
}

parse.nucmer.snps=function(dfs)
{
    if (!is.null(dfs)) {
        dfs = dfs[,c(1,4,2,3,11,12)]
        names(dfs) = c("base1", "base2", "nt1", "nt2", "name1", "name2")

        dfs$is.ref1 = !grepl("\\|", dfs$name1)
        dfs$is.ref2 = !grepl("\\|", dfs$name2)

        dfs$sample1 = ifelse(dfs$is.ref1, "ref", sapply(strsplit(dfs$name1, "\\|"), function(x) { x[1] }))
        dfs$cycle1 = ifelse(dfs$is.ref1, dfs$name1, sapply(strsplit(dfs$name1, "\\|"), function(x) { x[2] }))
        dfs$sample2 = ifelse(dfs$is.ref2, "ref", sapply(strsplit(dfs$name2, "\\|"), function(x) { x[1] }))
        dfs$cycle2 = ifelse(dfs$is.ref2, dfs$name2, sapply(strsplit(dfs$name2, "\\|"), function(x) { x[2] }))
        dfs$key1 = paste0(dfs$sample1, "_", dfs$cycle1)
        dfs$key2 = paste0(dfs$sample2, "_", dfs$cycle2)

        ## dfs$sample1 = sapply(strsplit(dfs$name1, "\\|"), function(x) { x[1] })
        ## dfs$cycle1 = sapply(strsplit(dfs$name1, "\\|"), function(x) { x[2] })
        ## dfs$sample2 = sapply(strsplit(dfs$name2, "\\|"), function(x) { x[1] })
        ## dfs$cycle2 = sapply(strsplit(dfs$name2, "\\|"), function(x) { x[2] })
        ## dfs$key1 = paste0(dfs$sample1, "_", dfs$cycle1)
        ## dfs$key2 = paste0(dfs$sample2, "_", dfs$cycle2)
    } else {
        dfs = setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("base1", "base2", "nt1", "nt2", "name1", "name2", "key1", "key2"))
    }
    dfs
}

parse.nucmer.coords=function(dfx)
{
    names(dfx) = c("start1", "end1", "start2", "end2", "alength1", "alength2", "identity", "name1", "name2")

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

align.summary=function(dfc, dfx, dfs)
{
    keys = dfc$key
    N = length(keys)
    result.table = NULL
    result.mat = matrix(rep(0, N^2), N, N)
    rownames(result.mat) = keys
    colnames(result.mat) = keys

    dfx = dfx[dfx$key1 != dfx$key2 & is.element(dfx$key1, keys) & is.element(dfx$key2, keys),]
    dfx$key = paste(dfx$key1, dfx$key2, sep=":")
    for (key in unique(dfx$key)) {
        dfx.ii = dfx[dfx$key == key,]
        dfx.ii$identity = dfx.ii$identity / 100

        key1 = dfx.ii$key1[1]
        key2 = dfx.ii$key2[1]

        length1 = dfc$length[match(key1, dfc$key)]
        length2 = dfc$length[match(key2, dfc$key)]

        get.total.alignment=function(start, end, len)
        {
            left = pmin(start, end)
            right = pmax(start, end)
            vv = rep(0, len)
            for (i in 1:length(left)) {
                vv[left[i]:right[i]] = 1
            }
            sum(vv)
        }

        # weights are by match length
        dfx.ii$match.length1 = pmin(length1, abs(dfx.ii$end1 - dfx.ii$start1) + 1)
        dfx.ii$match.length2 = pmin(length2, abs(dfx.ii$end2 - dfx.ii$start2) + 1)
        dfx.ii$weight1 = dfx.ii$match.length1 / length1
        dfx.ii$weight2 = dfx.ii$match.length2 / length2
        dfx.ii$weight = (dfx.ii$weight1 + dfx.ii$weight2) / 2
        dfx.ii$weight = dfx.ii$weight / sum(dfx.ii$weight)

        # get alignment length
        total.alignment1 = get.total.alignment(start=dfx.ii$start1, end=dfx.ii$end1, len=length1)
        total.alignment2 = get.total.alignment(start=dfx.ii$start2, end=dfx.ii$end2, len=length2)
        overlap.identity = 100 * sum(dfx.ii$identity * dfx.ii$weight)

        overlap.fraction = 100*(total.alignment1/length1 + total.alignment2/length2)/2

        # dfx
        ix = dfs$key1 == key1 & dfs$key2 == key2
        total.snps = sum(ix)

        df = data.frame(key1=key1, key2=key2, length1=length1, length2=length2,
                        overlap.fraction=overlap.fraction,
                        overlap.identity=overlap.identity,
                        snp.count=total.snps)
        result.table = rbind(result.table, df)

        i1 = match(key1, keys)
        i2 = match(key2, keys)
        result.mat[i1, i2] = overlap.identity
    }
    result.table$dist = 1 - result.table$overlap.fraction / 100 * result.table$overlap.identity / 100
    list(table=result.table, mat=result.mat)
}

order.cycles=function(mat)
{
    hh = hclust(as.dist(100-mat))
    colnames(mat)[hh$order]
}

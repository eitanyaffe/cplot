options(stringsAsFactors=F)

align.summary=function(dfc, dfx, dfs)
{
    keys = dfc$key
    N = length(keys)
    result.table = NULL
    result.mat = matrix(rep(0, N^2), N, N)
    rownames(result.mat) = keys
    colnames(result.mat) = keys
    for (i1 in 1:N) {
        for (i2 in 1:N) {
            key1 = keys[i1]
            key2 = keys[i2]

            # dfx
            ix = dfx$key1 == key1 & dfx$key2 == key2
            if (sum(ix) == 0)
                stop("no alignment found")
            dfx.ii = dfx[ix,]
            dfx.ii$match.length1 = abs(dfx.ii$end1 - dfx.ii$start1) + 1
            dfx.ii$match.length2 = abs(dfx.ii$end2 - dfx.ii$start2) + 1
            dfx.ii$weight1 = dfx.ii$match.length1 / dfx.ii$length1
            dfx.ii$weight2 = dfx.ii$match.length2 / dfx.ii$length2
            dfx.ii$weight = (dfx.ii$weight1 + dfx.ii$weight2) / 2
            dfx.ii$weight = dfx.ii$weight / sum(dfx.ii$weight)

            length1 = dfx.ii$length1[1]
            length2 = dfx.ii$length2[1]

            total.alignment1 = sum(dfx.ii$match.length1)
            total.alignment2 = sum(dfx.ii$match.length2)
            total.alignment = round((total.alignment1 + total.alignment2) / 2)
            total.identity = sum(dfx.ii$identity * dfx.ii$weight)

            total.alignment.percent = 100*(total.alignment1/length1 + total.alignment2/length2)/2

            # dfx
            ix = dfs$key1 == key1 & dfs$key2 == key2
            total.snps = sum(ix)

            df = data.frame(key1=key1, key2=key2, length1=length1, length2=length2,
                            match.length=total.alignment, match.percent=total.alignment.percent,
                            identity=total.identity, snp.count=total.snps)
            result.table = rbind(result.table, df)

            result.mat[i1, i2] = total.identity
        }
    }
    list(table=result.table, mat=result.mat)
}

order.cycles=function(mat)
{
    hh = hclust(as.dist(100-mat))
    colnames(mat)[hh$order]
}

plot.nucmer.cluster=function(dfc, set.id, set.title, cluster)
{
    idir.cluster = "/relman01/home/nshalon/work/pipe/sour/cluster"
    idir.nucmer = "/relman01/home/nshalon/work/pipe/sour/nucmer"

    dfx = read.cache(paste0(idir.nucmer, "/", set.id, "/CLUSTER_", cluster, ".coords"), function(fn) { read.delim(fn, header=F) })
    dfs = read.cache(paste0(idir.nucmer, "/", set.id, "/CLUSTER_", cluster, ".snps"), function(fn) { safe.read.delim(fn, header=F) })

    dfc = dfc[dfc$cluster == cluster,]
    dfc$key = paste0(dfc$sample, "_", dfc$cycle)

    # are these columns correct?
    names(dfx) = c("start1", "end1", "start2", "end2", "length1", "length2", "identity", "name1", "name2")
    dfx$sample1 = sapply(strsplit(dfx$name1, "\\|"), function(x) { x[1] })
    dfx$cycle1 = sapply(strsplit(dfx$name1, "\\|"), function(x) { x[2] })
    dfx$sample2 = sapply(strsplit(dfx$name2, "\\|"), function(x) { x[1] })
    dfx$cycle2 = sapply(strsplit(dfx$name2, "\\|"), function(x) { x[2] })
    dfx = dfx[,-match(c("name1", "name2", "length1", "length2"), names(dfx))]
    dfx$key1 = paste0(dfx$sample1, "_", dfx$cycle1)
    dfx$key2 = paste0(dfx$sample2, "_", dfx$cycle2)

    if (!is.null(dfs)) {
        dfs = dfs[,c(1,4,2,3,11,12)]
        names(dfs) = c("base1", "base2", "nt1", "nt2", "name1", "name2")
        dfs$sample1 = sapply(strsplit(dfs$name1, "\\|"), function(x) { x[1] })
        dfs$cycle1 = sapply(strsplit(dfs$name1, "\\|"), function(x) { x[2] })
        dfs$sample2 = sapply(strsplit(dfs$name2, "\\|"), function(x) { x[1] })
        dfs$cycle2 = sapply(strsplit(dfs$name2, "\\|"), function(x) { x[2] })
        dfs$key1 = paste0(dfs$sample1, "_", dfs$cycle1)
        dfs$key2 = paste0(dfs$sample2, "_", dfs$cycle2)
    } else {
        dfs = setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("base1", "base2", "nt1", "nt2", "name1", "name2", "key1", "key2"))
    }

    dfx = dfx[is.element(dfx$key1, dfc$key) & is.element(dfx$key2, dfc$key),]
    dfs = dfs[is.element(dfs$key1, dfc$key) & is.element(dfs$key2, dfc$key),]
    dfx$length1 = dfc$length[match(dfx$key1, dfc$key)]
    dfx$length2 = dfc$length[match(dfx$key2, dfc$key)]

    dll = align.summary(dfc=dfc, dfx=dfx, dfs=dfs)
    df = dll$table
    keys = order.cycles(mat=dll$mat)
    N = length(keys)

    odir = paste0("figures/nucmer/", set.id)
    system(paste("mkdir -p", odir))
    ofn = paste0(odir, "/", cluster, ".pdf")
    cat(sprintf("generating plot: %s\n", ofn))

    width = N + 0.75
    height = N + 0.75
    pdf(ofn, width=width, height=height)
    layout(matrix(1:(N*N), N, N))
    par(xaxs="i", yaxs="i")
    gap = 0.05
    for (i1 in 1:N) {
        for (i2 in 1:N) {
            key1 = keys[i1]
            key2 = keys[i2]
            df.ii = df[df$key1 == key1 & df$key2 == key2,]
            dfx.ii = dfx[dfx$key1 == key1 & dfx$key2 == key2,]
            dfs.ii = dfs[dfs$key1 == key1 & dfs$key2 == key2,]
            # diagonal
            if (i1 == i2) {
                par(mai=c(gap,gap,gap,gap))
                plot.new()
                plot.window(xlim=0:1, ylim=0:1)
                ix = match(key1, dfc$key)
                text(x=0.5, y=0.5, paste0(dfc$sample[ix], "\n", dfc$cycle[ix], "\n", df.ii$length1, "bp"), font=2, cex=1.25)
                box()
            } else if (i1 > i2) {
                par(mai=c(gap,gap,gap,gap))
                plot.new()
                xlim = c(0, df.ii$length1)
                ylim = c(0, df.ii$length2)
                plot.window(xlim=xlim, ylim=ylim)
                grid()
                points(x=dfs.ii$base1, y=dfs.ii$base2, cex=0.5, pch=3, col="red")
                segments(x0=dfx.ii$start1, x1=dfx.ii$end1, y0=dfx.ii$start2, y1=dfx.ii$end2)
                box()
            } else {
                par(mai=c(gap,gap,gap,gap))
                plot.new()
                plot.window(xlim=0:1, ylim=0:1)
                lab = paste0(df.ii$match.length, "bp (", round(df.ii$match.percent,3), "%)",
                             "\nD=", round(100-df.ii$identity,3), "%",
                             "\nsnps=", df.ii$snp.count)
                text(x=0.5, y=0.5, lab)
                box()
            }
        }
    }
    graphics.off()
}

append.cycle.data=function(dfc)
{
    dfc$key = paste0(dfc$sample, "_", dfc$cycle)
    samples = unique(dfc$sample)
    rr = NULL
    for (sample in samples) {
        idir = paste0("/oak/stanford/groups/relman/users/nshalon/pipe/cyc_find_out4/", sample, "_0.0")
        df = read.cache(paste0(idir, "/dominant_cycles_filter/cycle_summary_classification"))
        df$key  = paste0(sample, "_", df$cycle)
        rr = rbind(rr, df)
    }
    ix = match(dfc$key, rr$key)
    if (any(is.na(ix)))
        stop("not all cycles found")
    dfc$length = rr$length[ix]
    dfc
}

plot.nucmer=function(set.title="uniq_final_k21_cycles")
{
    idir.cluster = "/relman01/home/nshalon/work/pipe/sour/cluster"

    set.id = set.title
    dfc = read.cache(paste0(idir.cluster, "/", set.id))
    dfc = append.cycle.data(dfc)
    dfc = dfc[dfc$length >= 1000,]
    tt = sort(table(dfc$cluster), decreasing=T)
    tt = tt[tt>1]
    cat(sprintf("dataset table: %s\n", paste0(idir.cluster, "/", set.id)))
    cat(sprintf("number of clusters: %d\n", length(tt)))
    for (cluster in as.numeric(names(tt)))
        plot.nucmer.cluster(dfc=dfc, cluster=cluster, set.id=set.id, set.title=set.title)
}

plot.nucmer.all=function()
{
    title.ids = "uniq_final_k21_cycles"
    for (set.title in title.ids)
        plot.nucmer(set.title=set.title)
}

rl.align=function()
{
    source("cplot.r")
    source("palign.r")
}

ex.align=function()
{
    plot.nucmer.all()
}

options(stringsAsFactors=F)

.ver = "v2"

######################################################################################################
# util functions
######################################################################################################

get.cluster.rects=function(ct, hc)
{
    vv = ct[hc$order]
    N = length(ct)
    start = which(diff(c(-1,vv)) != 0) - 1
    endx = start
    data.frame(start=start, end=c(endx[-1], N))
}

get.min.dist=function(df, dfm)
{
    N = dim(df)[1]
    dfm = dfm[is.element(dfm$key1, df$key) &  is.element(dfm$key2, df$key),]
    dfm$i1 = match(dfm$key1, df$key)
    dfm$i2 = match(dfm$key2, df$key)
    mm = smatrix2matrix.sym(smatrix=dfm, dim=N, i.field="i1", j.field="i2", value.field="dist", default.value=1)
    diag(mm) = 1
    apply(mm, 1, min)
}

smatrix2matrix.sym=function(smatrix, dim, i.field="i", j.field="j", value.field="value", default.value=0)
{
  indices1 = smatrix[,i.field] + (smatrix[,j.field]-1) * dim
  indices2 = smatrix[,j.field] + (smatrix[,i.field]-1) * dim
  v = rep(default.value, dim*dim)
  v[indices1] = smatrix[,value.field]
  v[indices2] = smatrix[,value.field]
  matrix(v, dim, dim)
}

######################################################################################################
# generate basic files for plotting mat/dendrogram
######################################################################################################

cluster.cycs=function()
{
    dfm = read.delim("output/nucmer_summary.txt")
    df = read.delim("output/cycles_basic.txt")

    df$min.dist = get.min.dist(df=df, dfm=dfm)
    #df = df[df$min.dist <= 0.5,]
    N = dim(df)[1]

    dfm = dfm[is.element(dfm$key1, df$key) &  is.element(dfm$key2, df$key),]
    dfm$i1 = match(dfm$key1, df$key)
    dfm$i2 = match(dfm$key2, df$key)

    mm = smatrix2matrix.sym(smatrix=dfm, dim=N, i.field="i1", j.field="i2", value.field="dist", default.value=1)
    diag(mm) = 0
    hc = hclust(as.dist(mm), method="single")
    ct = cutree(hc, h=0.05)
    df$group = paste0("g", ct)
    df = df[hc$order,]

    # add prev cluster
    dfo = read.delim("data/alignment_stats.txt")
    dfo = dfo[,1:2]
    dfo$key = sapply(strsplit(dfo$order, "\\|"), function(x) { paste(x[1],x[2],sep="_") })
    ix = match(df$key, dfo$key)
    df$prev.cluster = ifelse(!is.na(ix), dfo$species[ix], "na")

    phix.group = df$group[match("gut_q_CYC115", df$key)]
    df$group[df$group == phix.group] = "PhiX"

    tt = table(df$group)
    tt = tt[tt>1]
    rr = data.frame(group=names(tt), count=as.numeric(tt))
    rr$prev.cluster = ""
    for (i in 1:dim(rr)[1]) {
        keys.i = df$key[df$group == rr$group[i]]
        if (length(keys.i) > 1) {
            dfo.i = dfo[is.element(dfo$key, keys.i),]
            prev.clusters = unique(dfo.i$species)
            if (length(prev.clusters) != 1)
                stop("expecting 1:1 mapping")
            rr$prev.cluster[i] = prev.clusters
        } else {
            rr$prev.cluster[i] = "na"
        }
    }
    rr = rr[rr$group != "PhiX",]
    rr = rr[order(rr$count, decreasing=T),]
    N = dim(rr)[1]
    rr$id = paste0("M", 1:N)

    df = df[is.element(df$group, rr$group),]
    df$id = rr$id[match(df$group, rr$group)]

    odir = paste0("output/clusters_", .ver)
    system(paste("mkdir -p", odir))
    cat(sprintf("saving clusters to: %s\n", odir))
    write.table(x=df, file=paste0(odir, "/cycles_grouped.txt"), quote=F, row.names=F, sep="\t")
    write.table(x=rr, file=paste0(odir, "/groups.txt"), quote=F, row.names=F, sep="\t")
}

######################################################################################################
# plot matrix/dendro
######################################################################################################

plot.clustering=function()
{
    dfm = read.delim("output/nucmer_summary.txt")
    df = read.delim("output/cycles_basic.txt")
    dfc = read.delim(paste0("output/clusters_", .ver, "/cycles_grouped.txt"))

    df$sample = label.sample(df$sample)
    df$id = dfc$id[match(df$key, dfc$key)]

    # add environment
    df$env = sapply(strsplit(df$sample, "_"), function(x) { x[1] })
    df$env[df$env == "pilot"] = "gut"

    # add cluster
    dfo = read.delim("data/alignment_stats.txt")
    dfo = dfo[,1:2]
    dfo$key = sapply(strsplit(dfo$order, "\\|"), function(x) { paste(x[1],x[2],sep="_") })
    ix = match(df$key, dfo$key)
    df$cluster = ifelse(!is.na(ix), dfo$species[ix], "na")

    # cycle index
    df$cycle.index = sapply(df$cycle, function(x) { gsub("CYC", "", x[1]) })

    df$min.dist = get.min.dist(df=df, dfm=dfm)
    df = df[df$min.dist <= 0.5,]
    N = dim(df)[1]

    df$label.detailed = paste(df$id, df$cluster, paste0(df$sample, " : c", df$cycle.index))
    df$label.clean = paste(paste0(df$sample, " : c", df$cycle.index))

    df$color = ifelse(df$env == "gut", "blue", ifelse(df$env == "sewage", "green", "orange"))
    dfm = dfm[is.element(dfm$key1, df$key) &  is.element(dfm$key2, df$key),]
    dfm$i1 = match(dfm$key1, df$key)
    dfm$i2 = match(dfm$key2, df$key)

    mm = smatrix2matrix.sym(smatrix=dfm, dim=N, i.field="i1", j.field="i2", value.field="dist", default.value=1)
    diag(mm) = 0
    rownames(mm) = df$label.clean
    colnames(mm) = df$label.clean
    hc = hclust(as.dist(mm), method="single")

    ct = cutree(hc, h=0.05)

    # keep df before reordering
    df.orig = df
    df = df[hc$order,]
    df$index = 1:N

    mmo = mm[hc$order,hc$order]
    sm = matrix2smatrix(mmo)
    lim = c(0,N)

    breaks = c(0, 0.5, 0.95, 0.9999999,1)
    cols = c("white", "blue", "darkblue", "red", "orange")
    legend.labels = c("0%", "50%", "95%", "100%")
    panel = make.color.panel(colors=cols)
    sm$col = panel[vals.to.cols(vals=1-sm$value, breaks=breaks)]

    odir = paste0("figures/mat_", .ver)
    system(paste("mkdir -p", odir))

    breaks.clean = c(0, 0.5, 0.95,1)
    cols.clean = c("white", "blue", "darkblue", "red")
    legend.labels = c("0%", "50%", "95%", "<100%")
    panel.clean = make.color.panel(colors=cols.clean)
    wlegend2(fdir=odir, breaks=breaks.clean, labels=legend.labels, panel=panel.clean, title="ANI",
             size=1, height=4, width=1, make.file=T, tick.size=0.3, border=NA, lwd=2)

    ###############################################################
    # matrix
    ###############################################################

    # get cluster rects
    bb = get.cluster.rects(ct=ct, hc=hc)

    pmat=function(clean, ofn)
    {
        cat(sprintf("plotting: %s\n", ofn))
        labels = if (clean) df$label.clean else df$label.detailed
        pdf(ofn, height=0.2*N, width=0.2*N)
        par(mai=c(2, 2, 0.5, 0.5))
        plot.new()
        plot.window(xlim=lim, ylim=lim)
        rect(xleft=sm$i-1, xright=sm$i, ybottom=sm$j-1, ytop=sm$j, col=sm$col, border=NA, lwd=0.5)
        abline(h=unique(c(bb$start, bb$end)), v=unique(c(bb$start, bb$end)), col="gray")
        rect(xleft=bb$start, xright=bb$end, ybottom=bb$start, ytop=bb$end, col=NA, border=1)
        mtext(side=1, text=labels, at=df$index-0.5, las=2, cex=0.75, line=-0.5)
        mtext(side=2, text=labels, at=df$index-0.5, las=2, cex=0.75, line=-0.5)
        graphics.off()
    }
    pmat(ofn=paste0(odir, "/matrix_clean.pdf") ,clean=T)
    pmat(ofn=paste0(odir, "/matrix_detailed.pdf") ,clean=F)

    ###############################################################
    # dendro
    ###############################################################

    pdend=function(style)
    {
        ofn =  paste0(odir, "/dendro_", style, ".pdf")
        cat(sprintf("plotting: %s\n", ofn))
        hc$labels = if (style == "clean") df.orig$label.clean else df.orig$label.detailed
        dd = as.dendrogram(hc)

        height = 0.2*N
        width = 5
        pdf(ofn, height=height, width=width)
        par(mai=c(1, 0.2, 0.5, 3.5))
        par(cex=0.7, las=2)
        plot(dd, horiz=T)
        graphics.off()
    }
    pdend(style="clean")
    pdend(style="detailed")

    colLab.clean=function(n){
        if(is.leaf(n)){
            a = attributes(n)
            ix = match(attributes(n)$label,df.orig$label.clean)
            col = df.orig$color[ix]
            attr(n,"nodePar") = c(a$nodePar,list(cex=0.65, lab.cex=0.7, pch=15, col=col, lab.col=1, lab.font=1))
        }
        return(n)
    }

    pdend.nodes=function(style)
    {
        ofn =  paste0(odir, "/dendro_nodes_", style, ".pdf")
        cat(sprintf("plotting: %s\n", ofn))
        hc$labels = df.orig$label.clean
        dd = as.dendrogram(hc)
        ddc = dendrapply(dd, colLab.clean)

        leaflab = if(style == "clean") "none" else "perpendicular"
        width = 4
        height = 0.12*N
        pdf(ofn, height=height, width=width)
        par(mai=c(2, 0.5, 0.5, 2))
        par(cex=0.7, las=2)
        plot(ddc, horiz=T, leaflab=leaflab)
        graphics.off()
    }
    pdend.nodes(style="clean")
    pdend.nodes(style="detailed")

}

rl=function()
{
    source("utils.r")
    source("cplot_utils.r")
    source("cluster_cycs.r")
}

ex=function()
{
    cluster.cycs()
}

pl=function()
{
    plot.clustering()
}


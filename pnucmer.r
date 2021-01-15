options(stringsAsFactors=F)

##############################################################################################################
# util functions
##############################################################################################################

smatrix2matrix=function(smatrix, dim, i.field="i", j.field="j", value.field="value", default.value=0)
{
  indices1 = smatrix[,i.field] + (smatrix[,j.field]-1) * dim
  indices2 = smatrix[,j.field] + (smatrix[,i.field]-1) * dim
  v = rep(default.value, dim*dim)
  v[indices1] = smatrix[,value.field]
  v[indices2] = smatrix[,value.field]
  matrix(v, dim, dim)
}

matrix2smatrix=function(matrix)
{
  dim1 = dim(matrix)[1]
  dim2 = dim(matrix)[2]
  v = as.vector(matrix)
  indices = 1:(dim1*dim2)
  i = (indices-1) %% dim1 + 1
  j = floor((indices-1) / dim1) + 1
  data.frame(i=i, j=j, value=v, stringsAsFactors=F)
}

make.color.panel=function(colors, ncols=256)
{
    library(gplots)
    panel = NULL
    for (i in 2:length(colors))
        panel = c(panel, colorpanel(ncols, colors[i-1], colors[i]))
    panel
}

vals.to.cols=function(vals, breaks, ncols=256)
{
  min = breaks[1]
  max = breaks[length(breaks)]
  vals = ifelse(vals < min, min, ifelse(vals>max, max, vals))
  n = length(breaks)-1
  cols = rep(-1, length(vals))
  for (i in 1:n)
  {
    ind = (breaks[i] <= vals) & (vals <= breaks[i+1])
    if (!any(ind))
      next
    # normalize to [0,1]
    cols[ind] = (vals[ind] - breaks[i]) / (breaks[i+1] - breaks[i])
    # normalize to [i*ncols,i*(ncols+1)]
    cols[ind] = (i-1)*ncols + cols[ind]*(ncols-1) + 1
    # round
    cols[ind] = floor(cols[ind])
  }
  return (cols)
}

##############################################################################################################
# figure functions
##############################################################################################################

get.min.dist=function(df, dfm)
{
    N = dim(df)[1]
    dfm = dfm[is.element(dfm$cyc1, df$cycle) &  is.element(dfm$cyc2, df$cycle),]
    dfm$i1 = match(dfm$cyc1, df$cycle)
    dfm$i2 = match(dfm$cyc2, df$cycle)
    dfm$dist = 1 - dfm$weighted_genome_dist
    mm = smatrix2matrix(smatrix=dfm, dim=N, i.field="i1", j.field="i2", value.field="dist", default.value=1)
    apply(mm, 1, min)
}

plot.scatter=function()
{
    idir = "/relman01/home/nshalon/work/pipe/sour/nucmer_summary"
    dfm = read.delim(paste0(idir, "/all_all_summary"))
    dfm = dfm[dfm$len1 >= 1000 & dfm$len2 >= 1000 & dfm$genome_dist > 0,]

    xlim = c(0.5, 1)
    xlim = range(dfm$genome_dist)
    ylim = range(dfm$ani)
    odir = "figures/cluster_nucmer/"
    system(paste("mkdir -p", odir))
    ofn = paste0(odir, "/scatter.pdf")
    cat(sprintf("plotting: %s\n", ofn))
    pdf(ofn, height=6, width=6)
    plot.new()
    plot.window(xlim=xlim, ylim=ylim)
    grid()

    axis(1)
    axis(2)
    title(xlab="aligned fraction", ylab="aligned identity")
    points(dfm$genome_dist, dfm$ani, pch=19, cex=0.5)

    graphics.off()


    xlim = c(0.94, 1)
    ylim = range(0.94, 1)
    ofn = paste0(odir, "/scatter_zoom.pdf")
    cat(sprintf("plotting: %s\n", ofn))
    pdf(ofn, height=6, width=6)
    plot.new()
    plot.window(xlim=xlim, ylim=ylim)
    grid()

    axis(1)
    axis(2)
    title(xlab="aligned fraction", ylab="aligned identity")
    points(dfm$genome_dist, dfm$ani, pch=19, cex=0.5)

    graphics.off()

    ofn = paste0(odir, "/hist_ani.pdf")
    cat(sprintf("plotting: %s\n", ofn))
    pdf(ofn, height=4, width=6)
    hist(dfm$ani,100)
    graphics.off()

    ofn = paste0(odir, "/hist_dist.pdf")
    cat(sprintf("plotting: %s\n", ofn))
    pdf(ofn, height=4, width=6)
    hist(dfm$genome_dist,100)
    graphics.off()
}


plot.nucmer=function()
{
    idir = "/relman01/home/nshalon/work/pipe/sour/nucmer_summary"
    dfm = read.delim(paste0(idir, "/all_all_summary"))

    df = data.frame(cycle=sort(unique(c(dfm$cyc1, dfm$cyc2))))
    df$length = dfm$len1[match(df$cycle, dfm$cyc1)]
    df$min.dist = get.min.dist(df=df, dfm=dfm)

    # filter for >1k
    df = df[df$length >= 1000,]

    # filter out isolated cyckes
    df = df[df$min.dist <= 0.5,]

    dfm$dist = 1 - dfm$weighted_genome_dist
    N = dim(df)[1]


    dfm = dfm[is.element(dfm$cyc1, df$cycle) &  is.element(dfm$cyc2, df$cycle),]
    dfm$i1 = match(dfm$cyc1, df$cycle)
    dfm$i2 = match(dfm$cyc2, df$cycle)

    mm = smatrix2matrix(smatrix=dfm, dim=N, i.field="i1", j.field="i2", value.field="dist", default.value=1)
    diag(mm) = 0
    hc = hclust(as.dist(mm), method="single")
    ct = cutree(hc, h=0.05)
    df$group = ct
    hc$labels = paste(df$group, df$cycle)

    df = df[hc$order,]
    df$index = 1:N
    df$label = paste(df$cycle, df$group)

    mmo = mm[hc$order,hc$order]
    sm = matrix2smatrix(mmo)
    lim = c(0,N)

    breaks = c(0, 0.95, 0.999, 1)
    cols = c("white", "darkblue", "red", "orange")
    panel = make.color.panel(colors=cols)
    sm$col = panel[vals.to.cols(vals=1-sm$value, breaks=breaks)]
    df$font = 1

    odir = "figures/cluster_nucmer/"
    system(paste("mkdir -p", odir))
    ofn = paste0(odir, "/matrix_all.pdf")
    cat(sprintf("plotting: %s\n", ofn))
    pdf(ofn, height=0.2*N, width=0.2*N)
    par(mai=c(2, 2, 0.5, 0.5))
    plot.new()
    plot.window(xlim=lim, ylim=lim)
    rect(xleft=sm$i-1, xright=sm$i, ybottom=sm$j-1, ytop=sm$j, col=sm$col, border=NA, lwd=0.5)

    mtext(side=1, text=df$label, at=df$index-0.5, las=2, cex=0.75, line=-0.5)
    mtext(side=2, text=df$label, at=df$index-0.5, las=2, cex=0.75, line=-0.5)
    graphics.off()

    ofn = paste0(odir, "/dendro_all.pdf")
    dd = as.dendrogram(hc)
    cat(sprintf("plotting: %s\n", ofn))
    pdf(ofn, height=0.2*N, width=5)
    par(mai=c(1, 0.2, 0.5, 3.5))
    par(cex=0.7, las=2)
    plot(dd, horiz=T)
    graphics.off()

    vv = dfm$weighted_genome_dist
    vv = sort(vv[vv>0])
    ofn = paste0(odir, "/pairwise_genomic.pdf")
    cat(sprintf("plotting: %s\n", ofn))
    pdf(ofn, height=4, width=8)
    vv = vv[vv>0.5]
    # plot(ecdf(vv), col.points=ifelse(vv>0.99, "red", ifelse(vv>0.95, "blue", "black")), cex.axis=0.7, las=2)
    plot(ecdf(vv), cex.axis=0.7, las=2, main="genomic")
    graphics.off()

    vv = dfm$ani
    vv = sort(vv[vv>0])
    ofn = paste0(odir, "/pairwise_ani.pdf")
    cat(sprintf("plotting: %s\n", ofn))
    pdf(ofn, height=4, width=8)
    vv = vv[vv>0.5]
    # cols = rainbow(length(vv))
    # cols = ifelse(vv>0.99, "red", ifelse(vv>0.95, "blue", "black"))
    # plot(ecdf(vv), do.points=T, col.points=cols, cex.axis=0.7, las=2)
    plot(ecdf(vv), cex.axis=0.7, las=2, main="ANI")
    graphics.off()

    vv = dfm$genome_dist
    vv = sort(vv[vv>0])
    ofn = paste0(odir, "/pairwise_breadth.pdf")
    cat(sprintf("plotting: %s\n", ofn))
    pdf(ofn, height=4, width=8)
    vv = vv[vv>0.5]
    # cols = rainbow(length(vv))
    # cols = ifelse(vv>0.99, "red", ifelse(vv>0.95, "blue", "black"))
    # plot(ecdf(vv), do.points=T, col.points=cols, cex.axis=0.7, las=2)
    plot(ecdf(vv), cex.axis=0.7, las=2, main="breadth")
    graphics.off()
}

rl=function()
{
    source("pnucmer.r")
}

ex=function()
{
    plot.nucmer()
}


##################################################################################################################
# cluster align profile
##################################################################################################################

segment.alignment=function(dfc, dfa, dfo, cols)
{
    dfa = dfa[dfa$key1 != dfa$key2,]
    dfc$key = paste(dfc$sample, dfc$cycle, sep="_")

    # add cycle length
    dfa$key1 = paste(dfa$sample1, dfa$cycle1, sep="_")
    dfa$key2 = paste(dfa$sample2, dfa$cycle2, sep="_")
    add.length=function(key) {
        ix = match(key, dfc$key)
        if (any(is.na(ix)))
            stop("not all cycles found")
        dfc$length[ix]
    }
    dfa$length1 = add.length(dfa$key1)
    dfa$length2 = add.length(dfa$key2)

    dfa$alen1 = abs(dfa$end1 - dfa$start1) + 1
    dfa$alen2 = abs(dfa$end2 - dfa$start2) + 1

    dfa$step1 = round(dfa$length1/length(cols))
    dfa$step2 = round(dfa$length2/length(cols))

    # breakdown into segments
    rr = NULL
    for (i in 1:dim(dfa)[1]) {

        if (dfa$alen2[i] <= dfa$step2[i]) {
            rr = rbind(rr, data.frame(dfa[i,c("start1", "end1", "start2", "end2", "sample1", "cycle1", "sample2", "cycle2", "length1", "length2")]))
            next
        }
        N = ceiling(abs(dfa$end2[i]-dfa$start2[i])/dfa$step2[i])+1

        gcoords=function(start, end) {
            strand = ifelse(end>start,+1,-1)
            coords = seq(start, end, length.out=N)
            rx = range(coords)
            dx = data.frame(start=coords[-N], end=coords[-1])
            force.range=function(vv) { ifelse(vv>rx[2],rx[2],ifelse(vv<rx[1],rx[1],vv)) }
            dx$start = force.range(dx$start)
            dx$end = force.range(dx$end)
            dx
        }
        ll1 = gcoords(dfa$start1[i], dfa$end1[i])
        ll2 = gcoords(dfa$start2[i], dfa$end2[i])

        dd = data.frame(start1=ll1$start, end1=ll1$end, start2=ll2$start, end2=ll2$end,
                        sample1=rep(dfa$sample1[i], N-1), cycle1=rep(dfa$cycle1[i], N-1),
                        sample2=rep(dfa$sample2[i], N-1), cycle2=rep(dfa$cycle2[i], N-1),
                        length1=rep(dfa$length1[i], N-1), length2=(rep(dfa$length2[i], N-1)))
        dd = dd[abs(dd$end1 - dd$start1)>0,]
        rr = rbind(rr, dd)
    }
    ix = match(paste(rr$sample1, rr$cycle1, rr$sample2, rr$cycle2, sep="_"), paste(dfo$key1, dfo$key2, sep="_"))
    if (any(is.na(ix)))
        stop("not all keys found")

    rr$frac1 = dfo$frac1[ix]
    rr$orig2 = dfo$coord2[ix]
    rr$strand1 = dfo$strand1[ix]
    rr$strand2 = dfo$strand2[ix]
    rr$center2 = (rr$end2 + rr$start2) / 2

    rr$frac2 = (rr$center2-1) / (rr$length2-1)
    rr$col.index = 1+ceiling((length(cols)-1)*rr$frac2)

    # rotate colors
    # rr$frac2 = ifelse(rr$strand2 == 1, (rr$center2-rr$orig2) / (rr$length2-1), (rr$orig2-rr$center2) / (rr$length2-1)) %% 1
    # rr$frac2.shift = ifelse(rr$strand1 == 1, rr$frac2 + rr$frac1, rr$frac2 - rr$frac1) %% 1
    # rr$col.index = 1+ceiling((length(cols)-1)*rr$frac2.shift)


    rr$col = cols[rr$col.index]
    rr
}

get.align.offset=function(df, df.cycs)
{
    df = df[df$key1 != df$key2,]
    ss = split(df, paste(df$key1, df$key2))
    rr = NULL
    for (i in 1:length(ss)) {
        dfx = ss[[i]]
        dfx$len2 = dfx$end2 - dfx$start2

        get.map=function(ii) {
            strand1 = ifelse(dfx$end1[ii] > dfx$start1[ii], 1, -1)
            strand2 = ifelse(dfx$end2[ii] > dfx$start2[ii], 1, -1)
            list(coord1=dfx$start1[ii], strand1=strand1,
                 coord2=dfx$start2[ii], strand2=strand2)
        }
        llf = get.map(which.min((dfx$start2)))
        llm = get.map(which.max(abs(dfx$len2)))

        rr = rbind(rr, data.frame(sample1=dfx$sample1[1], cycle1=dfx$cycle1[1],
                                  sample2=dfx$sample2[1], cycle2=dfx$cycle2[1],
                                  key1=dfx$key1[1], key2=dfx$key2[1],
                                  coord1=llm$coord1, strand1=llm$strand1,
                                  coord2=llm$coord2, strand2=llm$strand2,
                                  orig.coord=llf$coord1, orig.strand=llf$strand1))
    }
    rr$length1 = df.cycs$length[match(rr$key1, paste(df.cycs$sample, df.cycs$cycle, sep="_"))]
    rr$frac1 = ifelse(rr$strand1 == 1, (rr$coord1-1)/(rr$length1-1),(rr$length1-rr$coord1)/(rr$length1-1))
    rr
}

align.profile=function(base.height=0.05, title="coverage", dfc, table.cycles, table.align, table.snps, gap=0.15)
{
    n.members = dim(dfc)[1]-1
    height = base.height * n.members

    df.offset = get.align.offset(df=table.align, df.cycs=table.cycles)

    # align.cols = colorRampPalette(c("darkblue", colors()[43]))(100)
    align.cols = colorRampPalette(c(colors()[610], colors()[614]))(100)
    table.align.segs = segment.alignment(dfc=table.cycles, dfa=table.align, dfo=df.offset, cols=align.cols)

    nts = c("A", "C", "G", "T")
    cols = c("red", "green", "blue", "orange")
    ix = match(table.snps$nt2, nts)
    table.snps$col = ifelse(!is.na(ix), cols[ix], "gray")

    list(height=height, title=title, n.members=n.members, dfc=dfc,
         table.align=table.align, table.align.segs=table.align.segs, table.snps=table.snps, df.offset=df.offset,
         plot.f=function(cx, pr) {
             hh = (1-gap)*(pr$height / n.members)
             dfc = pr$dfc[!(pr$dfc$sample == cx$sample & pr$dfc$cycle == cx$cycle),]
             dfc$index = 1:dim(dfc)[1]
             dfc$ycenter = ((dfc$index-0.5) / pr$n.members) * pr$height
             restrict=function(table) {
                 if (is.null(table) || dim(table)[1] == 0) return (NULL)
                 ix = table$sample1 == cx$sample & table$cycle1 == cx$cycle & !(table$sample2 == cx$sample & table$cycle2 == cx$cycle)
                 if (sum(ix) == 0) return (NULL)
                 table = table[ix,]
                 table$key = paste0(table$sample2, "_", table$cycle2)
                 table$index = match(table$key, dfc$key)
                 table$ycenter = ((table$index-0.5) / pr$n.members) * pr$height
                 table$ybottom = table$ycenter - hh/2
                 table$ytop = table$ycenter + hh/2
                 table
             }

             # alignment
             table.align = restrict(pr$table.align)
             table.align.segs = restrict(pr$table.align.segs)
             df.offset = restrict(pr$df.offset)

             # snps
             table.snps = restrict(pr$table.snps)

             # color alignment
             for (i in 1:dim(table.align.segs)[1])
                 cplot.rect(cx=cx, xleft=table.align.segs$start1[i], xright=table.align.segs$end1[i],
                            ybottom=table.align.segs$ybottom[i], ytop=table.align.segs$ytop[i],
                            col=table.align.segs$col[i], border=table.align.segs$col[i], lwd=0.5)

             # rectangles around alignments
             for (i in 1:dim(table.align)[1]) {
                 ## cplot.rect(cx=cx, xleft=table.align$start1[i], xright=table.align$end1[i],
                 ##            ybottom=table.align$ybottom[i], ytop=table.align$ytop[i],
                 ##            col=NA, border=1, lwd=0.5)
                 cplot.vsegs(cx=cx, x=table.align$start1[i], y0=table.align$ybottom[i], y1=table.align$ytop[i],
                            col=1, lwd=0.5)
                 cplot.vsegs(cx=cx, x=table.align$end1[i], y0=table.align$ybottom[i], y1=table.align$ytop[i],
                            col=1, lwd=0.5)
             }

             # mark alignment origin
             ## for (i in 1:dim(df.offset)[1]) {
             ##     mark.width = if (cx$max.coord > 50000) cx$max.coord/200 else cx$max.coord/100
             ##     if (df.offset$orig.strand[i] == 1) {
             ##         mark.x = c(df.offset$orig.coord[i], df.offset$orig.coord[i], df.offset$orig.coord[i]+mark.width)
             ##         mark.y = c(df.offset$ybottom[i], df.offset$ytop[i], df.offset$ycenter[i])
             ##     } else {
             ##         mark.x = c(df.offset$orig.coord[i], df.offset$orig.coord[i], df.offset$orig.coord[i]-mark.width)
             ##         mark.y = c(df.offset$ybottom[i], df.offset$ytop[i], df.offset$ycenter[i])
             ##     }
             ##     cplot.polygon(cx=cx, x=mark.x, y=mark.y, col=1, border=NA)
             ## }

             if (!is.null(table.snps))
                 for (i in 1:dim(table.snps)[1])
                     cplot.vsegs(cx=cx, x=table.snps$base1[i], y0=table.snps$ybottom[i], y1=table.snps$ytop[i], col=table.snps$col[i], lwd=0.5)


             for (i in 1:n.members)
                 cplot.title(cx=cx, pr=pr, title=dfc$key[i], height=dfc$ycenter[i], cex=0.3)
    })
}

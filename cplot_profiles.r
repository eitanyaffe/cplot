##################################################################################################################
# profile implementations
##################################################################################################################

# each profile is a list that must contain:
# height: profile height
# title: string title
# plot.f(): plots profile using context

##################################################################################################################
# empty profile
##################################################################################################################

# plots horizontal line in the middle of profile
empty.profile=function(height=1, title="", lty=1)
{
    list(title=title, height=height, lty=lty, plot.f=function(cx, pr) { })
}

##################################################################################################################
# v and h line profile
##################################################################################################################

vlines.profile=function(title="vline.profile", lty=2, table, field)
{
    list(height=0, title=title, table=table, field=field, lty=lty, plot.f=function(cx, pr) {
        coords = pr$table[pr$table$sample==cx$sample & pr$table$cycle==cx$cycle,pr$field]
        if (length(coords) > 0)
            cplot.vlines(cx=cx, coords=coords, lty=pr$lty)
    })
}

##################################################################################################################
# gene profile
##################################################################################################################

# gene table must contain:
# start, end: coords
# strand: +1 or -1, to plot above or below
# col.field: color field of gene (if not supplied we color by strand)
# label: text label

gene.profile=function(height=1, title="line.profile", table, col.field, add.label, is.top, plot.trig=F)
{
    list(title=title, height=height, table=table, is.top=is.top, plot.trig=plot.trig,
         plot.f=function(cx, pr) {
             pr$table = pr$table[pr$table$sample==cx$sample & pr$table$cycle == cx$cycle,]


             trig.width = if (cx$max.coord > 50000) cx$max.coord/720 else cx$max.coord/360
             label.coords = (pr$table$start  + pr$table$end)/2

             cplot.rect(cx=cx, xleft=0, xright=cx$max.coord, ybottom=0, ytop=pr$height, col="lightgray", border=NA)
             if (cx$type == "circle" && !cx$multi) cplot.title(cx=cx, pr=pr, title=pr$title, cex=0.4)
             if (dim(pr$table)[1] == 0)
                 return (NULL)

             if (is.element(col.field, names(pr$table)))
                 pr$table$col = pr$table[,col.field]
             else
                 pr$table$col = ifelse(pr$table$strand == "+", "red", "blue")
             for (i in 1:dim(pr$table)[1]) {
                 cplot.rect(cx=cx, xleft=pr$table$start[i], xright=pr$table$end[i], ybottom=0, ytop=pr$height,
                            col=pr$table$col[i], border=NA)
                 if (pr$table$strand[i] == "+") {
                     trig.x = c(pr$table$start[i], pr$table$start[i], pr$table$start[i]+trig.width)
                     trig.y = c(0, height, height/2)
                 } else {
                     trig.x = c(pr$table$end[i], pr$table$end[i], pr$table$end[i]-trig.width)
                     trig.y = c(0, height, height/2)
                 }
                 if (plot.trig)
                     cplot.polygon(cx=cx, x=trig.x, y=trig.y, col=1, border=NA)
             }
             if (add.label && !cx$multi) {
                 labels = pr$table$label
                 ix = grepl("hypothetical", labels, ignore.case=T) |
                     grepl("uncharacterized", labels, ignore.case=T) |
                     grepl("no_hit", labels, ignore.case=T)
                 labels[ix] = ""
                 cplot.margin.text(cx=cx, pr=pr, x=label.coords, labels=labels, is.top=pr$is.top)
             }
         })
}

gene.identity.profile=function(height, table, add.label=F, is.top=T, plot.trig=F, odir.legend=NA)
{
    breaks = c(0, 20, 70, 98, 100)
    cols = c("gray", colorRampPalette(c(colors()[430], colors()[461]))(3))
    # cols = colorpanel(4, "white", "darkblue")
    # cols = c("gray", "darkgreen", "darkblue", "darkred", "orange")
    table$identity.col = vals2colors(vals=table$identity, breaks=breaks, cols=cols)
    if (!is.na(odir.legend))
        plot.legend.breaks(odir=odir.legend, cols=cols, breaks=breaks, title="uniref_idenity")

    gene.profile(height=height, title="identity", table=table, col.field="identity.col", plot.trig=plot.trig,
                 add.label=add.label, is.top=is.top)
}

gene.class.profile=function(height, table, cclass, col.list, add.label=F, is.top=T, plot.trig=F)
{
    table$class.col = ifelse(grepl(cclass, table$class), col.list[[cclass]], "gray")
    gene.profile(height=height, title=cclass, table=table, col.field="class.col",
                 add.label=add.label, is.top=is.top)
}

gene.strand.profile=function(height, table, add.label=F, is.top=T, plot.trig=F)
{
    table$strand.col = ifelse(table$strand == "+", "red", "blue")
    gene.profile(height=height, title="strand", table=table, col.field="strand.col", plot.trig=plot.trig,
                 add.label=add.label, is.top=is.top)
}

##################################################################################################################
# points profile
##################################################################################################################

# table must contain:
# x,y: vector of coords

plot.cov.f=function(cx, pr)
{
    MDC = pr$cycle.table$MDC[pr$cycle.table$sample==cx$sample & pr$cycle.table$cycle == cx$cycle]
    table = pr$cov.table
    table = table[table$sample==cx$sample & table$cycle == cx$cycle,]
    fields = pr$fields

    val.range = c(0, 1.1*max(table[,fields]))
    for (i in 1:length(fields)) {
        field = fields[i]
                 yfield = paste0(field, ".y")
        table[,yfield] = val2y(vals=table[,field], val.range=val.range, height=pr$height)
    }
    grid.vals = get.grid.at(val.range=val.range, grid.nlines=pr$grid.nlines)
    grid.at = val2y(vals=grid.vals, val.range=val.range, height=pr$height)

    MDC.y = val2y(vals=MDC, val.range=val.range, height=pr$height)
    cplot.hline(cx=cx, height=MDC.y, lty=pr$MDC.lty, col=pr$MDC.col)

    max.coord = max(table$base)
    binsize = if(max.coord > 10000) 100 else 10
    breaks = unique(c(seq(from=0, to=max.coord, by=binsize), max.coord))
    ccut = cut(table$base, breaks=breaks)
    N = length(breaks)
    rr = data.frame(from=breaks[-N], to=breaks[-1])
    rr$x = (rr$from + rr$to) / 2
    for (i in 1:length(fields))
        rr[,fields[i]] = pr$height * sapply(split(table[,fields[i]], ccut), median) / val.range[2]

    cplot.hgrid(cx=cx, pr=pr, lty=1, col="gray",
                vals=grid.vals, at=grid.at, add.label=pr$grid.label, cex=pr$axis.cex)

    for (i in 1:length(fields))
        cplot.lines(cx=cx, table=rr, yfield=fields[i], col=pr$cols[i], lwd=pr$lwd)
    cplot.box(cx=cx, pr=pr, col=1)
}

cov.profile=function(height=1, title="coverage", cov.table, cycle.table,
                     grid.label=F, grid.nlines, axis.cex=0.4, lwd=1,
                     MDC.lty=1, MDC.col="red")
{
    fields = c("support",
               "out_cov", "in_cov",
               "out_paired_weird", "in_paired_weird",
               "out_singleton", "in_singleton")
    cols = c("black", "blue", "darkblue", "red", "darkred", "green", "darkgreen")

    list(height=height, title=title, cycle.table=cycle.table, cov.table=cov.table, lwd=lwd,
         fields=fields, cols=cols, axis.cex=axis.cex, MDC.lty=MDC.lty, MDC.col=MDC.col,
         grid.label=grid.label, grid.nlines=grid.nlines,
         plot.f=plot.cov.f)
}

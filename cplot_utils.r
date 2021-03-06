
##################################################################################################################
# file util functions
##################################################################################################################

safe.read.delim=function(ifn, ...) {
    tryCatch(read.delim(ifn, ...), error = function(c) { NULL } )
}


read.cache=function(fn, read.f=read.delim)
{
    if (!exists(".cache")) .cache <<- list()
    if (!is.element(fn, names(.cache)))
        .cache[[fn]] <<- read.f(fn)
    .cache[[fn]]
}

##################################################################################################################
# general util functions
##################################################################################################################

label.sample=function(x)
{
    df = read.delim("data/sample_table.txt")
    ix = match(x, df$sample)
    if (any(is.na(ix) & x != "ref"))
        stop("some samples not found in labeling table")
    ifelse(!is.na(ix), df$label[ix], x)
}

sample.rect.points=function(xleft, xright, ybottom, ytop, npoints)
{
    xx = seq(xleft, xright, length.out=npoints)
    data.frame(x=c(xx, rev(xx)), y=c(rep(ybottom, npoints), rep(ytop, npoints)))
}

coord2angle=function(cx, x)
{
    -2*pi*(x/(cx$max.coord*cx$gap.factor)) + pi/2
}

rect2circle=function(cx, x, y, base.radius)
{
    angle = coord2angle(cx=cx, x=x)
    deg = (angle*180)/pi
    radius = base.radius + y
    px = cx$center.x + radius * cos(angle)
    py = cx$center.y + radius * sin(angle)
    data.frame(x=px, y=py, deg=deg)
}

val2y=function(vals, val.range, height)
{
    height * (vals - val.range[1]) / diff(val.range)
}

get.grid.at=function(val.range, grid.nlines)
{
    axis.vals = sort(c(10^(0:5), 2*10^(0:5), 5*10^(0:5)))
    ii = findInterval(diff(val.range) / grid.nlines, axis.vals)
    grid.delta = axis.vals[ii]
    grid.vals = seq(from=val.range[1], to=val.range[2], by=grid.delta)
    grid.vals[grid.vals<val.range[2]]
}

vals2colors=function(vals, cols, breaks, all.inside=T, ...)
{
    ii = findInterval(vals, breaks, all.inside=all.inside, ...)
    cols[ii]
}

restrict.table=function(table, cx)
{
    if (is.element("cycle", names(table)))
        table = table[table$cycle == cx$cycle,]
    if (is.element("sample", names(sample)))
        table = table[table$sample == cx$sample,]
    table
}

##################################################################################################################
# util plotting functions
##################################################################################################################

plot.legend=function(odir, cols, names, title)
{
    ofn = paste0(odir, "/", title, ".pdf")
    pdf(ofn, width=2, height=1+0.2*length(cols))
    par(mai=c(0,0,0.6,0))
    plot.new()
    title(main=title)
    legend("top", fill=cols, legend=names, border=NA, cex=1, box.lwd=NA)
    graphics.off()
}

plot.legend.breaks=function(odir, cols, breaks, title, max.bin.offset=0)
{
    names = paste(breaks[-length(breaks)], breaks[-1]-max.bin.offset, sep="-")

    # single value
    names = ifelse(breaks[-length(breaks)] !=  breaks[-1]-max.bin.offset, names, breaks[-length(breaks)])
    names[length(breaks)-1] = paste0(">=", breaks[length(breaks)-1])

    plot.legend(odir=odir, cols=cols, names=names, title=title)
}

cplot.points=function(cx, table, yfield, pch, col)
{
    hh = cx$current.height
    switch (cx$type,
            circle={
                pp = rect2circle(cx=cx, x=table$x, y=table[,yfield], base.radius=cx$base.rad+cx$current.height)
                points(x=pp$x, y=pp$y, pch=pch, col=col)
            },
            rect={
                points(x=table$x, y=hh+table[,yfield], pch=pch, col=col)
            },
            stop("unknown type"))
}

cplot.lines=function(cx, table, yfield, lwd, col)
{
    hh = cx$current.height
    switch (cx$type,
            circle={
                pp = rect2circle(cx=cx, x=table$x, y=table[,yfield], base.radius=cx$base.rad+cx$current.height)
                lines(x=pp$x, y=pp$y, col=col, lwd=lwd)
            },
            rect={
                lines(x=table$x, y=hh+table[,yfield], col=col, lwd=lwd)
            },
            stop("unknown type"))
}

cplot.vlines=function(cx, coords, lty=1)
{
    switch (cx$type,
            circle={
                pstart = rect2circle(cx=cx, x=coords, y=0, base.radius=cx$base.rad)
                pend = rect2circle(cx=cx, x=coords, y=cx$heights, base.radius=cx$base.rad)
                segments(x0=pstart$x, x1=pend$x, y0=pstart$y, y1=pend$y, lty=lty)
            },
            rect={
                segments(x0=coords, x1=coords, y0=0, y1=cx$heights, lty=lty)
            },
            stop("unknown type"))
}

cplot.vsegs=function(cx, x, y0, y1, ...)
{
    hh = cx$current.height
    switch (cx$type,
            circle={
                pstart = rect2circle(cx=cx, x=x, y=hh+y0, base.radius=cx$base.rad)
                pend = rect2circle(cx=cx, x=x, y=hh+y1, base.radius=cx$base.rad)
                segments(x0=pstart$x, x1=pend$x, y0=pstart$y, y1=pend$y, ...)
            },
            rect={
                segments(x0=x, x1=x, y0=hh+y0, y1=hh+y1, ...)
            },
            stop("unknown type"))
}


# horizontal line
cplot.hline=function(cx, height=0, lty=1, col=1)
{
    hh = cx$current.height + height
    switch (cx$type,
            circle=draw.arc(x=cx$center.x, y=cx$center.y, radius=cx$base.rad + hh,, n=cx$arc.n,
                            angle1=coord2angle(cx=cx, x=0), angle2=coord2angle(cx=cx, x=cx$max.coord),
                            lty=lty, col=col),
            rect=segments(x0=0, x1=cx$max.coord, y0=hh, y1=hh, lty=lty, col=col),
            stop("unknown type"))
}

cplot.hgrid=function(cx, pr, vals, at, add.label, lty=1, col=1, cex=0.5, tick.size=0.2)
{
    at = at + cx$current.height
    for (i in 1:length(at)) {
        at.i = at[i]
        switch (cx$type,
                circle=draw.arc(x=cx$center.x, y=cx$center.y, radius=cx$base.rad+at.i, n=cx$arc.n,
                            angle1=coord2angle(cx=cx, x=0), angle2=coord2angle(cx=cx, x=cx$max.coord),
                            lty=lty, col=col),
                rect=segments(x0=0, x1=cx$max.coord, y0=at.i, y1=at.i, lty=lty, col=col),
                stop("unknown type"))
    }
    if (add.label) {
        switch (cx$type,
                circle={
                    pp = rect2circle(cx=cx, x=0, y=at, base.radius=cx$base.rad)
                    if (!cx$multi) {
                        text(x=pp$x, y=pp$y, pos=2, labels=vals, offset=.2, cex=cex, adj=1)
                        segments(x0=pp$x, x1=pp$x-tick.size, y0=pp$y, y1=pp$y)
                    }
                },
                rect={
                    if (!cx$multi) mtext(text=vals, side=4, at=at, cex=cex, adj=0, las=1)
                    segments(x0=cx$max.coord, x1=cx$max.coord*1.01, y0=at, y1=at)

                },
                stop("unknown type"))
    }
}

cplot.title=function(cx, pr, title, cex=0.4, height=NULL, ...)
{
    hh = if (is.null(height)) pr$height/2 else height
    switch (cx$type,
            circle={
                pp = rect2circle(cx=cx, x=0, y=hh, base.radius=cx$base.rad+cx$current.height)
                text(x=pp$x, y=pp$y, pos=2, labels=title, offset=.2, cex=cex, adj=1, ...)
            },
            rect={
                mtext(text=title, side=4, at=cx$current.height+hh, cex=cex*1.5, adj=0, las=1, ...)

            },
            stop("unknown type"))
}

# plot polygon
cplot.polygon=function(cx, x, y, ...)
{
    hh = cx$current.height
    switch (cx$type,
            circle={
                pp = rect2circle(cx=cx, x=x, y=y, base.radius=(cx$base.rad+cx$current.height))
                polygon(x=pp$x, y=pp$y, ...)
            },
            rect=polygon(x=x, y=y+hh, ...),
            stop("unknown type"))
}

cplot.margin.text=function(cx, pr, x, is.top, labels, cex=0.25)
{
    hh = cx$current.height
    switch (cx$type,
            circle={
                yval = if (is.top) 2*pr$height else -1.5*pr$height
                pp = rect2circle(cx=cx, x=x, y=yval, base.radius=(cx$base.rad+cx$current.height))
                pp$adj = ifelse(pp$deg < -90, 1, 0)
                pp$deg = ifelse(pp$deg < -90, pp$deg + 180, pp$deg)
                for (i in 1:dim(pp)[1])
                    text(x=pp$x[i], y=pp$y[i], labels=labels[i], srt=pp$deg[i], xpd=NA, cex=cex, adj=pp$adj[i])
            },
            rect={
                pos = if (is.top) 3 else 1
                adj = if (is.top) -0.5 else +0.5
                text(x=x, y=cx$current.height+pr$height, labels=labels, adj=adj,
                     srt=45, xpd=NA, cex=cex)
            },
            stop("unknown type"))
}

# box around profile
cplot.box=function(cx, pr, col=1)
{
    cplot.rect(cx=cx, xleft=0, xright=cx$max.coord, ybottom=0, ytop=pr$height, border=col)
}

# plot rectangle
cplot.rect=function(cx, xleft, xright, ybottom, ytop, ...)
{
    hh = cx$current.height
    switch (cx$type,
            circle={
                degrees = min(360, (360 * abs(xright - xleft)) / cx$max.coord)
                npoints = 2 * max(1,round(degrees))
                pp = sample.rect.points(xleft=xleft, xright=xright, ybottom=ybottom, ytop=ytop, npoints=npoints)
                cplot.polygon(cx=cx, x=pp$x, y=pp$y, ...)
            },
            rect={
                rect(xleft=xleft, xright=xright, ybottom=ybottom+hh, ytop=ytop+hh, ...)
            },
            stop("unknown type"))
}

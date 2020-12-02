library("plotrix")
source("cplot_utils.r")
source("cplot_profiles.r")
options(stringsAsFactors=F)

##################################################################################################################
# core
##################################################################################################################

get.height.sum=function(profiles)
{
    rr = 0
    for (pr in profiles) {
        rr = rr + pr$height
    }
    rr
}

get.height.vec=function(profiles)
{
    rr = NULL
    for (pr in profiles)
        rr = c(rr, pr$height)
    rr
}

cplot=function(cx, profiles)
{
    cx$current.height = 0
    for (i in 1:length(profiles)) {
        pr = profiles[[i]]
        # cat(sprintf("plotting profile: %s\n", pr$title))

        if (cx$type == "rect" && pr$height>0)
            text(x=0, y=cx$current.height+pr$height/2, labels=pr$title, pos=2, xpd=NA)

        pr$plot.f(cx=cx, pr=pr)
        cx$current.height = cx$current.height + pr$height
    }

}

##################################################################################################################
# main
##################################################################################################################

cplot.on.circle=function(sample, cycle, ofn, label, base.rad, max.coord, profiles, width, height, extra=1.5)
{
    total.rad = base.rad + get.height.sum(profiles=profiles)

    xlim = c(-total.rad*extra,total.rad*extra)
    ylim = c(-total.rad*extra,total.rad*extra)

    pdf(ofn, width=width, height=height)
    plot.new()
    par(mai=c(0,0,0,0))
    plot.window(xlim=xlim, ylim=ylim)
    text(x=0, y=0, labels=paste(label, collapse="\n"), cex=0.75)

    cx = list(sample=sample, cycle=cycle, max.coord=max.coord, center.x=0, center.y=0,
              heights=get.height.sum(profiles=profiles), base.rad=base.rad,
              gap.factor=1.05, arc.n=0.01,
              type="circle", border.col="gray", npoints=360*2, multi=F)
    cplot(cx=cx, profiles=profiles)

    if (cx$type == "circle") {
        cplot.vlines(cx=cx, coords=0, lty=1)
        cplot.vlines(cx=cx, coords=cx$max.coord, lty=1)
    }

    graphics.off()
}

cplot.on.rect=function(sample, cycle, ofn, label, max.coord, profiles, width, height)
{
    total.heights = get.height.sum(profiles=profiles)

    xlim = c(0, max.coord)
    ylim = c(0, total.heights)

    pdf(ofn, width=width, height=height)
    plot.new()
    par(mai=c(0.2,0.8,1.6,0.2))
    plot.window(xlim=xlim, ylim=ylim)
    title(main=paste(label, collapse=", "))

    cx = list(sample=sample, cycle=cycle, max.coord=max.coord, type="rect", border.col="gray",
              heights=get.height.sum(profiles=profiles), npoints=360*2, multi=F)
    cplot(cx=cx, profiles=profiles)

    graphics.off()
}

cplot.singles=function(profiles, df, cc, base.rad, odir.circle, odir.rect)
{
    cat(sprintf("plotting %d cycles into %s and %s\n", dim(df)[1], odir.circle, odir.rect))
    for (i in dim(df)[1]:1) {
        sample = df$sample[i]
        id = df$id[i]
        cycle = df$cycle[i]
        len = df$length[i]
        cclass = df$classification[i]
        ofn.circle = paste0(odir.circle, "/", id, "_", cycle, ".pdf")
        ofn.rect = paste0(odir.rect, "/", id, "_", cycle, ".pdf")
        ncontigs = sum(cc$sample == sample & cc$cycle == cycle)+1

        label = c(sample, cycle, cclass,
                  paste0(len, "bp"),
                  paste0(ncontigs, " contigs"),
                  paste0("bottle=", df$bottleneck[i]),
                  paste0("score=", df$score[i]))

        cplot.on.circle(sample=sample, cycle=cycle, ofn=ofn.circle, label=label, max.coord=len,
                        base.rad=base.rad, profiles=profiles,
                        width=5, height=5)
        cplot.on.rect(sample=sample, cycle=cycle, ofn=ofn.rect, label=label, max.coord=len, profiles=profiles,
                      width=8, height=4)
    }

}

cplot.multi=function(profiles, df, cc, base.rad, plot.height, extra=1.1, ofn)
{
    cat(sprintf("plotting %d cycles in single plot: %s\n", dim(df)[1], ofn))
    # df$factor = df$length / min(df$length)
    # df$base.rad = base.rad * df$factor

    ring.rad = get.height.sum(profiles=profiles)
    total.rad = extra*(base.rad+ring.rad)

    NN = dim(df)[1]
    Nx = ceiling(sqrt(NN))
    Ny = ceiling(NN/Nx)
    df$index = 1:NN
    df$loc.x = (df$index-1) %% Nx + 1
    df$loc.y = (df$index-df$loc.x)/Ny+1
    df$center.x = 2 * (df$loc.x-1) * total.rad + total.rad
    df$center.y = 2 * (Ny - df$loc.y) * total.rad + total.rad

    xlim = c(0, 2*Nx*total.rad)
    ylim = c(0, 2*Ny*total.rad)
    # xlim = c(0, max(df$center.x+df$total.rad))
    # ylim = c(-max(df$total.rad), max(df$total.rad))

    #df$center.x = df$total.rad + c(0, cumsum(2*df$total.rad[-dim(df)[1]]))
    # df$center.y = 0
    # xlim = c(0, max(df$center.x+df$total.rad))
    # ylim = c(-max(df$total.rad), max(df$total.rad))

    plot.ratio = diff(xlim) / diff(ylim)
    pdf(ofn, width=plot.height*plot.ratio, height=plot.height)
    plot.new()
    par(mai=c(0,0,0,0))
    plot.window(xlim=xlim, ylim=ylim)

    for (i in 1:dim(df)[1]) {
        cycle = df$cycle[i]
        id = df$id[i]
        sample = df$sample[i]
        MDC = df$MDC[i]
        len = df$length[i]
        cclass = df$classification[i]
        ncontigs = sum(cc$cycle == cycle)+1

        label = c(id, paste0(round(MDC),"x", collapse=""), paste0(round(len/1000),"kb", collapse=""))

        text(x=df$center.x[i], y=df$center.y[i], labels=paste(label, collapse="\n"), cex=0.5)

        cx = list(sample=sample, cycle=cycle,
                  max.coord=len, center.x=df$center.x[i], center.y=df$center.y[i],
                  heights=get.height.sum(profiles=profiles), base.rad=base.rad,
                  gap.factor=1, arc.n=0.01,
                  type="circle", border.col="gray", npoints=360*2, multi=T)
        cplot(cx=cx, profiles=profiles)
    }
    graphics.off()
}

##################################################################################################################
# unit tests
##################################################################################################################

# basic test function
pt.test=function()
{
    genes.df = data.frame(start=c(5, 20, 50, 80), end=c(10, 30, 70, 85), strand=c("+","-","+","+"),
                          label=paste0("g", 1:4))
    line.df = data.frame(x=c(5, 20, 50)+5, y=c(2, 5, 9))

    profiles = list(
        p1=points.profile(height=10, title="vals", table=line.df),
        p2=empty.profile(height=1),
        p3=gene.profile(height=1, title="genes", table=genes.df),
        lines=vlines.profile(coords=c(genes.df$start,genes.df$end), lty=2))

    cplot.main(ofn.rect="test_rect.pdf", ofn.circle="test_circle.pdf", label="test", max.coord=90, profiles=profiles, width.circ=5, height.circ=5, width.rect=5, height.rect=2)
}

rl.test=function()
{
    source("cplot.r")
}

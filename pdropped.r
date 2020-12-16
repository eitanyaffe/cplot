options(stringsAsFactors=F)

dropped.label.f=function(df, i)
{
    c(df$cycle[i], paste0(df$length[i], "bp", sep=""),
      paste0("global=", round(df$score[i],1), ", P=", df$pval[i]),
      paste0("local=", round(df$singleton_score[i],1), ", P=", df$singleton_score_pval[i]))
}

plot.dropped=function(idir="/oak/stanford/groups/relman/users/nshalon/pipe/cyc_find_out4/cipro_clean_0.0")
{
    cat(sprintf("input dir: %s\n", idir))

    # load tables
    dff = read.delim("data/dropped.df")
    df = read.delim("/oak/stanford/groups/relman/users/nshalon/pipe/cyc_find_out4/cipro_clean_0.0/dominant_cycles_filter/cycle_stats_singletons")
    df = df[is.element(df$cycle, dff$cycle),]
    cc = read.delim("data/dropped_cc")
    covs = read.cache(paste0(idir, "/cycle_covs_long_adjusted"))

    sample = "AAB"
    cc$sample = sample
    covs$sample = sample
    df$sample = sample

    cc = cc[cc$cum_sum != 1,]
    df$id = sample

    # internal radius of circles
    base.rad = 3

    odir.base = "figures"
    odir.sample = paste0(odir.base, "/dropped/", sample)
    odir.circle = paste0(odir.sample, "/circles")
    odir.rect = paste0(odir.sample, "/rects")
    odir.legend = paste0(odir.base, "/legends")
    system(paste("mkdir -p", odir.circle, odir.rect, odir.legend))
    class.col.list = list(mobile="darkgreen", plasmid="darkred", phage="orange")
    plot.legend(odir=odir.legend, cols=unlist(class.col.list), names=names(class.col.list), title="gene_class")

    cov.col.list = list(support="black", out_cov="blue", in_cov="darkblue", out_paired_weird="red", in_paired_weird="darkred", out_singleton="green", in_singleton="darkgreen")
    plot.legend(odir=odir.legend, cols=unlist(cov.col.list), names=names(cov.col.list), title="coverage_class")
    profiles = list(
        cov.profile(height=3, title="cov", cov.table=covs,
                    cycle.table=df, grid.nlines=5, grid.label=T, plot.MDC=F, tick.size=0.1, lwd=2),
        vlines.profile(table=cc, field="cum_sum", lty=2))
    cplot.singles(profiles=profiles, df=df, cc=cc, base.rad=base.rad,
                  odir.circle=odir.circle, odir.rect=odir.rect, label.f=dropped.label.f)

    ofn = paste0(odir.base, "/dropped/all.pdf")
    cplot.multi(profiles=profiles, df=df, cc=cc, base.rad=base.rad, label.f=dropped.label.f,
                plot.height.per.cycle=4, style="r", extra=1.1, ofn=ofn, multi.flag=F, gap.factor=1.05)
}

rl.dropped=function()
{
    source("cplot.r")
    source("pdropped.r")
}

ex.dropped=function()
{
    plot.dropped()
}

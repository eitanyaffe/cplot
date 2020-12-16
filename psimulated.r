options(stringsAsFactors=F)

sim.label.f=function(df, i)
{
    c(paste0("coverage: ", round(df$median_support[i]),"x", collapse=""),
      paste0("bottleneck: ", df$bottleneck[i]),
      paste0("score: ", df$score[i]))
}

plot.simulated=function(idir="data/fig3_ex_plasmid", sample="sim.plasmid")
{
    cat(sprintf("input dir: %s\n", idir))

    # load tables
    df = read.cache(paste0(idir, "/dominant_cycle_stats_singletons"))
    cc = read.cache(paste0(idir, "/cycle_contig_table"))
    covs = read.cache(paste0(idir, "/cycle_covs_long"))

    cc$sample = sample
    covs$sample = sample
    df$sample = sample

    cc = cc[cc$cum_sum != 1,]
    df$id = sample

    # internal radius of circles
    base.rad = 3

    odir.base = "figures"
    odir.sample = paste0(odir.base, "/example/", sample)
    odir.circle = paste0(odir.sample, "/circles")
    odir.rect = paste0(odir.sample, "/rects")
    odir.legend = paste0(odir.base, "/legends")
    system(paste("mkdir -p", odir.circle, odir.rect, odir.legend))
    class.col.list = list(mobile="darkgreen", plasmid="darkred", phage="orange")
    plot.legend(odir=odir.legend, cols=unlist(class.col.list), names=names(class.col.list), title="gene_class")

    covs$non_support = covs$out_cov + covs$in_cov + covs$out_paired_weird + covs$in_paired_weird + covs$out_singleton + covs$in_singleton
    profiles = list(
        cov.profile(height=3, title="cov", cov.table=covs,
                    cols.list = list(support="black", non_support="red"),
                    fields = c("support", "non_support"),
                    cycle.table=df, grid.nlines=5, grid.label=T, plot.MDC=F, tick.size=0.1, lwd=2),
        vlines.profile(table=cc, field="cum_sum", lty=2))
    cplot.singles(profiles=profiles, df=df, cc=cc, base.rad=base.rad,
                  odir.circle=odir.circle, odir.rect=odir.rect, label.f=sim.label.f)
}

rl.sim=function()
{
    source("cplot.r")
    source("psimulated.r")
}

ex.sim=function()
{
    plot.simulated()
}

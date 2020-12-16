options(stringsAsFactors=F)

plot.sample=function(sample="cipro_clean", sample.title="AAB", length.filter="only_long")
{
    idir = paste0("/oak/stanford/groups/relman/users/nshalon/pipe/cyc_find_out4/", sample, "_0.0")
    cat(sprintf("loading sample %s from dir: %s\n", sample, idir))

    # load tables
    df = read.cache(paste0(idir, "/dominant_cycles_filter/cycle_summary_classification"))
    cc = read.cache(paste0(idir, "/cycle_contig_table"))

    if (file.exists(paste0(idir, "/dominant_cycles_filter/prodigal/v1/gene.tab"))) {
        gene.df = read.cache(paste0(idir, "/dominant_cycles_filter/prodigal/v1/gene.tab"))
        uniref = read.cache(paste0(idir, "/dominant_cycles_filter/uniref/v1/2020_07/table_uniq_taxa"))
    } else {
        gene.df = read.cache(paste0(idir, "/dominant_cycles/prodigal/v1/gene.tab"))
        uniref = read.cache(paste0(idir, "/dominant_cycles/uniref/v1/2020_07/table_uniq_taxa"))
    }

    gene.class = read.cache(paste0(idir, "/dominant_cycles/gene_classification"))

    # For Nitan: Please add to the pipeline a step that creates a version of this table that is
    # restricted to dominant cycles, to enhance performance
    covs = read.cache(paste0(idir, "/cycle_covs_long_adjusted"))

    cc$sample = sample
    gene.df$sample = sample
    gene.df$cycle = gene.df$contig
    covs$sample = sample

    # !!! overriding value in table. For example, in sample is SENB in /oak/stanford/groups/relman/users/nshalon/pipe/cyc_find_out4/sewage_b_0.0/dominant_cycles_filter/cycle_summary_classification
    df$sample = sample
    df = df[order(df$MDC, decreasing=T),]

    # uniref
    ix = match(gene.df$gene, uniref$gene)
    gene.df$desc = ifelse(!is.na(ix), uniref$prot_desc[ix], "no_hit")
    gene.df$taxa = ifelse(!is.na(ix), uniref$tax[ix], "no_hit")
    gene.df$identity = ifelse(!is.na(ix), uniref$identity[ix], 0)
    gene.df$uniref_count = ifelse(!is.na(ix), uniref$uniref_count[ix], 0)
    gene.df$label = gene.df$gene

    # class
    gene.df$class = ""
    for (i in 1:dim(gene.df)[1]) {
        gene.df$class[i] = paste(sort(unique(gene.class$classification[gene.class$gene == gene.df$gene[i]])),
                                 collapse=",")
    }

    cc = cc[cc$cum_sum != 1,]

    df = switch(length.filter,
                only_short=df[df$length < 1000,],
                only_long=df[df$length >= 1000,],
                all=df)
    df$id = sample.title

    # internal radius of circles
    base.rad = 5

    odir.base = "figures"
    odir.sample = paste0(odir.base, "/samples/", sample.title, "/", length.filter)
    odir.grouped = paste0(odir.base, "/samples/sets")
    odir.circle = paste0(odir.sample, "/circles")
    odir.rect = paste0(odir.sample, "/rects")
    odir.legend = paste0(odir.base, "/legends")
    system(paste("mkdir -p", odir.circle, odir.rect, odir.legend, odir.grouped))
    class.col.list = list(mobile="darkgreen", plasmid="darkred", phage="orange")
    plot.legend(odir=odir.legend, cols=unlist(class.col.list), names=names(class.col.list), title="gene_class")

    #######################################################################################
    # plot on cycle per plot
    #######################################################################################

    profiles.detailed = list(
        cov.profile(height=3, title="cov", cov.table=covs, cycle.table=df, grid.nlines=5, grid.label=T),
        empty.profile(height=0.5),
        gene.identity.profile(height=0.4, table=gene.df, odir.legend=odir.legend, plot.trig=T),
        empty.profile(height=0.1),
        gene.uniref.count.profile(height=0.4, table=gene.df, odir.legend=odir.legend),
        empty.profile(height=0.5),
        gene.class.profile(height=0.4, table=gene.df, cclass="mobile", col.list=class.col.list),
        empty.profile(height=0.1),
        gene.class.profile(height=0.4, table=gene.df, cclass="plasmid", col.list=class.col.list),
        empty.profile(height=0.1),
        gene.class.profile(height=0.4, table=gene.df, cclass="phage", col.list=class.col.list, add.label=T),
        vlines.profile(table=cc, field="cum_sum", lty=2))
    cplot.singles(profiles=profiles.detailed, df=df, cc=cc, base.rad=base.rad,
                  odir.circle=odir.circle, odir.rect=odir.rect)

    #######################################################################################
    # plot all cycles of sample together
    #######################################################################################

    profiles.multi = list(
        cov.profile(height=3, title="cov", cov.table=covs, cycle.table=df, grid.nlines=0),
        empty.profile(height=0.5),
        gene.identity.profile(height=0.4, table=gene.df, plot.trig=F),
        empty.profile(height=0.1),
        gene.uniref.count.profile(height=0.4, table=gene.df, odir.legend=odir.legend),
        empty.profile(height=0.5),
        gene.class.profile(height=0.4, table=gene.df, cclass="mobile", col.list=class.col.list),
        empty.profile(height=0.1),
        gene.class.profile(height=0.4, table=gene.df, cclass="plasmid", col.list=class.col.list),
        empty.profile(height=0.1),
        gene.class.profile(height=0.4, table=gene.df, cclass="phage", col.list=class.col.list, add.label=T),
        vlines.profile(table=cc, field="cum_sum", lty=2))
    ofn.sample = paste0(odir.grouped, "/", sample.title, "_", length.filter, ".pdf")
    cplot.multi(profiles=profiles.multi, df=df, cc=cc,
                base.rad=base.rad, plot.height.per.cycle=1.5, extra=1.1, style="r",
                ofn=ofn.sample)
}

plot.samples=function()
{
    length.filter = "all"
    ll = list(AAB.1="cipro_clean", AAB.2="cipro_b", FP.1="FP", FP.2="FP_B")
    for (i in 1:length(ll))
        plot.sample(sample.title=names(ll)[i], sample=ll[[i]], length.filter=length.filter)
}

rls=function()
{
    source("cplot.r")
    source("psamples.r")
}

exs=function()
{
    plot.samples()
}



plot.sample=function(sample.id, df, covs, cc, gene.df, uniref, gene.class, odir)
{
    # internal radius of circles
    base.rad = 5

    odir.circle = paste0(odir, "/", sample.id, "/circles")
    odir.rect = paste0(odir, "/", sample.id, "/rects")
    odir.grouped = paste0(odir, "/sets")
    odir.legend = paste0(odir, "/legends")
    system(paste("mkdir -p", odir.circle, odir.rect, odir.legend, odir.grouped))

    class.col.list = list(mobile="darkgreen", plasmid="darkred", phage="orange")
    plot.legend(odir=odir.legend, cols=unlist(class.col.list), names=names(class.col.list), title="gene_class")

    #######################################################################################
    # plot each cycle in plot
    #######################################################################################

    profiles.detailed = list(
        cov.profile(height=3, title="cov", cov.table=covs, cycle.table=df, grid.nlines=5, grid.label=T),
        empty.profile(height=0.1),
        gene.identity.profile(height=0.4, table=gene.df, odir.legend=odir.legend, plot.trig=T),
        empty.profile(height=0.1),
        gene.class.profile(height=0.4, table=gene.df, cclass="mobile", col.list=class.col.list),
        empty.profile(height=0.05),
        gene.class.profile(height=0.4, table=gene.df, cclass="plasmid", col.list=class.col.list),
        empty.profile(height=0.05),
        gene.class.profile(height=0.4, table=gene.df, cclass="phage", col.list=class.col.list, add.label=T),
        vlines.profile(table=cc, field="cum_sum", lty=2))
    cplot.singles(profiles=profiles.detailed, df=df, cc=cc, base.rad=base.rad,
                  odir.circle=odir.circle, odir.rect=odir.rect)

    #######################################################################################
    # plot all cycles of sample together
    #######################################################################################

    profiles.multi = list(
        cov.profile(height=3, title="cov", cov.table=covs, cycle.table=df, grid.nlines=0),
        empty.profile(height=0.1),
        gene.identity.profile(height=0.6, table=gene.df, plot.trig=F),
        empty.profile(height=0.1),
        gene.class.profile(height=0.6, table=gene.df, cclass="mobile", col.list=class.col.list),
        empty.profile(height=0.05),
        gene.class.profile(height=0.6, table=gene.df, cclass="plasmid", col.list=class.col.list),
        empty.profile(height=0.05),
        gene.class.profile(height=0.6, table=gene.df, cclass="phage", col.list=class.col.list, add.label=T),
        vlines.profile(table=cc, field="cum_sum", lty=2))
    ofn.sample = paste0(odir.grouped, "/", sample.id, ".pdf")
    cplot.multi(profiles=profiles.multi, df=df, cc=cc,
                base.rad=base.rad, plot.height=8, extra=1.1,
                ofn=ofn.sample)
}

plot.main=function(sample.id="aab")
{
    df = read.delim("data/aab/cycle_summary_classification")
    covs = read.delim("data/aab/cycle_covs_doms")
    cc = read.delim("data/aab/cycle_contig_table")
    gene.df = read.delim("data/aab/gene.tab")
    uniref = read.delim("data/aab/table_uniq_taxa")
    gene.class = read.delim("data/aab/gene_classification")

    df = df[order(df$length, decreasing=T),]
    gene.df$cycle = gene.df$contig

    # uniref
    ix = match(gene.df$gene, uniref$gene)
    gene.df$desc = ifelse(!is.na(ix), uniref$prot_desc[ix], "no_hit")
    gene.df$taxa = ifelse(!is.na(ix), uniref$tax[ix], "no_hit")
    gene.df$identity = ifelse(!is.na(ix), uniref$identity[ix], 0)
    #gene.df$label = gene.df$gene
    gene.df$label = paste(gene.df$gene, gene.df$desc)

    # class
    gene.df$class = ""
    for (i in 1:dim(gene.df)[1]) {
        gene.df$class[i] = paste(sort(unique(gene.class$classification[gene.class$gene == gene.df$gene[i]])),
                                      collapse=",")
    }

    cc = cc[cc$cum_sum != 1,]

    # should be sample-specific
    df$sample = sample.id
    gene.df$sample = sample.id
    covs$sample = sample.id

    # id should be shorter than the CYC names, used for plotting throughout
    df$id = paste0("c", 1:dim(df)[1])

    plot.sample(sample.id=sample.id, df=df, covs=covs, cc=cc, gene.df=gene.df, uniref=uniref, gene.class=gene.class, odir="figures")
}

rl=function()
{
    source("cplot.r")
    source("cplot_main.r")
}

ex=function()
{
    plot.main()
}

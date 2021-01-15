options(stringsAsFactors=F)
source("align_utils.r")
.ver = "v2"

cluster.label.f=function(df, i)
{
    c(df$cluster.id[i], df$label[i])
    ## c(df$cluster.id[i], df$prev.cluster[i],
    ##   paste0(df$label[i]),
    ##   df$classification[i],
    ##   paste0(round(df$MDC[i]),"x", collapse=""),
    ##   paste0(df$length[i], "bp"))
}

plot.cluster.internal=function(cluster.id, prev.cluster, dfc, ll, plot.style, odir)
{
    # internal radius of circles
    base.rad = 6

    odir.legend = paste0(odir, "/legends")
    system(paste("mkdir -p", odir.legend))

    class.col.list = list(mobile="darkgreen", plasmid="darkred", phage="orange")
    plot.legend(odir=odir.legend, cols=unlist(class.col.list), names=names(class.col.list), title="gene_class")

    ll$df$prev.cluster = prev.cluster
    ll$df$cluster.id = cluster.id

    is.clean = (grepl("clean", plot.style))
    pgap = 0.1

    profiles.short = list(
        cov.profile(height=1.2, cov.table=ll$covs, cycle.table=ll$df, grid.nlines=4, grid.label=T),
        empty.profile(height=pgap),
        align.profile(base.height=0.15, dfc=dfc, table.cycles=ll$df, table.align=ll$dfx, table.snps=ll$dfs,
                      odir.legend=odir.legend, plot.titles=!is.clean, plot.segments=F),
        empty.profile(height=pgap),
        gene.arrow.profile(height=0.1, table=ll$gene.df, add.label=!is.clean),
        vlines.profile(table=ll$cc, field="cum_sum", lty=2))

    profiles.details = list(
        cov.profile(height=1.2, cov.table=ll$covs, cycle.table=ll$df, grid.nlines=4, grid.label=T),
        empty.profile(height=pgap),
        align.profile(base.height=0.1, dfc=dfc, table.cycles=ll$df, table.align=ll$dfx, table.snps=ll$dfs,
                      odir.legend=odir.legend, plot.titles=T, plot.segments=T),
        empty.profile(height=pgap),
        gene.identity.profile(height=0.25, table=ll$gene.df, plot.trig=F, add.label=F, is.top=F),
        empty.profile(height=pgap),
        gene.uniref.count.profile(height=0.25, table=ll$gene.df, odir.legend=odir.legend),
        empty.profile(height=pgap),
        gene.class.profile(height=0.25, table=ll$gene.df, cclass="mobile", col.list=class.col.list),
        empty.profile(height=0.05),
        gene.class.profile(height=0.25, table=ll$gene.df, cclass="plasmid", col.list=class.col.list),
        empty.profile(height=0.05),
        gene.class.profile(height=0.25, table=ll$gene.df, cclass="phage", col.list=class.col.list),
        empty.profile(height=pgap),
        gene.arrow.profile(height=0.1, table=ll$gene.df, add.label=!is.clean),
        vlines.profile(table=ll$cc, field="cum_sum", lty=2))

    odir = paste0(odir, "/", plot.style)
    cat(sprintf("plotting %s to odir: %s\n", plot.style, odir))

    if (!grepl("rep_", plot.style)) {
        odir.base = paste0(odir, "/", cluster.id)
    } else {
        odir.base = paste0(odir)
    }
    odir.circle = paste0(odir.base, "/circles")
    odir.rect = paste0(odir.base, "/rects")
    system(paste("mkdir -p", odir.circle, odir.rect))

    rect.mai = c(0.2, 0.75, 0.75, 2.8)
    rect.height.factor = 0.75
    rect.width = 1.2 * (median(dfc$length)/2000) + rect.mai[2] + rect.mai[4]
    if (plot.style == "detailed") {
        cplot.singles(profiles=profiles.details, df=ll$df, cc=ll$cc, base.rad=base.rad, label.f=cluster.label.f,
                      odir.circle=odir.circle, odir.rect=odir.rect, circle.inch=7,
                      rect.width=rect.width, rect.height.factor=rect.height.factor, rect.mai=rect.mai)
    } else if (plot.style == "clean" || plot.style == "label") {
        cplot.singles(profiles=profiles.short, df=ll$df, cc=ll$cc, base.rad=base.rad, label.f=cluster.label.f,
                      odir.circle=odir.circle, odir.rect=odir.rect, circle.inch=7,
                      rect.width=rect.width, rect.height.factor=rect.height.factor, rect.mai=rect.mai)
    } else if (plot.style == "rep_label" || plot.style == "rep_clean") {
        if (any(ll$df$sample == "cipro_clean"))
            df = ll$df[ll$df$sample == "cipro_clean",]
        else
            df = ll$df[which.max(ll$df$length),]
        df$id = cluster.id
        cplot.singles(profiles=profiles.short, df=df, cc=ll$cc, base.rad=base.rad, label.f=cluster.label.f,
                      odir.circle=odir.circle, odir.rect=odir.rect, circle.inch=7,
                      rect.width=rect.width, rect.height.factor=rect.height.factor, rect.mai=rect.mai)
    } else {
        stop("unknown plot style")
    }
}

plot.cluster=function(dfc, prev.cluster, cluster.id)
{
    idir.cluster = "/relman01/home/nshalon/work/pipe/sour/cluster"
    idir.nucmer = "/relman01/home/nshalon/work/pipe/sour/nucmer"

    dfc = dfc[dfc$cluster.id == cluster.id,]
    dfc$key = paste0(dfc$sample, "_", dfc$cycle)

    dfx = read.cache(paste0(idir.nucmer, "/circ_refs/circ_refs.coords"), function(fn) { read.delim(fn, header=F) })
    dfx = parse.nucmer.coords(dfx)
    dfx = dfx[is.element(dfx$key1, dfc$key) & is.element(dfx$key2, dfc$key),]

    dfs = read.cache(paste0(idir.nucmer, "/circ_refs/circ_refs.snps"), function(fn) { safe.read.delim(fn, header=F) })
    dfs = parse.nucmer.snps(dfs)
    dfs = dfs[is.element(dfs$key1, dfc$key) & is.element(dfs$key2, dfc$key),]

    ll = list(dfx=dfx, dfs=dfs, df=NULL, cc=NULL, gene.df=NULL, gene.class=NULL, covs=NULL)

    samples = unique(dfc$sample)
    cat(sprintf("loading %d samples for cluster %s\n", length(samples), as.character(cluster.id)))
    for (sample in samples) {
        if (sample == "ref") next

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

        covs = read.cache(paste0(idir, "/cycle_covs_long_adjusted"))

        cc$sample = sample
        gene.df$sample = sample
        gene.df$cycle = gene.df$contig
        covs$sample = sample

        df$sample = sample
        df$id = sample

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

        restrict=function(table) { table[is.element(paste0(table$sample, "_", table$cycle), dfc$key),] }
        is.found = (dim(restrict(df))[1] > 0)
        if (!is.found) next
        ll$df = rbind(ll$df, restrict(df))
        ll$cc = rbind(ll$cc, restrict(cc))
        ll$gene.df = rbind(ll$gene.df, restrict(gene.df))
        ll$uniref = rbind(ll$uniref, restrict(uniref))
        ll$gene.class = rbind(ll$gene.class, restrict(gene.class))
        ll$covs = rbind(ll$covs, restrict(covs))
    }
    ll$df$label = dfc$label[match(paste(ll$df$sample,ll$df$cycle),paste(dfc$sample,dfc$cycle))]
    odir = paste0("figures/cplots/clusters_", .ver)

    dll = align.summary(dfc=dfc, dfx=dfx, dfs=dfs)
    keys = order.cycles(mat=dll$mat)
    dfc = dfc[match(keys, dfc$key),]

    if (is.null(ll$df))
        return (NULL)

    styles = c("detailed", "rep_clean", "rep_label", "clean", "label")
#    styles = c("rep_label")
    for (style in styles)
        plot.cluster.internal(cluster.id=cluster.id, prev.cluster=prev.cluster, dfc=dfc, ll=ll, odir=odir, plot.style=style)
}

plot.clusters=function()
{
    get.label=function(df)
    {
        fsamples = label.sample(df$sample)
        ifelse(fsamples != "ref", paste(fsamples, ":", df$cycle, sep=""), df$cycle)
    }

    get.ref.label=function(df)
    {
        dfr = read.delim("data/refs_anns.txt")
        ix = match(df$cycle, dfr$accession)
        dfr$sample_location[dfr$sample_location == ""] = "Location NA"
        dfr$collection_year[is.na(dfr$collection_year)] = "Year NA"
        ifelse(!is.na(ix), paste(dfr$host_name[ix], ", ", dfr$sample_location[ix], ", ", dfr$collection_year[ix], " (", dfr$accession[ix], ")", sep=""), "")
    }

    df.clusts = read.delim(paste0("output/clusters_", .ver, "/groups.txt"))
    dfg = read.delim(paste0("output/clusters_", .ver, "/cycles_grouped.txt"))
    refs = read.delim(paste0("output/clusters_", .ver, "/associated_refs_details.txt"))

    dfc = data.frame(cluster.id=dfg$id, sample=dfg$sample, cycle=dfg$cycle, length=dfg$length)
    dfc = rbind(dfc, data.frame(cluster.id=refs$cluster.id, sample="ref", cycle=refs$ACC_NUCCORE, length=refs$Length_NUCCORE))

    dfc$prev.cluster = dfg$prev.cluster[match(dfc$cluster.id, dfg$id)]
    dfc$base.label = get.label(dfc)

    dfc$ref.label = get.ref.label(dfc)
    dfc$is.ref = dfc$sample == "ref"
    dfc$label = ifelse(dfc$is.ref, dfc$ref.label, dfc$base.label)

    cat(sprintf("number of clusters: %d\n", dim(df.clusts)[1]))
    for (i in rev(1:dim(df.clusts)[1])) {
        prev.cluster = df.clusts$prev.cluster[i]
        cluster.id = df.clusts$id[i]
        plot.cluster(dfc=dfc, cluster.id=cluster.id, prev.cluster=prev.cluster)
    }
}

rl=function()
{
    source("cplot.r")
    source("pcluster.r")
}

ex=function()
{
    plot.clusters()
}

options(stringsAsFactors=F)

plot.cluster.internal=function(cluster, dfc, ll, plot.rep, odir)
{
    # internal radius of circles
    base.rad = 6

    odir.cluster = paste0(odir, "/", cluster)
    odir.legend = paste0(odir.cluster, "/legends")
    odir.circle = paste0(odir.cluster, "/circles")
    odir.rect = paste0(odir.cluster, "/rects")
    odir.sets = paste0(odir, "/sets")
    system(paste("mkdir -p", odir.cluster, odir.legend, odir.circle, odir.rect, odir.sets))

    class.col.list = list(mobile="darkgreen", plasmid="darkred", phage="orange")
    plot.legend(odir=odir.legend, cols=unlist(class.col.list), names=names(class.col.list), title="gene_class")

    ll$df$id = ll$df$sample
    profiles.cluster = list(
        cov.profile(height=3, title="cov", cov.table=ll$covs, cycle.table=ll$df, grid.nlines=4),
        empty.profile(height=0.5),
        gene.identity.profile(height=0.5, table=ll$gene.df, plot.trig=T, add.label=F, is.top=F),
        empty.profile(height=0.1),
        gene.uniref.count.profile(height=0.4, table=gene.df, odir.legend=odir.legend),
        empty.profile(height=0.5),
        gene.class.profile(height=0.5, table=ll$gene.df, cclass="mobile", col.list=class.col.list),
        empty.profile(height=0.1),
        gene.class.profile(height=0.5, table=ll$gene.df, cclass="plasmid", col.list=class.col.list),
        empty.profile(height=0.1),
        gene.class.profile(height=0.5, table=ll$gene.df, cclass="phage", col.list=class.col.list),
        empty.profile(height=0.5),
        align.profile(base.height=0.75, dfc=dfc, table.cycles=ll$df, table.align=ll$dfx, table.snps=ll$dfs),
        vlines.profile(table=ll$cc, field="cum_sum", lty=2))

    ofn = paste0(odir.sets, "/", cluster, ".pdf")

    if (!plot.rep) {
        cplot.multi(profiles=profiles.cluster, df=ll$df, cc=ll$cc, base.rad=base.rad,
                    plot.height.per.cycle=4, style="r", extra=1.1, ofn=ofn)
    } else {
        odir.circle = paste0(odir, "/reps/circles")
        odir.rect = paste0(odir, "/reps/rects")
        system(paste("mkdir -p", odir.circle, odir.rect))

        if (any(ll$df$sample == "cipro_clean"))
            df = ll$df[ll$df$sample == "cipro_clean",]
        else
            df = ll$df[which.max(ll$df$length),]

        df$id = cluster

        cplot.singles(profiles=profiles.cluster, df=df, cc=ll$cc, base.rad=base.rad,
                      odir.circle=odir.circle, odir.rect=odir.rect, circle.inch=7, rect.inch=4)
    }
}

plot.cluster=function(dfc, set.id, set.title, cluster, plot.rep)
{
    idir.cluster = "/relman01/home/nshalon/work/pipe/sour/cluster"
    idir.nucmer = "/relman01/home/nshalon/work/pipe/sour/nucmer"

    dfx = read.cache(paste0(idir.nucmer, "/", set.id, "/CLUSTER_", cluster, ".coords"), function(fn) { read.delim(fn, header=F) })
    dfs = read.cache(paste0(idir.nucmer, "/", set.id, "/CLUSTER_", cluster, ".snps"), function(fn) { safe.read.delim(fn, header=F) })

    dfc = dfc[dfc$cluster == cluster,]
    dfc$key = paste0(dfc$sample, "_", dfc$cycle)

    # are these columns correct?
    names(dfx) = c("start1", "end1", "start2", "end2", "length1", "length2", "identity", "name1", "name2")
    dfx$sample1 = sapply(strsplit(dfx$name1, "\\|"), function(x) { x[1] })
    dfx$cycle1 = sapply(strsplit(dfx$name1, "\\|"), function(x) { x[2] })
    dfx$sample2 = sapply(strsplit(dfx$name2, "\\|"), function(x) { x[1] })
    dfx$cycle2 = sapply(strsplit(dfx$name2, "\\|"), function(x) { x[2] })
    dfx = dfx[,-match(c("name1", "name2", "length1", "length2"), names(dfx))]
    dfx$key1 = paste0(dfx$sample1, "_", dfx$cycle1)
    dfx$key2 = paste0(dfx$sample2, "_", dfx$cycle2)

    if (!is.null(dfs)) {
        dfs = dfs[,c(1,4,2,3,11,12)]
        names(dfs) = c("base1", "base2", "nt1", "nt2", "name1", "name2")
        dfs$sample1 = sapply(strsplit(dfs$name1, "\\|"), function(x) { x[1] })
        dfs$cycle1 = sapply(strsplit(dfs$name1, "\\|"), function(x) { x[2] })
        dfs$sample2 = sapply(strsplit(dfs$name2, "\\|"), function(x) { x[1] })
        dfs$cycle2 = sapply(strsplit(dfs$name2, "\\|"), function(x) { x[2] })
        dfs$key1 = paste0(dfs$sample1, "_", dfs$cycle1)
        dfs$key2 = paste0(dfs$sample2, "_", dfs$cycle2)
    } else {
        dfs = setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("base1", "base2", "nt1", "nt2", "name1", "name2", "key1", "key2"))
    }

    dfx = dfx[is.element(dfx$key1, dfc$key) & is.element(dfx$key2, dfc$key),]
    dfs = dfs[is.element(dfs$key1, dfc$key) & is.element(dfs$key2, dfc$key),]

    ll = list(dfx=dfx, dfs=dfs, df=NULL, cc=NULL, gene.df=NULL, gene.class=NULL, covs=NULL)

    samples = unique(dfc$sample)
    cat(sprintf("loading %d samples for cluster %d\n", length(samples), cluster))
    for (sample in samples) {
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

        # uniref
        ix = match(gene.df$gene, uniref$gene)
        gene.df$desc = ifelse(!is.na(ix), uniref$prot_desc[ix], "no_hit")
        gene.df$taxa = ifelse(!is.na(ix), uniref$tax[ix], "no_hit")
        gene.df$identity = ifelse(!is.na(ix), uniref$identity[ix], 0)
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
    odir = paste0("figures/clusters/", set.title)
    if (!is.null(ll$df))
        plot.cluster.internal(cluster=cluster, dfc=dfc, ll=ll, odir=odir, plot.rep=plot.rep)
}

append.cycle.data=function(dfc)
{
    dfc$key = paste0(dfc$sample, "_", dfc$cycle)
    samples = unique(dfc$sample)
    rr = NULL
    for (sample in samples) {
        idir = paste0("/oak/stanford/groups/relman/users/nshalon/pipe/cyc_find_out4/", sample, "_0.0")
        df = read.cache(paste0(idir, "/dominant_cycles_filter/cycle_summary_classification"))
        df$key  = paste0(sample, "_", df$cycle)
        rr = rbind(rr, df)
    }
    ix = match(dfc$key, rr$key)
    if (any(is.na(ix)))
        stop("not all cycles found")
    dfc$length = rr$length[ix]
    dfc
}

plot.set=function(set.title="uniq_new_guts_cycles", plot.rep=T)
{
    idir.cluster = "/relman01/home/nshalon/work/pipe/sour/cluster"

    set.id = set.title
    if (set.title == "repeats")
        set.id = "ecologies_cycles"

    if (set.title == "shorts")
        set.id = "uniq_cycles"

    dfc = read.cache(paste0(idir.cluster, "/", set.id))
    dfc = append.cycle.data(dfc)

    if (set.title != "shorts")
        dfc = dfc[dfc$length >= 1000,]
    else
        dfc = dfc[dfc$length < 1000,]

    if (set.title == "repeats") {
        keep.samples = paste("gut", c("b", "d", "j"), sep="_")
        dfc = dfc[ is.element(dfc$sample, keep.samples),]
    }

    if (set.id == "uniq_cycles") {
        remove.samples = c("FP_B", "cipro_b", "gut_e")
        dfc = dfc[ !is.element(dfc$sample, remove.samples),]
    }

    tt = sort(table(dfc$cluster), decreasing=T)
    tt = tt[tt>1]
    cat(sprintf("Dataset table: %s\n", paste0(idir.cluster, "/", set.id)))
    cat(sprintf("number of clusters: %d\n", length(tt)))
    for (cluster in as.numeric(names(tt)))
        plot.cluster(dfc=dfc, cluster=cluster, set.id=set.id, set.title=set.title, plot.rep=plot.rep)
}

plot.sets=function()
{
    plot.rep = T
    title.ids = c("aab_long_cycles", "fp_long_cycles", "uniq_cycles", "repeats", "shorts")
    # title.ids = c("uniq_new_guts_cycles")
    for (set.title in title.ids)
        plot.set(set.title=set.title, plot.rep=plot.rep)
}

rlc=function()
{
    source("cplot.r")
    source("pcluster.r")
}

exc=function()
{
    plot.sets()
}

options(stringsAsFactors=F)

publish.tables=function()
{
    dfc = read.delim("data/sample_table.txt")
    rr.cycles = NULL
    rr.genes = NULL
    for (i in 1:dim(dfc)[1]) {
        sample = dfc$sample[i]
        sample.label = dfc$label[i]
        idir = paste0("/oak/stanford/groups/relman/users/nshalon/pipe/cyc_find_out4/", sample, "_0.0")
        cat(sprintf("loading sample %s from dir: %s\n", sample, idir))

        ifn = paste0(idir, "/dominant_cycles_filter/cycle_summary_classification")
        if (!file.exists(ifn)) {
            cat(sprintf("file missing: %s\nskipping sample: %s\n", ifn, sample))
            next
        }
        df = read.cache(ifn)
        df$sample = sample.label

        # !!! k was skipped
        cycs = df$cycle
        if (file.exists(paste0(idir, "/dominant_cycles_filter/prodigal/v1/gene.tab"))) {
            genes = read.cache(paste0(idir, "/dominant_cycles_filter/prodigal/v1/gene.tab"))
            uniref = read.cache(paste0(idir, "/dominant_cycles_filter/uniref/v1/2020_07/table_uniq_taxa"))
        } else {
            genes = read.cache(paste0(idir, "/dominant_cycles/prodigal/v1/gene.tab"))
            uniref = read.cache(paste0(idir, "/dominant_cycles/uniref/v1/2020_07/table_uniq_taxa"))
        }

        restrict=function(x) {
            if (any(is.element(x$cycle, cycs))) {
                x = x[is.element(x$cycle, cycs),]
                N = dim(x)[1]
                ix = match(c("cycle"), names(x))
                return (data.frame(sample=rep(sample.label, N), cycle=x$cycle, x[,-ix]))
            } else {
                return (NULL)
            }
        }

        uniref$cycle = genes$contig[match(uniref$gene, genes$gene)]
        uniref = restrict(uniref)

        genes$cycle = genes$contig
        genes = restrict(genes)

        ix = match(genes$gene, uniref$gene)
        genes$uniref = ifelse(!is.na(ix), uniref$uniref[ix], "no_hit")
        genes$desc = ifelse(!is.na(ix), uniref$prot_desc[ix], "no_hit")
        genes$taxa = ifelse(!is.na(ix), uniref$tax[ix], "no_hit")
        genes$identity = ifelse(!is.na(ix), uniref$identity[ix], 0)
        genes$coverage = ifelse(!is.na(ix), uniref$coverage[ix], 0)
        genes$uniref_count = ifelse(!is.na(ix), uniref$uniref_count[ix], 0)

        gene.class = read.cache(paste0(idir, "/dominant_cycles/gene_classification"))
        ix = match(genes$gene, gene.class$gene)
        genes$class = ifelse(!is.na(ix), gene.class$classification[ix], "none")
        genes$class.desc = ifelse(!is.na(ix), gene.class$description[ix], "none")

        # fix gene length bug
        genes$identity = round(genes$identity * genes$aa_length / (genes$aa_length-1),1)
        genes$length = genes$length - 3

        fields = c("sample", "cycle", "gene", "start", "end", "strand", "length", "uniref",
                   "desc", "taxa", "identity", "uniref_count", "class", "class.desc")
        rr.cycles = rbind(rr.cycles, df)
        rr.genes = rbind(rr.genes, genes[,fields])
    }

    # all cycles and genes
    odir = "output/all_cycles"
    system(paste("mkdir -p", odir))
    cat(sprintf("saving all cycles to: %s\n", odir))
    write.table(x=rr.cycles, file=paste0(odir, "/cycles.txt"), quote=F, row.names=F, sep="\t")
    write.table(x=rr.genes, file=paste0(odir, "/genes.txt"), quote=F, row.names=F, sep="\t")

    # circulating only
    df = read.delim("output/clusters_v2/cycles_grouped.txt")
    df$sample = label.sample(df$sample)
    df$key = paste(df$sample, df$cycle)
    restrict2=function(x) {
        x$key = paste(x$sample, x$cycle)
        x = x[is.element(x$key,df$key),]
        id = df$id[match(x$key,df$key)]
        rr = cbind(id, x[,-match("key",names(x))])
        rr[order(as.numeric(gsub("M", "", rr$id))),]
    }
    odir = "output/circulating_cycles"
    system(paste("mkdir -p", odir))
    cat(sprintf("saving circ cycles to: %s\n", odir))
    write.table(x=restrict2(rr.cycles), file=paste0(odir, "/cycles.txt"), quote=F, row.names=F, sep="\t")
    write.table(x=restrict2(rr.genes), file=paste0(odir, "/genes.txt"), quote=F, row.names=F, sep="\t")
}

publish.sequences=function()
{
    dfc = read.delim("data/sample_table.txt")

    cp=function(ifn, ofn) {
        if (!file.exists(ifn))
            stop(sprintf("file not found: %s", ifn))
        command = sprintf("cp %s %s", ifn, ofn)
        if (system(command) != 0)
            stop("cp failed")
    }

    # cycle fasta
    odir = "output/seq/cycles"
    system(paste("mkdir -p", odir))
    for (i in 1:dim(dfc)[1]) {
        sample = dfc$sample[i]
        sample.label = dfc$label[i]
        idir = paste0("/oak/stanford/groups/relman/users/nshalon/pipe/cyc_find_out4/", sample, "_0.0")
        ifn = paste0(idir, "/dominant_cycles_filter/cycles.fasta")
        ofn = paste0(odir, "/", sample.label, ".fasta")
        cp(ifn, ofn)
    }


    # gene faa and fna
    odir = "output/seq/genes"
    system(paste("mkdir -p", odir))
    for (i in 1:dim(dfc)[1]) {
        sample = dfc$sample[i]
        sample.label = dfc$label[i]
        idir = paste0("/oak/stanford/groups/relman/users/nshalon/pipe/cyc_find_out4/", sample, "_0.0")

        if (file.exists(paste0(idir, "/dominant_cycles/prodigal/v1/gene.tab"))) {
            idir.genes = paste0(idir, "/dominant_cycles/prodigal/v1")
        } else {
            idir.genes = paste0(idir, "/dominant_cycles_filter/prodigal/v1")
        }

        ifn.faa = paste0(idir.genes, "/genes_final.faa")
        ofn.faa = paste0(odir, "/", sample.label, ".faa")
        cp(ifn.faa, ofn.faa)

        ifn.fna = paste0(idir.genes, "/genes_final.fna")
        ofn.fna = paste0(odir, "/", sample.label, ".fna")
        cp(ifn.fna, ofn.fna)
    }
}

rl=function()
{
    source("cplot_utils.r")
    source("publish.r")
}

ex=function()
{
    publish.tables()
}

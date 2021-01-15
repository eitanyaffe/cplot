options(stringsAsFactors=F)
source("align_utils.r")

######################################################################################################
# util functions
######################################################################################################

append.cycle.data=function(dfc)
{

    dfc$key = paste0(dfc$sample, "_", dfc$cycle)
    samples = unique(dfc$sample)
    rr = NULL
    for (sample in samples) {
        if (sample == "ref") next
        idir = paste0("/oak/stanford/groups/relman/users/nshalon/pipe/cyc_find_out4/", sample, "_0.0")
        df = read.cache(paste0(idir, "/dominant_cycles_filter/cycle_summary_classification"))
        df$key  = paste0(sample, "_", df$cycle)
        rr = rbind(rr, df)
    }
    ix = match(dfc$key, rr$key)
    dfc$length = ifelse(!is.na(ix), rr$length[ix], 0)
    dfc
}


######################################################################################################
# summarize all vs nucmer results nucmer
######################################################################################################

basic.summary=function()
{
    dfx = parse.nucmer.coords(read.delim("/relman01/home/nshalon/work/pipe/sour/nucmer/all_all/all_all.coords", header=F))
    dfs = parse.nucmer.snps(read.delim("/relman01/home/nshalon/work/pipe/sour/nucmer/all_all/all_all.snps", header=F))

    dd = dfx[dfx$key1 == dfx$key2,]
    dd = dd[!duplicated(dd$key1),]
    dfc = data.frame(sample=dd$sample1, cycle=dd$cycle1, key=dd$key1)
    dfc = append.cycle.data(dfc=dfc)
    dfc = dfc[dfc$length >= 1000,]

    dll = align.summary(dfc=dfc, dfx=dfx, dfs=dfs)
    rr = dll$table
    rr = rr[rr$key1 != rr$key2 & rr$overlap.fraction > 0,]

    write.table(x=rr, file="output/nucmer_summary.txt", quote=F, row.names=F, sep="\t")
    write.table(x=dfc, file="output/cycles_basic.txt", quote=F, row.names=F, sep="\t")
}

ref.summary=function()
{
    dfx = parse.nucmer.coords(read.delim("/relman01/home/nshalon/work/pipe/sour/nucmer/circ_refs/circ_refs.coords",
                                         header=F))
    dfs = parse.nucmer.snps(read.delim("/relman01/home/nshalon/work/pipe/sour/nucmer/circ_refs/circ_refs.snps",
                                       header=F))

    dd = dfx[dfx$key1 == dfx$key2,]
    dd = dd[!duplicated(dd$key1),]
    dfc = data.frame(sample=dd$sample1, cycle=dd$cycle1, key=dd$key1)
    dfc = append.cycle.data(dfc=dfc)

    # add ref lengths
    raw = read.delim("data/plsdb_raw.txt")
    iy = match(dfc$key, paste("ref", raw$ACC_NUCCORE, sep="_"))
    if (any(is.na(iy) & dfc$length == 0))
        stop("missing ref in plsdb table")
    dfc$length = ifelse(dfc$length == 0, raw$Length_NUCCORE[iy], dfc$length)

    dll = align.summary(dfc=dfc, dfx=dfx, dfs=dfs)
    rr = dll$table
    rr = rr[rr$key1 != rr$key2 & rr$overlap.fraction > 0,]

    write.table(x=rr, file="output/refs_nucmer_summary.txt", quote=F, row.names=F, sep="\t")
    write.table(x=dfc, file="output/refs_cycles.txt", quote=F, row.names=F, sep="\t")
}


rl=function()
{
    source("cplot.r")
    source("nucmer_summary.r")
}

ex=function()
{
    basic.summary()
    ref.summary()
}

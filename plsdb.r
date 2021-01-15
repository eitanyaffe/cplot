options(stringsAsFactors=F)
.ver = "v2"

get.id=function(dfc, key)
{
    ix = match(key, dfc$key)
    ifelse(!is.na(ix), dfc$id[ix], "na")
}

associate.refs=function()
{
    max.dist = 0.1

    wdir = paste0("output/clusters_", .ver)

    # dfm = read.delim("/relman01/home/nshalon/work/pipe/sour//strain_species/circ_refs_distances")

    dfm = read.delim(paste0("output/refs_nucmer_summary.txt"))
    dfc = read.delim(paste0(wdir, "/cycles_grouped.txt"))

    dfm$cluster.id1 = get.id(dfc=dfc, key=dfm$key1)
    dfm$cluster.id2 = get.id(dfc=dfc, key=dfm$key2)
    dfm$is.ref1 = grepl("^ref", dfm$key1)
    dfm$is.ref2 = grepl("^ref", dfm$key2)

    xx1 = dfm[dfm$cluster.id1 != "na" & dfm$is.ref2,]
    xx2 = dfm[dfm$cluster.id2 != "na" & dfm$is.ref1,]
    dd = data.frame(
        cluster.id=c(xx1$cluster.id1, xx2$cluster.id2),
        nearest.member=c(xx1$key1, xx2$key2),
        ref=c(xx1$key2, xx2$key1),
        dist=c(xx1$dist, xx2$dist))
    dd = dd[dd$dist <= max.dist,]

    dd$ref = sub("ref_", "", dd$ref)

    rr = NULL
    for (cluster.id in unique(dfc$id)) {
        if (!any(dd$cluster.id == cluster.id)) next
        ddx = dd[dd$cluster.id == cluster.id,]
        ddx = ddx[!duplicated(ddx$ref),]
        rr = rbind(rr, ddx)
    }
    rr$prev.cluster = dfc$prev.cluster[match(rr$cluster.id, dfc$id)]

    cat(sprintf("working directory: %s\n", wdir))
    write.table(x=rr, file=paste0(wdir, "/associated_refs_",max.dist, ".txt"), quote=F, row.names=F, sep="\t")

    # annotate
    df = read.delim("data/plsdb_raw.txt")
    df = df[,-match(c("query", "distance", "pvalue", "shared_hashes", "recid"), names(df))]
    ix = match(rr$ref, df$ACC_NUCCORE)
    df.rr = cbind(rr[,-match("ref", names(rr))], df[ix,])
    write.table(x=df.rr, file=paste0(wdir, "/associated_refs_details_",max.dist, ".txt"), quote=F, row.names=F, sep="\t")

    # compare to previous results
    # zz = read.delim("data/alignment_stats.txt")
    # zz = zz[zz$species != "phix",]
    # zz$cluster.id = dfc$cluster.id[match(zz$order, dfc$cycle)]
}


rl=function()
{
    source("plsdb.r")
}

ex=function()
{
    associate.refs()
}

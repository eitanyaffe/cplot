options(stringsAsFactors=F)
make.table1=function()
{
    dfi = read.delim("data/table_1_v2.txt")
    df = read.delim("output/clusters_v2/groups.txt")
    dfr = read.delim("output/clusters_v2/associated_refs_details.txt")

    ix = match(df$prev.cluster, dfi$Internal.Cluster.Name)
    df$length = ifelse(dfi$Min.Length[ix] == dfi$Min.Length.1[ix], paste0(dfi$Min.Length[ix],"bp"), paste0(dfi$Min.Length[ix],"-", dfi$Min.Length.1[ix], "bp"))
    df$environment = dfi$Environment[ix]
    df$identity = round(100*(1-dfi$Genomic.Distance[ix]),2)
    df$aligned.fraction = round(100*dfi$Aligned.Fraction[ix],2)
    df$aligned.identity = round(100*dfi$Aligned.Identity[ix],2)
    df$snps = round(dfi$SNPS[ix])
    df$AP = dfi$Mean.AP[ix]
    df$gene.count = dfi$X..Genes[ix]
    df$class = dfi$Class[ix]

    df$ref95.count = 0
    df$ref99.9.count = 0
    df$nearest.ref.identity = 0
    df$nearest.ref = ""
    taxas = c("taxon_species_name", "taxon_genus_name", "taxon_family_name", "taxon_order_name", "taxon_class_name", "taxon_phylum_name")
    for (taxa in taxas)
        df[,taxa] = ""

    dfr$identity = 100*(1-dfr$dist)
    for (id in unique(dfr$cluster.id)) {
        dfri = dfr[dfr$cluster.id == id,]
        iy = which.max(dfri$identity)
        ix = match(id, df$id)

        df$ref95.count[ix] = sum(dfri$identity>=95)
        df$ref99.9.count[ix] = sum(dfri$identity>=99.9)
        df$nearest.ref.identity[ix] = round(dfri$identity[iy],2)
        df$nearest.ref[ix] = dfri$ACC_NUCCORE[iy]

        taxas = c("taxon_species_name", "taxon_genus_name", "taxon_family_name", "taxon_order_name", "taxon_class_name", "taxon_phylum_name")
        for (taxa in taxas)
            df[ix,taxa] = paste(unique(dfri[,taxa]), sep=" ", collapse="")
    }
    write.table(x=df, file="output/table1.txt", quote=F, row.names=F, sep="\t")
}

rl=function()
{
    source("table1.r")
}

ex=function()
{
    make.table1()
}

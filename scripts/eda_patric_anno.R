## selected gene functions
# sialidase
func.sialidase <- fasta.superkey[grep("Sialidase", fasta.superkey$full_id, ignore.case = TRUE),]
gv.no.sialidase <- unique(fasta.superkey$strain_id[!(fasta.superkey$strain_id %in% func.sialidase$strain_id)])
clus.sialidase.id <- unique(func.sialidase$clust_id)

# check for unannotated sialidase in same clusters
# clus.sialidase <- fasta.superkey[fasta.superkey$clust_id %in% clus.sialidase.id,] # it's the same

# split by sialidase cluster
func.sialidase.A <- func.sialidase[func.sialidase$clust_id=="clust_0917",]
func.sialidase.B <- func.sialidase[func.sialidase$clust_id=="clust_1160",]
func.sialidase.C <- func.sialidase[func.sialidase$clust_id=="clust_1950",]
func.sialidase.D <- func.sialidase[func.sialidase$clust_id=="clust_2656",]

# amylase
func.amylase <- fasta.superkey[grep("amylase", fasta.superkey$full_id, ignore.case = TRUE),]

# vaginolysin --> none
# func.vaginolysin <- fasta.superkey[grep("vaginolysin", fasta.superkey$full_id, ignore.case = TRUE),]

# lysin
func.lysin <- fasta.superkey[grep("lysin", fasta.superkey$full_id, ignore.case = TRUE),]
func.lysin <- func.lysin[-grep("Lysine", func.lysin$full_id),]

# cytolysin
func.cytolysin <- fasta.superkey[grep("cytolysin", fasta.superkey$full_id, ignore.case = TRUE),]

# biofilm --> none
# func.biofilm <- fasta.superkey[grep("biofilm", fasta.superkey$full_id, ignore.case = TRUE),]

# toxin
func.toxin <- fasta.superkey[grep("toxin", fasta.superkey$full_id, ignore.case = TRUE),]

# plasmid
func.plasmid <- fasta.superkey[grep("plasmid", fasta.superkey$full_id, ignore.case = TRUE),]

# mucin --> none
# func.mucin <- fasta.superkey[grep("mucin", fasta.superkey$full_id, ignore.case = TRUE),]

# neuraminidase --> none
# func.neuraminidase <- fasta.superkey[grep("neuraminidase", fasta.superkey$full_id, ignore.case = TRUE),]


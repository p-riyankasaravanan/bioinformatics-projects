# ============================================================
# Vaginal Microbiome 16S rRNA Analysis
# Pilot DADA2 workflow
# ============================================================

library(dada2)

# ---- Project paths ----

project_path <- "E:/Vaginal Microbiome Project"

raw_path <- file.path(
  project_path,
  "data",
  "raw_fastq"
)

filtered_path <- file.path(
  project_path,
  "data",
  "filtered"
)

results_path <- file.path(
  project_path,
  "results"
)

taxonomy_path <- file.path(
  project_path,
  "taxonomy"
)

dir.create(filtered_path, showWarnings = FALSE, recursive = TRUE)
dir.create(results_path, showWarnings = FALSE, recursive = TRUE)

# ---- Pilot sample ----

sample_id <- "SRR12976656"

fnF <- file.path(
  raw_path,
  paste0(sample_id, "_1.fastq")
)

fnR <- file.path(
  raw_path,
  paste0(sample_id, "_2.fastq")
)

filtF <- file.path(
  filtered_path,
  paste0(sample_id, "_1_filtered.fastq.gz")
)

filtR <- file.path(
  filtered_path,
  paste0(sample_id, "_2_filtered.fastq.gz")
)

# ---- Quality filtering ----

filt <- filterAndTrim(
  fnF, filtF,
  fnR, filtR,
  truncLen = c(280, 230),
  maxEE = c(2, 2),
  maxN = 0,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = FALSE
)

# ---- Learn sequencing error rates ----

errF <- learnErrors(
  filtF,
  multithread = FALSE
)

errR <- learnErrors(
  filtR,
  multithread = FALSE
)

# ---- Infer ASVs ----

dadaF <- dada(
  filtF,
  err = errF,
  multithread = FALSE
)

dadaR <- dada(
  filtR,
  err = errR,
  multithread = FALSE
)

# ---- Merge paired reads ----

merger <- mergePairs(
  dadaF,
  filtF,
  dadaR,
  filtR,
  verbose = TRUE
)

# ---- Construct sequence table ----

seqtab <- makeSequenceTable(merger)

# ---- Remove chimeric sequences ----

seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = FALSE,
  verbose = TRUE
)

# Save non-chimeric ASV table

saveRDS(
  seqtab.nochim,
  file.path(results_path, "seqtab_nochim.rds")
)

# ---- Taxonomic assignment ----

tax_train <- file.path(
  taxonomy_path,
  "silva_nr99_v138.2_toGenus_trainset.fa.gz"
)

taxa <- assignTaxonomy(
  seqtab.nochim,
  tax_train,
  multithread = FALSE
)

saveRDS(
  taxa,
  file.path(results_path, "taxonomy.rds")
)

# ---- Create ASV taxonomy and abundance table ----

asv_counts <- as.numeric(seqtab.nochim[1, ])

asv_table <- data.frame(
  ASV = paste0("ASV_", seq_along(asv_counts)),
  Reads = asv_counts,
  taxa,
  check.names = FALSE
)

asv_table$Relative_abundance <-
  100 * asv_table$Reads / sum(asv_table$Reads)

asv_table <- asv_table[
  order(asv_table$Reads, decreasing = TRUE),
]

write.csv(
  asv_table,
  file.path(
    results_path,
    "pilot_ASV_taxonomy_abundance.csv"
  ),
  row.names = FALSE
)

# ---- Taxonomic richness ----

taxonomic_richness <- data.frame(
  Rank = c(
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus"
  ),
  Richness = c(
    length(unique(na.omit(taxa[, "Phylum"]))),
    length(unique(na.omit(taxa[, "Class"]))),
    length(unique(na.omit(taxa[, "Order"]))),
    length(unique(na.omit(taxa[, "Family"]))),
    length(unique(na.omit(taxa[, "Genus"])))
  )
)

write.csv(
  taxonomic_richness,
  file.path(
    results_path,
    "taxonomic_richness.csv"
  ),
  row.names = FALSE
)
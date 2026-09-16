# ============================================================
# Vaginal Microbiome 16S rRNA Analysis
# Pilot taxonomy, abundance and diversity analysis
# ============================================================

library(ggplot2)

# ---- Project paths ----

project_path <- "E:/Vaginal Microbiome Project"

results_path <- file.path(
  project_path,
  "results"
)

figures_path <- file.path(
  project_path,
  "figures"
)

dir.create(
  figures_path,
  showWarnings = FALSE,
  recursive = TRUE
)

# ---- Load pilot results ----

seqtab.nochim <- readRDS(
  file.path(results_path, "seqtab_nochim.rds")
)

taxa <- readRDS(
  file.path(results_path, "taxonomy.rds")
)

# ---- ASV abundance table ----

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

# ---- Phylum abundance ----

phylum_counts <- tapply(
  asv_table$Reads,
  asv_table$Phylum,
  sum,
  na.rm = TRUE
)

phylum_table <- data.frame(
  Phylum = names(phylum_counts),
  Reads = as.numeric(phylum_counts)
)

phylum_table$Relative_abundance <-
  100 * phylum_table$Reads / sum(phylum_table$Reads)

phylum_table <- phylum_table[
  order(phylum_table$Relative_abundance, decreasing = TRUE),
]

write.csv(
  phylum_table,
  file.path(results_path, "phylum_abundance.csv"),
  row.names = FALSE
)

# ---- Genus abundance ----

genus <- as.character(asv_table$Genus)

genus[is.na(genus) | genus == ""] <- "Unclassified"

genus_counts <- tapply(
  asv_table$Reads,
  genus,
  sum
)

genus_table <- data.frame(
  Genus = names(genus_counts),
  Reads = as.numeric(genus_counts)
)

genus_table$Relative_abundance <-
  100 * genus_table$Reads / sum(genus_table$Reads)

genus_table <- genus_table[
  order(genus_table$Relative_abundance, decreasing = TRUE),
]

write.csv(
  genus_table,
  file.path(results_path, "genus_abundance.csv"),
  row.names = FALSE
)

# ---- Family abundance ----

family <- as.character(asv_table$Family)

family[is.na(family) | family == ""] <- "Unclassified"

family_counts <- tapply(
  asv_table$Reads,
  family,
  sum
)

family_table <- data.frame(
  Family = names(family_counts),
  Reads = as.numeric(family_counts)
)

family_table$Relative_abundance <-
  100 * family_table$Reads / sum(family_table$Reads)

family_table <- family_table[
  order(family_table$Relative_abundance, decreasing = TRUE),
]

write.csv(
  family_table,
  file.path(results_path, "family_abundance.csv"),
  row.names = FALSE
)

# ---- ASV dominance ----

dominance <- data.frame(
  ASVs = c(1, 5, 10, 20),
  Relative_abundance = c(
    sum(head(asv_table$Relative_abundance, 1)),
    sum(head(asv_table$Relative_abundance, 5)),
    sum(head(asv_table$Relative_abundance, 10)),
    sum(head(asv_table$Relative_abundance, 20))
  )
)

write.csv(
  dominance,
  file.path(results_path, "pilot_ASV_dominance.csv"),
  row.names = FALSE
)

# ---- Alpha diversity ----

counts <- asv_table$Reads
counts <- counts[counts > 0]

p <- counts / sum(counts)

shannon <- -sum(p * log(p))

simpson <- 1 - sum(p^2)

alpha_diversity <- data.frame(
  Sample = "SRR12976656",
  Observed_ASVs = length(counts),
  Shannon = shannon,
  Simpson = simpson
)

write.csv(
  alpha_diversity,
  file.path(results_path, "pilot_alpha_diversity.csv"),
  row.names = FALSE
)

# ---- Pastel palette ----

pastel_palette <- c(
  "#B8E0D2", "#F7C8E0", "#CDB4DB",
  "#A9D6E5", "#FFE5B4", "#CDEAC0",
  "#BDE0FE", "#F4B6C2", "#D8E2DC",
  "#E8C7B8", "#C6D8E4"
)

# ---- Phylum plot ----

p_phylum <- ggplot(
  phylum_table,
  aes(
    x = reorder(Phylum, Relative_abundance),
    y = Relative_abundance,
    fill = Phylum
  )
) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = pastel_palette
  ) +
  labs(
    title = "Bacterial composition at phylum level",
    subtitle = "SRR12976656 | Pilot sample",
    x = NULL,
    y = "Relative abundance (%)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

ggsave(
  file.path(figures_path, "pilot_phylum_composition.png"),
  p_phylum,
  width = 8,
  height = 6,
  dpi = 300
)

# ---- Top 15 genus plot ----

top15_genus <- head(genus_table, 15)

other_genus <- data.frame(
  Genus = "Other",
  Reads = sum(genus_table$Reads) -
    sum(top15_genus$Reads)
)

other_genus$Relative_abundance <-
  100 * other_genus$Reads / sum(genus_table$Reads)

genus_plot <- rbind(
  top15_genus,
  other_genus
)

p_genus <- ggplot(
  genus_plot,
  aes(
    x = reorder(Genus, Relative_abundance),
    y = Relative_abundance,
    fill = Genus
  )
) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = pastel_palette
  ) +
  labs(
    title = "Dominant bacterial genera",
    subtitle = "SRR12976656 | Pilot sample",
    x = NULL,
    y = "Relative abundance (%)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

ggsave(
  file.path(figures_path, "pilot_genus_composition.png"),
  p_genus,
  width = 8,
  height = 7,
  dpi = 300
)

# ---- Top 10 family plot ----

top10_family <- head(family_table, 10)

other_family <- data.frame(
  Family = "Other",
  Reads = sum(family_table$Reads) -
    sum(top10_family$Reads)
)

other_family$Relative_abundance <-
  100 * other_family$Reads / sum(family_table$Reads)

family_plot <- rbind(
  top10_family,
  other_family
)

p_family <- ggplot(
  family_plot,
  aes(
    x = reorder(Family, Relative_abundance),
    y = Relative_abundance,
    fill = Family
  )
) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = pastel_palette
  ) +
  labs(
    title = "Dominant bacterial families",
    subtitle = "SRR12976656 | Pilot sample",
    x = NULL,
    y = "Relative abundance (%)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

ggsave(
  file.path(figures_path, "pilot_family_composition.png"),
  p_family,
  width = 8,
  height = 6,
  dpi = 300
)
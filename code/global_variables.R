# ==============================================================================
# Global Variables
# ==============================================================================

# --- Statistical thresholds ---------------------------------------------------
p_val <- 0.05
p_val_strict <- 0.01
FC <- 0.25 # Log2 fold-change threshold for filtering

# --- Output file extensions ---------------------------------------------------
tiff_extension <- ".tiff"
pdf_extension <- ".pdf"

# --- Species / organism -------------------------------------------------------
species <- 9606 # NCBI Taxonomy ID for Homo sapiens #9606 for Human, 10090 for mouse, 7955 for zebrafish
organism <- "org.Hs.eg.db" # Bioconductor annotation package #"org.Dr.eg.db". "org.Hs.eg.db"
kegg_organism <- "hsa" # KEGG organism code for Homo sapiens #dre for Danio rerio, hsa for Homo sapiens, mmu for Mus musculus
keyType <- "SYMBOL"
KEGGkeyType <- "kegg"


# --- QC thresholds (standard human PBMC values) --------------------------------
# Based on: https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
# and standard 10X PBMC recommendations (PBMCs are cytoplasm-rich → higher MT tolerance)
QC_MIN_FEATURES <- 200 # Minimum genes detected per cell
QC_MAX_FEATURES <- 5000 # Maximum genes per cell (high = likely doublet)
QC_MIN_COUNTS <- 500 # Minimum UMI counts per cell
QC_MAX_COUNTS <- 25000 # Maximum UMI counts per cell (high = likely doublet)
QC_MAX_MT <- 20 # Maximum % mitochondrial reads
# NOTE: PBMCs have cytoplasm → higher baseline MT than nuclei
# Paper used snRNA from brain (1% MT), here 20% is standard
QC_MIN_COMPLEXITY <- 0.80 # Minimum log10(genes / UMI) novelty score
QC_MIN_RIBO <- 0 # Minimum % ribosomal reads (no hard lower bound)
QC_MAX_RIBO <- 60 # Maximum % ribosomal reads


# --- Integration and clustering parameters ------------------------------------
N_INTEGRATION_FEATURES <- 2000 # Lower than default 3000 to save RAM
N_PCS_MAX <- 50 # Maximum PCs for ElbowPlot
N_PCS_INTEGRATION <- 20 # PCs for RPCA integration
CLUSTERING_DIMS <- c(10, 20, 30) # Dims to test
CLUSTERING_RESOLUTIONS <- c(0.1, 0.2, 0.4, 0.8) # Resolutions to test

# --- Memory / hardware settings (i7-7560U: 2 physical / 4 logical cores, 16 GB RAM) --
# Physical cores = 2; sequential plan avoids fork-overhead and protects RAM
# Swap available: 16 GB — can absorb moderate overflow but slows analysis
PARALLEL_WORKERS <- 1 # Sequential: safer than parallel on 2-core laptop
FUTURE_GLOBALS_MAX_MB <- 6000 # 6 GB global size limit for {future} (conservative)
options(future.globals.maxSize = FUTURE_GLOBALS_MAX_MB * 1024^2)

# --- Checkpoint naming --------------------------------------------------------
# Checkpoints are stored as: output_path/RData/checkpoint_<NAME>.rds
# This allows recovery from any step if the session crashes
CHECKPOINT_PREFIX <- "checkpoint_"

# --- Filtering patterns --------------------------------------------------------
mito_pattern <- "^mt-" # Mitochondrial genes start with MT- in human and mouse; in mouse they start with mt-
ribo_pattern <- "^Rp[ls]" # Ribosomal protein genes start with RPS or RPL; in mouse they start with Rps or Rpl


# --- Directory creation -------------------------------------------------------
create_directories <- function(base_path) {
  dirs <- c(
    "",
    "tables/",
    "figures/",
    "RData/",
    "RData/sct_objects/", # Per-library SCT checkpoints
    "RData/doublets/", # Per-library DoubletFinder results
    "figures/QC/",
    "figures/integration/",
    "figures/clustering/",
    "figures/annotation/",
    "figures/DE/",
    "figures/proportions/"
  )
  for (subdir in dirs) {
    path <- paste0(base_path, subdir)
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
      message("  [DIR] Created: ", path)
    }
  }
  message("[OK] Directory structure ready at: ", base_path)
}

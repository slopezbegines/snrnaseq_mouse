# ==============================================================================
# Package Loading — GSE194315 PBMC CITE-seq Analysis
# All packages required for the full pipeline
# ==============================================================================

# Check pak package
if ("pak" %in% installed.packages()) {
  library("pak")
} else {
  install.packages("pak")
}

# Helper: install if missing -----------------------------------------------
install_if_missing <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("  [INSTALL] Installing: ", pkg)
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        pak::pkg_install("BiocManager")
      }
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      pak::pkg_install(pkg)
    }
  }
}

# Helper: install from GitHub if missing (DoubletFinder) --------------------
install_github_if_missing <- function(pkg_name, github_repo) {
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    message("  [INSTALL] Installing from GitHub: ", github_repo)
    if (!requireNamespace("remotes", quietly = TRUE)) {
      pak::pkg_install("remotes")
    }
    remotes::install_github(github_repo, upgrade = "never")
  }
}

# =============================================================================
# CRAN packages
# =============================================================================
message("[PACKAGES] Loading CRAN packages...")

cran_packages <- c(
  # Core data manipulation
  "dplyr", "tidyr", "tibble", "stringr", "purrr", "data.table",
  # Visualization
  "ggplot2", "patchwork", "ggrepel", "scales",
  # I/O
  "readxl", "writexl",
  # Utilities
  "tictoc", "remotes", "future", "future.apply",
  # Clustering diagnostics
  "clustree"
)

for (pkg in cran_packages) install_if_missing(pkg)
invisible(lapply(cran_packages, library, character.only = TRUE, warn.conflicts = TRUE, quietly = TRUE))

# =============================================================================
# Bioconductor packages
# =============================================================================
message("[PACKAGES] Loading Bioconductor packages...")

bioc_packages <- c(
  # Core single-cell
  "Seurat", # Main scRNA-seq analysis framework
  "SeuratObject", # Seurat data structures
  # Cell type annotation
  "SingleR", # Automated cell type annotation
  "celldex", # Reference datasets for SingleR
  # Gene annotations (Human)
  "org.Hs.eg.db",
  # Bioconductor infrastructure
  "BiocParallel",
  "SummarizedExperiment",
  "S4Vectors",
  # Optional: faster SCTransform backend
  "glmGamPoi",
  # Optional: faster DE testing
  "MAST",
  # Visualization:
  "EnhancedVolcano",
  # Doublet detection:
  "scDblFinder"
)

for (pkg in bioc_packages) install_if_missing(pkg, bioc = TRUE)
invisible(lapply(bioc_packages, library, character.only = TRUE, warn.conflicts = TRUE, quietly = TRUE))

# =============================================================================
# GitHub packages
# =============================================================================
message("[PACKAGES] Checking GitHub packages...")

# DoubletFinder — doublet detection (not on CRAN/Bioconductor)
install_github_if_missing("DoubletFinder", "chris-mcginnis-ucsf/DoubletFinder")
library(DoubletFinder)

# SeuratDisk — optional, for H5Seurat save/load
# Uncomment if needed:
install_github_if_missing("SeuratDisk", "mojaveazure/seurat-disk")
if (requireNamespace("SeuratDisk", quietly = TRUE)) library(SeuratDisk)

# SeuratData —
# Uncomment if needed:
install_github_if_missing("SeuratData", "satijalab/seurat-data")
if (requireNamespace("SeuratData", quietly = TRUE)) library(SeuratData)

# SeuratWrappers —
# Uncomment if needed:
install_github_if_missing("SeuratWrappers", "satijalab/seurat-wrappers")
if (requireNamespace("SeuratWrappers", quietly = TRUE)) library(SeuratWrappers)


# =============================================================================
# Session info
# =============================================================================
message("\n[PACKAGES] All packages loaded. Session info:")
cat("R version:", R.version$version.string, "\n")
cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")
cat("SingleR version:", as.character(packageVersion("SingleR")), "\n")
cat("DoubletFinder loaded:", requireNamespace("DoubletFinder", quietly = TRUE), "\n")
cat("glmGamPoi available:", requireNamespace("glmGamPoi", quietly = TRUE), "\n")

# Clean enviroment
rm(install_if_missing, install_github_if_missing, cran_packages, bioc_packages)

# renv::snapshot(type = "all") # Uncomment to save package versions to renv.lock

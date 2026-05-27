# ==============================================================================
# Package Loading — snRNAseq Mouse Pipeline
# Uses pak for all installations: parallel, handles Bioc + GitHub, fast.
# First-time use: run `sudo apt-get install -y libuv1-dev` before this script.
# ==============================================================================

# Bootstrap pak itself --------------------------------------------------------
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak", repos = sprintf(
    "https://r-lib.github.io/p/pak/stable/%s/%s/%s",
    .Platform$pkgType, R.Version()$os, R.Version()$arch
  ))
}
library(pak)

# Print system requirements for any missing packages (informational) ----------
check_sysreqs <- function(pkgs) {
  tryCatch(
    {
      reqs <- pak::pkg_sysreqs(pkgs)
      if (length(reqs$packages) > 0) {
        message("[SYSREQS] Missing system packages detected:")
        message("  Run: sudo apt-get install -y ", paste(reqs$packages, collapse = " "))
      }
    },
    error = function(e) invisible(NULL)
  )
}

# Helper: install + load, skipping already-installed -------------------------
# pkgs can use pak prefixes (bioc::Pkg, user/repo) for install;
# library() needs bare names, so strip everything up to and including "::".
load_pkgs <- function(pkgs) {
  bare <- sub("^[^:]+::", "", pkgs)
  missing_idx <- !vapply(bare, requireNamespace, logical(1), quietly = TRUE)
  if (any(missing_idx)) {
    message("[INSTALL] Installing: ", paste(pkgs[missing_idx], collapse = ", "))
    pak::pkg_install(pkgs[missing_idx], ask = FALSE, upgrade = FALSE)
  }
  invisible(lapply(bare, function(p) {
    suppressPackageStartupMessages(suppressWarnings(
      library(p, character.only = TRUE, warn.conflicts = FALSE, quietly = TRUE)
    ))
  }))
}

# =============================================================================
# CRAN packages  (includes Seurat + SeuratObject — they live on CRAN, not Bioc)
# =============================================================================
message("[PACKAGES] Loading CRAN packages...")

cran_packages <- c(
  # Single-cell core (CRAN)
  "Seurat", "SeuratObject", "Signac",
  # Core data manipulation
  "dplyr", "tidyr", "tibble", "stringr", "purrr", "data.table",
  # Visualization
  "ggplot2", "patchwork", "ggrepel", "scales", "scCustomize",
  # I/O
  "readxl", "writexl",
  # Utilities
  "tictoc", "remotes", "future", "future.apply",
  # Clustering diagnostics
  "clustree",
  # Counting time
  "tictoc",
  # Harmony integration
  "harmony",
  #
  "scrapper"
)

load_pkgs(cran_packages)

# =============================================================================
# Bioconductor packages  (bioc:: prefix — these are genuine Bioc packages)
# =============================================================================
message("[PACKAGES] Loading Bioconductor packages...")

bioc_packages <- paste0("bioc::", c(
  "SingleR", # Automated cell type annotation
  "celldex", # Reference datasets for SingleR
  "org.Hs.eg.db", # Gene annotations (Human)
  "BiocParallel",
  "SummarizedExperiment",
  "S4Vectors",
  "glmGamPoi", # Faster SCTransform backend
  "MAST", # Faster DE testing
  "EnhancedVolcano",
  "scDblFinder", # Doublet detection
  "clusterProfiler", # GO enrichment
  "STRINGdb", # STRING PPI database
  "enrichplot", # GO enrichment visualization
  "org.Hs.eg.db", # Gene annotations (Human)
  "org.Mm.eg.db" # Gene annotations (Mouse)
))
#
load_pkgs(bioc_packages)

# =============================================================================
# GitHub / r-universe packages
# names  = bare package name used by requireNamespace() and library()
# values = pak installation spec (user/repo or URL::pkg)
# =============================================================================
message("[PACKAGES] Loading GitHub packages...")

gh_packages <- c(
  "DoubletFinder"  = "chris-mcginnis-ucsf/DoubletFinder",
  "SeuratDisk"     = "mojaveazure/seurat-disk",
  "SeuratData"     = "satijalab/seurat-data",
  "SeuratWrappers" = "satijalab/seurat-wrappers",
  "BPCells"        = "bnprks/BPCells",
  "presto"         = "immunogenomics/presto"
)

for (pkg_name in names(gh_packages)) {
  pkg_spec <- gh_packages[[pkg_name]]
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    message("[INSTALL] Installing from GitHub: ", pkg_spec)
    pak::pkg_install(pkg_spec, ask = FALSE, upgrade = FALSE)
  }
  if (requireNamespace(pkg_name, quietly = TRUE)) {
    suppressPackageStartupMessages(suppressWarnings(
      library(pkg_name, character.only = TRUE, warn.conflicts = FALSE, quietly = TRUE)
    ))
  }
}
#  remotes::install_github("satijalab/azimuth", ref = "master")
# devtools::install_github("satijalab/AzimuthAPI")
# =============================================================================
# Session info
# =============================================================================
message("\n[PACKAGES] All packages loaded.")
cat("R version:", R.version$version.string, "\n")
cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")
cat("SingleR version:", as.character(packageVersion("SingleR")), "\n")
cat("DoubletFinder loaded:", requireNamespace("DoubletFinder", quietly = TRUE), "\n")
cat("glmGamPoi available:", requireNamespace("glmGamPoi", quietly = TRUE), "\n")
cat("SeuratDisk available:", requireNamespace("SeuratDisk", quietly = TRUE), "\n")
cat("BPCells available:", requireNamespace("BPCells", quietly = TRUE), "\n")
cat("SeuratWrappers available:", requireNamespace("SeuratWrappers", quietly = TRUE), "\n")

# Clean environment
rm(check_sysreqs, load_pkgs, cran_packages, bioc_packages, gh_packages, pkg_spec, pkg_name)

# renv::snapshot(type = "all") # Uncomment to save package versions to renv.lock

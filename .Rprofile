local({
  # Named repos required by renv::snapshot() validation.
  # pak and BiocManager can leave unnamed entries — setting this here
  # ensures every session starts clean regardless of what packages do later.
  options(repos = c(
    CRAN     = "https://cloud.r-project.org",
    BioCsoft = "https://bioconductor.org/packages/3.22/bioc",
    BioCann  = "https://bioconductor.org/packages/3.22/data/annotation",
    BioCexp  = "https://bioconductor.org/packages/3.22/data/experiment"
  ))

  # Activate renv if the library has been restored (renv/ is gitignored;
  # run renv::restore() once after cloning to recreate it).
  if (file.exists("renv/activate.R")) source("renv/activate.R")
})

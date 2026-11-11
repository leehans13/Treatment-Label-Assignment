options(repos = c(CRAN = 'https://cloud.r-project.org'))
if (!requireNamespace('nebula', quietly = TRUE)) install.packages('nebula')
if (!requireNamespace('nebula', quietly = TRUE)) stop('Installation failed: nebula')

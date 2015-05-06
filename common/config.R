# Place all the common options here

func.filter.datadir <- function(dir.path) {
  files <- list.files(dir.path, full.names = TRUE)
  subset(files, grepl("^dataset_[.0-9]+", basename(files), perl = TRUE))
}

func.plink.filename <- function(plink.home) {
  sub(".bed", "", list.files(plink.files, pattern = "*.bed")[1])
}

# all input/output files root
root.dir <- "D:\\work\\bio\\workdir"

# common scripts root (normally it's a parent of this file)
lib.dir <- "D:\\work\\bio\\rlib\\common"

# libs aliases
lib.reader <- file.path(lib.dir, "readers.R")
lib.summary <- file.path(lib.dir, "summarizers.R")

# where all raw input files are stored
base <- file.path(root.dir, "raw")

# absolute paths to raw input 
raw.files <- func.filter.datadir(base)

# reader options
phe.as.fam <- FALSE
remove.dups <- TRUE
maf.thresh <- 0.01
ignore_cols <- NULL
use_cols <- NULL

# output folders homes 
home.lasso <- file.path(root.dir, "lasso_out")
home.gemma <- file.path(root.dir, "gemma_out")
home.moss <- file.path(root.dir, "moss_out")
home.corr <- file.path(root.dir, "correlation")
home.rf <- file.path(root.dir, "rf_out")
home.rbm <- file.path(root.dir, "rbm_out")
home.tmp <- "D:\\work\\bio\\tmp"

# set output precission
options("scipen" = 100, "digits" = 4)

# extra options for script (can be ignored)
extra <- ""

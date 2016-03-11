source("libs/functions.R")

echo("BacGWAS command line: ", commandArgs(trailingOnly = TRUE))
options("scipen" = 100, "digits" = 4)

.htmlOptions <- c("smartypants", "base64_images", "toc")
.startWd <- getwd()
.usageString <- "
Usage:
  bacgwas.R --plugin=NAME --input=DIR --output=DIR [--html]

Options:
  --plugin=NAME  Name of plugin to use. Corresponds to dir name within 
                 ./plugins folder.
  --input=DIR    Path to input directory. Input itself is plugin-specific, 
                 please read plugin docs for details.
  --output=DIR   Path to output directory. Actual output files are specific to
                 selected plugin type.
  --html         Switch plugin to html mode (check if plugin supports this 
                 mode before trying this option).
"

includer(c("knitr", "Cairo", "markdown", "docopt"))

source("conf/config.R")

.opt <- docopt(.usageString)
.plugin_home <- normalizePath(
  paste("plugins", .opt$plugin, sep = "/"),
  winslash = "/",
  mustWork = TRUE
)
.input <- normalizePath(.opt$input, winslash = "/", mustWork = TRUE)
.output <- normalizePath(.opt$output, winslash = "/", mustWork = FALSE)

if (!file.exists(.output))
  dir.create(.opt$output, recursive = TRUE)

tryCatch({
  setwd(.plugin_home)
  cparams <- c(config$common, config[[.opt$plugin]])
  
  if (.opt$html) {
    .markdownFile <- tempfile(
      pattern = "temp",
      tmpdir = tempdir(),
      fileext = ".md"
    )
    .markdownFile <- normalizePath(
      .markdownFile,
      winslash = "/",
      mustWork = FALSE
    )
    .opt$picsdir <- paste0(tempdir(), "/figure")
    .picsdir <- normalizePath(.opt$picsdir, winslash = "/", mustWork = FALSE)
    opts_chunk$set(
      dev = "png", 
      self.contained = TRUE, 
      dpi = 96,
      dev.args = list(type = "cairo"),
      fig.path = sub("([^/])$", "\\1/", .picsdir)
    )
    .report <- normalizePath(
      paste(.output, "result.html", sep = "/"), 
      winslash = "/", 
      mustWork = FALSE
    )
    .template <- normalizePath("init.rmd", winslash = "/", mustWork = TRUE)
    tryCatch({
      inject_args(cparams)
      knit(.template, .markdownFile, quiet = TRUE)
      markdownToHTML(
        .markdownFile, 
        output = .report,
        options = .htmlOptions,
        fragment.only = FALSE
      )
    }, finally = {
      if (file.exists(.markdownFile)) {
        echo("Remove intermediate markdown file: ", .markdownFile)
        file.remove(.markdownFile)
      }
      if (file.exists(.picsdir)) {
        echo("Remove pictures dir: ", .picsdir)
        unlink(.picsdir, recursive = TRUE)
      }
    })
  } else {
    source("init.R")
    if (exists("plugin_do")) {
      inject_args(cparams)
      plugin_do(.input, .output)
    } else {
      echo("Incorrect plugin: function plugin_do not defined!")
    }
  }
}, finally = {
  setwd(.startWd)
})

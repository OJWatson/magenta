## Build the shipped vignettes from their sources.
##
## magenta vignettes are expensive to run, so they are *pre-evaluated*. The
## editable sources live in `vignettes_src/` and contain live knitr chunks.
## This script renders each source, and writes a static `.Rmd` into
## `vignettes/` in which every chunk has already been replaced by its code and
## its output. `R CMD build` and pkgdown then only have to parse markdown, so
## building the package never re-runs a simulation.
##
## Usage, from the package root:
##
##     Rscript vignettes_src/build_vignettes.R              # build everything
##     Rscript vignettes_src/build_vignettes.R Resistance   # build selected ones
##
## Figures are written to `vignettes/<Name>_files/figure-html/` and referenced
## with a *relative* path, so that they resolve both in the installed vignette
## and on the pkgdown site.

build_vignette <- function(src, src_dir = "vignettes_src", out_dir = "vignettes") {

  stopifnot(file.exists(src))
  name <- sub("\\.Rmd$", "", basename(src))
  message("--- building ", name)

  rmarkdown::render(
    input = src,
    output_dir = out_dir,
    output_options = list("keep_md" = TRUE),
    quiet = TRUE
  )

  md <- file.path(out_dir, paste0(name, ".md"))
  rmd <- file.path(out_dir, paste0(name, ".Rmd"))
  if (!file.exists(md)) {
    stop("no markdown produced for ", name, " - is keep_md set?")
  }
  file.rename(md, rmd)

  # the rendered html is only an intermediate; the shipped vignette is the
  # static Rmd, which R CMD build and pkgdown render themselves
  unlink(file.path(out_dir, paste0(name, ".html")))

  lines <- readLines(rmd, warn = FALSE)

  # rmarkdown writes absolute figure paths when output_dir differs from the
  # source directory. Make them relative to the vignette instead.
  lines <- gsub(paste0(normalizePath(out_dir), "/"), "", lines, fixed = TRUE)
  lines <- gsub(paste0(file.path(getwd(), out_dir), "/"), "", lines, fixed = TRUE)

  # drop the intermediate frames of the progress bar, keeping the final one
  lines <- drop_progress_frames(lines)

  writeLines(lines, rmd)

  # pkgdown will only parse the static vignette once it has been knitted
  knitr::knit(rmd, output = tempfile(), quiet = TRUE)

  invisible(rmd)
}

## The simulation progress bar prints one line per refresh. Keep only the
## completed frame for each bar so the vignette is not full of noise.
drop_progress_frames <- function(lines) {
  is_frame <- grepl("^##\\s*Running:.*\\]\\s*\\d+%", lines)
  if (!any(is_frame)) {
    return(lines)
  }
  is_complete <- grepl("\\]\\s*100%", lines) & is_frame
  # within each run of consecutive frames, keep the last complete one only
  run_id <- cumsum(!is_frame)
  keep <- !is_frame
  for (r in unique(run_id[is_frame])) {
    idx <- which(is_frame & run_id == r)
    done <- idx[is_complete[idx]]
    keep[if (length(done)) max(done) else max(idx)] <- TRUE
  }
  lines[keep]
}

main <- function(selected = character(0)) {

  if (!file.exists("DESCRIPTION")) {
    stop("run this from the magenta package root")
  }

  srcs <- list.files("vignettes_src", pattern = "\\.Rmd$", full.names = TRUE)
  if (length(selected)) {
    keep <- sub("\\.Rmd$", "", basename(srcs)) %in% sub("\\.Rmd$", "", selected)
    if (!any(keep)) {
      stop("no vignette source matched: ", paste(selected, collapse = ", "))
    }
    srcs <- srcs[keep]
  }

  for (src in srcs) {
    build_vignette(src)
  }
  message("--- done (", length(srcs), " vignette(s))")
}

if (!interactive() && identical(environment(), globalenv())) {
  main(commandArgs(trailingOnly = TRUE))
}

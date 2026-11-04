# Builds the pkgdown site in docs/ (served by GitHub Pages) and makes sure
# MathJax is loaded from the public CDN (on some systems pkgdown links a local
# copy, which does not work online). Run from the package root:
#   source("dev/build_site.R")
options(pkgdown.internet = FALSE)
pkgdown::build_site(preview = FALSE, install = FALSE, new_process = FALSE)
for (f in list.files("docs", pattern = "\\.html$", recursive = TRUE, full.names = TRUE)) {
  x <- readLines(f, warn = FALSE)
  y <- gsub('<script src="file:///[^"]*MathJax.js"[^>]*>',
            '<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js" crossorigin="anonymous">', x)
  if (!identical(x, y)) writeLines(y, f)
}
file.create("docs/.nojekyll")

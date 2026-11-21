# Compare interrater against other implementations: the R tools.
#
# Reads data/<dataset>.csv (written by measure.py) and computes what irrCAC
# and irr compute on it. Writes r_tools.csv, one row per (dataset, weights,
# coefficient, quantity, tool), with the tool's version.
#
# What each tool covers:
#
# - irrCAC: Fleiss, AC1, Brennan-Prediger and Krippendorff on all the data,
#   missing ratings included, unweighted, linear and quadratic. Conger on
#   all the data. Cohen and Scott for two raters, and Light as the mean of
#   the pairwise Cohen's kappas, each on the subjects the raters in question
#   all rated. irrCAC reports estimates rounded to 5 decimals.
#   With weights, gwet.ac1.raw computes AC2, the weighted AC1, which is
#   the formula the package's AC2 follows.
# - irr: kripp.alpha, nominal and interval (alpha with quadratic weights on
#   1..K codes). kappam.fleiss (Fleiss, and Conger with exact = TRUE),
#   kappam.light and kappa2 (Cohen, unweighted, linear and quadratic), on
#   complete data only, since irr drops incomplete subjects. The null
#   standard errors of Fleiss and Cohen are value / statistic, the z
#   statistic irr reports.
#
# Fleiss (1971) has made-up rater slots, so Scott, Cohen, Conger and Light
# are skipped for it, as in measure.py. The weightings for a dataset are the
# ones interrater_values.csv has for it (none for Fleiss 1971, whose
# categories are unordered).
#
# Run after measure.py. The order is in that file's docstring.
#
# Dependencies: R, irrCAC, irr. Not the package.

suppressMessages({
  library(irrCAC)
  library(irr)
})

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
package <- read.csv(file.path(here, "interrater_values.csv"), stringsAsFactors = FALSE)
n_categories <- tapply(package$K, package$dataset, `[`, 1)

versions <- c(irrCAC = as.character(packageVersion("irrCAC")),
              irr = as.character(packageVersion("irr")))
rows <- list()
add <- function(dataset, weights, coef, quantity, tool, value) {
  rows[[length(rows) + 1]] <<- data.frame(
    dataset = dataset, weights = weights, coef = coef, quantity = quantity,
    tool = tool, tool_version = versions[[tool]], value = as.numeric(value)
  )
}
irrcac <- function(res) res$est$coeff.val[1]
irr_weight <- c(unweighted = "unweighted", linear = "equal", quadratic = "squared")

for (dataset in names(n_categories)) {
  K <- n_categories[[dataset]]
  d <- read.csv(file.path(here, "data", paste0(dataset, ".csv")), na.strings = "")
  complete <- d[complete.cases(d), ]
  is_complete <- nrow(complete) == nrow(d)
  identity <- dataset != "fleiss1971"
  labels <- seq_len(K)
  weightings <- unique(package$weights[package$dataset == dataset])

  # irrCAC
  for (w in weightings) {
    add(dataset, w, "fleiss", "estimate", "irrCAC", irrcac(fleiss.kappa.raw(d, weights = w, categ.labels = labels)))
    add(dataset, w, "ac1", "estimate", "irrCAC", irrcac(gwet.ac1.raw(d, weights = w, categ.labels = labels)))
    add(dataset, w, "bp", "estimate", "irrCAC", irrcac(bp.coeff.raw(d, weights = w, categ.labels = labels)))
    add(dataset, w, "kripp", "estimate", "irrCAC", irrcac(krippen.alpha.raw(d, weights = w, categ.labels = labels)))
    if (!identity) next
    add(dataset, w, "conger", "estimate", "irrCAC", irrcac(conger.kappa.raw(d, weights = w, categ.labels = labels)))
    if (ncol(d) == 2) {
      add(dataset, w, "scott", "estimate", "irrCAC", irrcac(fleiss.kappa.raw(d, weights = w, categ.labels = labels)))
      add(dataset, w, "cohen", "estimate", "irrCAC", irrcac(conger.kappa.raw(complete, weights = w, categ.labels = labels)))
    }
    pairs <- combn(ncol(d), 2)
    add(dataset, w, "light", "estimate", "irrCAC", mean(apply(pairs, 2, function(j) {
      pair <- d[, j]
      irrcac(conger.kappa.raw(pair[complete.cases(pair), ], weights = w, categ.labels = labels))
    })))
  }

  # irr
  add(dataset, "unweighted", "kripp", "estimate", "irr", kripp.alpha(t(as.matrix(d)), "nominal")$value)
  if ("quadratic" %in% weightings) add(dataset, "quadratic", "kripp", "estimate", "irr", kripp.alpha(t(as.matrix(d)), "interval")$value)
  if (is_complete) {
    f <- kappam.fleiss(complete)
    add(dataset, "unweighted", "fleiss", "estimate", "irr", f$value)
    add(dataset, "unweighted", "fleiss", "null_se", "irr", f$value / f$statistic)
    if (identity) {
      add(dataset, "unweighted", "conger", "estimate", "irr", kappam.fleiss(complete, exact = TRUE)$value)
      if (ncol(d) > 2) add(dataset, "unweighted", "light", "estimate", "irr", kappam.light(complete)$value)
    }
  }
  if (identity && ncol(d) == 2) {
    for (w in weightings) {
      k <- kappa2(complete, weight = irr_weight[[w]])
      add(dataset, w, "cohen", "estimate", "irr", k$value)
      add(dataset, w, "cohen", "null_se", "irr", k$value / k$statistic)
    }
  }
}

out <- do.call(rbind, rows)
write.csv(out, file.path(here, "r_tools.csv"), row.names = FALSE)
cat(sprintf("r_tools.csv: %d rows (R %s, irrCAC %s, irr %s)\n",
            nrow(out), getRversion(), versions[["irrCAC"]], versions[["irr"]]))

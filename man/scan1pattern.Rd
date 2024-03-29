% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggplot_scan1pattern.R, R/scan1pattern.R
\name{ggplot_scan1pattern}
\alias{ggplot_scan1pattern}
\alias{autoplot.scan1pattern}
\alias{scan1pattern}
\alias{summary.scan1pattern}
\title{Plot scan pattern usign ggplot2}
\usage{
ggplot_scan1pattern(
  object,
  map,
  plot_type = c("lod", "coef", "coef_and_lod"),
  patterns = object$patterns$founders,
  columns = 1:3,
  min_lod = 3,
  lodcolumn = seq_along(patterns),
  facet = "pheno",
  ...
)

\method{autoplot}{scan1pattern}(object, ...)

scan1pattern(
  probs1,
  phe,
  K = NULL,
  covar = NULL,
  map,
  patterns,
  condense_patterns = TRUE,
  blups = FALSE,
  do_scans = TRUE
)

\method{summary}{scan1pattern}(object, map, ...)
}
\arguments{
\item{object}{object of class \code{\link{scan1pattern}}}

\item{map}{genome map}

\item{plot_type}{type of plot from \code{c("lod","coef")}}

\item{patterns}{data frame of pattern information}

\item{columns}{columns used for coef plot}

\item{min_lod}{minimum LOD peak for contrast to be retained}

\item{lodcolumn}{columns used for scan1 plot (default all \code{patterns})}

\item{facet}{Plot facets if multiple phenotypes and patterns provided (default = \code{"pheno"}).}

\item{...}{additional parameters passed on to other methods}

\item{probs1}{object of class \code{\link[qtl2]{calc_genoprob}}}

\item{phe}{data frame with one phenotype}

\item{K}{kinship matrix}

\item{covar}{covariate matrix}

\item{condense_patterns}{remove snp_action from contrasts if TRUE}

\item{blups}{Create BLUPs if \code{TRUE}}

\item{do_scans}{Do scans if \code{TRUE}.}
}
\value{
object of class \code{\link[ggplot2]{ggplot}}

List containing:
\itemize{
\item{patterns} Data frame of summary for top patterns (column \code{founders} has pattern)
\item{dip_set} Diplotype sets for contrasts
\item{group} Group for each founder pattern
\item{scan} Object of class \code{\link[qtl2]{scan1}}.
\item{coef} Object of class \code{listof_scan1coef}. See package 'qtl2ggplot'.
}
}
\description{
Plot scan pattern usign ggplot2

Genome scan by pattern set
}
\examples{
dirpath <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex"

# Read DOex example cross from 'qtl2data'
DOex <- subset(qtl2::read_cross2(file.path(dirpath, "DOex.zip")), chr = "2")

\donttest{
# Download genotype probabilities
tmpfile <- tempfile()
download.file(file.path(dirpath, "DOex_genoprobs_2.rds"), tmpfile, quiet=TRUE)
pr <- readRDS(tmpfile)
unlink(tmpfile)

# Download SNP info for DOex from web and read as RDS.
tmpfile <- tempfile()
download.file(file.path(dirpath, "c2_snpinfo.rds"), tmpfile, quiet=TRUE)
snpinfo <- readRDS(tmpfile)
unlink(tmpfile)
snpinfo <- dplyr::rename(snpinfo, pos = pos_Mbp)

# Convert to SNP probabilities
snpinfo <- qtl2::index_snps(DOex$pmap, snpinfo)
snppr <- qtl2::genoprob_to_snpprob(pr, snpinfo)

# Scan SNPs
scan_snppr <- qtl2::scan1(snppr, DOex$pheno)
top_snps_tbl <- top_snps_pattern(scan_snppr, snpinfo)

# Summarize to find top patterns
patterns <- dplyr::arrange(summary(top_snps_tbl), dplyr::desc(max_lod))

# Scan using patterns.
scan_pat <- scan1pattern(pr, DOex$pheno, map = DOex$pmap, patterns = patterns)

# Summary of scan1pattern.
summary(scan_pat, DOex$pmap)
}

}
\author{
Brian S Yandell, \email{brian.yandell@wisc.edu}
}
\keyword{utilities}

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 7, fig.height = 5)

## -----------------------------------------------------------------------------
library(qtl2pattern)

## -----------------------------------------------------------------------------
dirpath <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex"

## -----------------------------------------------------------------------------
DOex <-
  subset(
    qtl2::read_cross2(file.path(dirpath, "DOex.zip")),
    chr = "2")

## -----------------------------------------------------------------------------
tmpfile <- tempfile()
download.file(file.path(dirpath, "DOex_genoprobs_2.rds"), tmpfile, quiet=TRUE)
pr <- readRDS(tmpfile)
unlink(tmpfile)

## ----eval = FALSE-------------------------------------------------------------
#  pr <- qtl2::calc_genoprob(DOex, error_prob=0.002)

## -----------------------------------------------------------------------------
tmpfile <- tempfile()
download.file(file.path(dirpath, "c2_snpinfo.rds"), tmpfile, quiet=TRUE)
snpinfo <- readRDS(tmpfile)
unlink(tmpfile)

## -----------------------------------------------------------------------------
snpinfo <- dplyr::rename(snpinfo, pos = pos_Mbp)
snpinfo <- qtl2::index_snps(DOex$pmap, snpinfo)

## -----------------------------------------------------------------------------
apr <- qtl2::genoprob_to_alleleprob(pr)
snpapr <- qtl2::genoprob_to_snpprob(apr, snpinfo)
dim(snpapr[["2"]])
dimnames(snpapr[["2"]])[[2]]

## -----------------------------------------------------------------------------
snppr <- qtl2::genoprob_to_snpprob(pr, snpinfo)
dim(snppr[["2"]])
dimnames(snppr[["2"]])[[2]]

## -----------------------------------------------------------------------------
rm(snpapr)

## -----------------------------------------------------------------------------
scan_snppr <- qtl2::scan1(snppr, DOex$pheno)

## -----------------------------------------------------------------------------
qtl2::find_peaks(scan_snppr, snpinfo)

## -----------------------------------------------------------------------------
top_snps_tbl <- top_snps_pattern(scan_snppr, snpinfo)

## -----------------------------------------------------------------------------
(patterns <- summary(top_snps_tbl))

## -----------------------------------------------------------------------------
head(summary(top_snps_tbl, "best"))

## -----------------------------------------------------------------------------
scan_pat <- scan1pattern(pr, DOex$pheno,
                         map = DOex$pmap,
                         patterns = patterns)

## -----------------------------------------------------------------------------
summary(scan_pat, DOex$pmap)

## -----------------------------------------------------------------------------
pat_names <- paste(c(52,43,16,20), colnames(scan_pat$scan), sep = "_")

## -----------------------------------------------------------------------------
plot(scan_pat$scan, DOex$pmap, lodcolumn = 2, col = "red", xlim = c(90,110), type = "b")
abline(v = c(96.5, 98.5), col = "darkgray", lwd = 2, lty = 2)
plot(scan_pat$scan, DOex$pmap, lodcolumn = 1, add = TRUE, col = "blue", type = "b")
plot(scan_pat$scan, DOex$pmap, lodcolumn = 3, add = TRUE, col = "purple", type = "b")
plot(scan_pat$scan, DOex$pmap, lodcolumn = 4, add = TRUE, col = "green", type = "b")
title("Scans for SDP 52 (blue), 43 (red), 16 (purple), 20 (green)")

## -----------------------------------------------------------------------------
plot(scan_pat$scan, DOex$pmap, lodcolumn = 2, col = "red")
abline(v = c(96.5, 98.5), col = "darkgray", lwd = 2, lty = 2)
plot(scan_pat$scan, DOex$pmap, lodcolumn = 1, add = TRUE, col = "blue")
plot(scan_pat$scan, DOex$pmap, lodcolumn = 3, add = TRUE, col = "purple")
plot(scan_pat$scan, DOex$pmap, lodcolumn = 4, add = TRUE, col = "green")
title("Scans for SDP 52 (blue), 43 (red), 16 (purple), 20 (green)")

## -----------------------------------------------------------------------------
oldpar <- par(mfrow = c(2,2))
cols <- c("blue","purple","red")
for(i in 1:4) {
  plot(scan_pat$coef[[i]], DOex$pmap, columns = 1:3, col = cols, xlim = c(90,110), type = "b")
  title(paste("Coefficients for", pat_names[i]))
  if(i == 3) {
    abline(v = c(96.5, 98.5), col = "darkgray", lwd = 2, lty = 2)
    legend(104, -3, legend = c("ref","het","alt"), 
       col = cols, lty = 1, lwd = 2)
  }
}
par(oldpar)

## -----------------------------------------------------------------------------
par(mfrow = c(2,2))
cols <- c("blue","purple","red")
for(i in 1:4) {
  plot(scan_pat$coef[[i]], DOex$pmap, columns = 1:3, col = cols)
  title(paste("Coefficients for", pat_names[i]))
  if(i == 3) {
    abline(v = c(96.5, 98.5), col = "darkgray", lwd = 2, lty = 2)
    legend(125, -3, legend = c("ref","het","alt"), 
       col = cols, lty = 1, lwd = 2)
  }
}

## -----------------------------------------------------------------------------
coefs2 <- qtl2::scan1coef(pr, DOex$pheno[,"OF_immobile_pct"], zerosum = FALSE)

## -----------------------------------------------------------------------------
x <- LETTERS[1:8]
ref <- outer(x,x,paste0)
ref <- ref[upper.tri(ref, diag=TRUE)]
ss <- 2 - sapply(stringr::str_split(ref, ""), function(x) x[1] == x[2])
alt <- ref[!stringr::str_detect(ref, c("A|B|C|D|G|H"))]
het <- ref[stringr::str_detect(ref, c("A|B|C|D|G|H")) & stringr::str_detect(ref, c("E|F"))]
altw <- ss[match(alt, ref, nomatch = 0)]
hetw <- ss[match(het, ref, nomatch = 0)]
refw <- ss[!stringr::str_detect(ref, c("E|F"))]
ref <- ref[!stringr::str_detect(ref, c("E|F"))]

## -----------------------------------------------------------------------------
coefsum <- coefs2
coefsum[,"AA"] <- apply(coefsum[,ref], 1, weighted.mean, w = refw)
coefsum[,"AB"] <- apply(coefsum[,het], 1, weighted.mean, w = hetw)
coefsum[,"BB"] <- apply(coefsum[,alt], 1, weighted.mean, w = altw)
colnames(coefsum)[1:3] <- c("ref","het","alt")

## -----------------------------------------------------------------------------
plot(coefsum, DOex$pmap, c("ref", "het", "alt", alt), xlim = c(90,110), ylim = c(-500, 500),
     col = c(1,8,7,2:4), type = "b")
abline(h = mean(DOex$pheno[,"OF_immobile_pct"]), lwd = 2, col = "darkgrey", lty = 2)
abline(v = c(96.5, 98.5), col = "darkgray", lwd = 2, lty = 2)
title("Allele Coefficients for OF_immobile_pct")
legend(105, 80, legend = c("ref", "het", "alt", alt), 
       col = c(1,8,7,2:4), lty = 1, lwd = 2)

## -----------------------------------------------------------------------------
plot(coefsum, DOex$pmap, c("ref", "het", alt[-1]), xlim = c(90,110), ylim = c(55,90),
     col = c(1,8,3:4), type = "b")
abline(h = mean(DOex$pheno[,"OF_immobile_pct"]), lwd = 2, col = "darkgrey", lty = 2)
abline(v = c(96.5, 98.5), col = "darkgray", lwd = 2, lty = 2)
title("Allele Pair Coefficients for OF_immobile_pct")
legend(107, 75, legend = c("ref", "het", alt[-1]), 
       col = c(1,8,3:4), lty = 1, lwd = 2)

## -----------------------------------------------------------------------------
scan_apr <- qtl2::scan1(apr, DOex$pheno)
coefs <- qtl2::scan1coef(apr, DOex$pheno)
coefs2 <- qtl2::scan1coef(pr, DOex$pheno)

## -----------------------------------------------------------------------------
alleles <- allele1(probD = pr,
                   scanH = scan_apr, coefH = coefs, coefD = coefs2, 
                   scan_pat = scan_pat, map = DOex$pmap, alt = "E")

## -----------------------------------------------------------------------------
ggplot2::autoplot(alleles, scan_apr, DOex$pmap, frame = FALSE)

## -----------------------------------------------------------------------------
aa <- subset(alleles, sources = levels(alleles$source)[c(1,4:6)])
ggplot2::autoplot(aa, scan_apr, DOex$pmap, frame = FALSE)

## -----------------------------------------------------------------------------
pos = patterns$max_pos[patterns$sdp == 48]
wh <- which.min(abs(pos - DOex$pmap[[1]]))
geno <- apply(snppr[[1]][,,wh], 1, function(x) which.max(x))
geno <- factor(c("ref","het","alt")[geno], c("ref","het","alt"))
dat <- data.frame(DOex$pheno, geno = geno)
dat$x <- jitter(rep(1, nrow(dat)))

## -----------------------------------------------------------------------------
ggplot2::ggplot(dat) +
  ggplot2::aes(x = x, y = OF_immobile_pct, col = geno, label = geno) +
  ggplot2::geom_boxplot() +
  ggplot2::geom_point() +
  ggplot2::facet_wrap(~ geno)

## -----------------------------------------------------------------------------
tmpfile <- tempfile()
download.file(file.path(dirpath, "c2_genes.rds"), tmpfile, quiet=TRUE)
gene_tbl <- readRDS(tmpfile)
unlink(tmpfile)
class(gene_tbl) <- c("feature_tbl", class(gene_tbl))

## -----------------------------------------------------------------------------
out <- merge_feature(top_snps_tbl, snpinfo, scan_snppr, exons = gene_tbl)
summary(out, "pattern")

## -----------------------------------------------------------------------------
out$snp_type <- sample(c(rep("intron", 40), rep("exon", nrow(out) - 40)))

## -----------------------------------------------------------------------------
ggplot2::autoplot(out, "OF_immobile_pct")

## -----------------------------------------------------------------------------
ggplot2::autoplot(out, "OF_immobile_pct", "consequence")

## ----eval = FALSE-------------------------------------------------------------
#  # Not quite right.
#  gene <- "Lrrc4c"
#  geneinfo <- dplyr::filter(gene_tbl, Name == gene)
#  parts <- geneinfo$start + diff(unlist(geneinfo[,c("start","stop")])) * (0:5) / 5
#  geneinfo <- geneinfo[rep(1,3),]
#  geneinfo$type <- "exon"
#  geneinfo$start <- parts[c(1,3,5)]
#  geneinfo$stop <- parts[c(2,4,6)]
#  gene_tbl <- rbind(gene_tbl[1,], geneinfo, gene_tbl[-1,])

## -----------------------------------------------------------------------------
qtl2::plot_genes(gene_tbl, xlim = c(96,99))

## -----------------------------------------------------------------------------
ggplot2::autoplot(gene_tbl)

## -----------------------------------------------------------------------------
ggplot2::autoplot(gene_tbl, top_snps_tbl = top_snps_tbl)

## -----------------------------------------------------------------------------
# Get Gene exon information.
out <- gene_exon(top_snps_tbl, gene_tbl)
summary(out, gene = out$gene[1])

## -----------------------------------------------------------------------------
ggplot2::autoplot(out, top_snps_tbl)

## -----------------------------------------------------------------------------
query_variants <- 
  qtl2::create_variant_query_func(
    file.path("qtl2shinyData", "CCmouse", "cc_variants.sqlite"))

## -----------------------------------------------------------------------------
query_genes <- 
  qtl2::create_gene_query_func(
    file.path("qtl2shinyData", "CCmouse", "mouse_genes.sqlite"))

## -----------------------------------------------------------------------------
objects(envir = environment(query_variants))

## -----------------------------------------------------------------------------
get("dbfile", envir = environment(query_variants))

## ----eval = FALSE-------------------------------------------------------------
#  saveRDS(query_variants, "query_variants.rds")
#  #
#  # other code
#  #
#  query_variants <- readRDS("query_variants.rds")

## ----eval = FALSE-------------------------------------------------------------
#  query_probs <-
#    create_probs_query_func(
#      file.path("qtl2shinyData", "CCmouse", "Recla"))


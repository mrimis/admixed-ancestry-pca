#!/usr/bin/env Rscript
# 04_ancestry.R
# Estimate global ancestry for every cohort sample in two independent ways,
# then write one covariate table.
#   Method 1: position in the reference PCA (fast, geometric).
#   Method 2: maximum likelihood with fixed reference allele frequencies.
#             Same likelihood as ADMIXTURE, with reference frequencies held fixed.
# Run from the repository root after 03_pca.sh:  Rscript scripts/04_ancestry.R

suppressPackageStartupMessages(library(data.table))
dir.create("results", showWarnings = FALSE)

refs <- c(AFR = "YRI", EUR = "CEU", EAS = "CHB")   # ancestry label = reference group

# ---------------------------------------------------------------------------
# Sample information
# ---------------------------------------------------------------------------
panel <- fread("data/affy6_samples.panel", header = TRUE, select = 1:2,
               col.names = c("IID", "POP"))
unrel <- fread("work/unrel.king.cutoff.in.id", header = TRUE)[[1]]
panel[, UNRELATED := as.integer(IID %in% unrel)]

read_sscore <- function(path) {
  x <- fread(path)
  setnames(x, 1, "IID")
  pcs <- grep("^PC[0-9]+_AVG$", names(x), value = TRUE)
  x <- x[, c("IID", pcs), with = FALSE]
  setnames(x, pcs, sub("_AVG$", "", pcs))
  x
}

# ---------------------------------------------------------------------------
# Method 1: ancestry from the reference PCA
# In PC space an admixed sample sits between its source groups. Its position
# between the three reference centres gives three weights that sum to 1.
# ---------------------------------------------------------------------------
ref_proj <- merge(read_sscore("work/ref_proj.sscore"), panel, by = "IID")

centres <- ref_proj[POP %in% refs & UNRELATED == 1,
                    .(PC1 = mean(PC1), PC2 = mean(PC2)), by = POP]
centres <- centres[match(refs, POP)]

A <- rbind(centres$PC1, centres$PC2, 1)             # 3 x 3 system
pc_weights <- function(pc1, pc2) t(solve(A, rbind(pc1, pc2, 1)))

cohort <- ref_proj[POP %in% c("ACB", "ASW")]
w <- pc_weights(cohort$PC1, cohort$PC2)
colnames(w) <- paste0(names(refs), "_PC")
pc_est <- cbind(cohort[, .(IID)], as.data.table(w))

# ---------------------------------------------------------------------------
# Method 2: maximum likelihood with fixed reference frequencies
# ---------------------------------------------------------------------------
read_freq <- function(pop) {
  f <- fread(sprintf("work/freq_%s.afreq", pop))
  setnames(f, 1, "CHROM")
  freq_col <- grep("^ALT.*FREQ", names(f), value = TRUE)[1]
  # Add half an allele to each count so no frequency is exactly 0 or 1.
  f[, .(ID, REF, ALT, P = (get(freq_col) * OBS_CT + 0.5) / (OBS_CT + 1))]
}
freqs <- lapply(refs, read_freq)
stopifnot(all(sapply(freqs, function(f) identical(f$ID, freqs[[1]]$ID))))
P <- sapply(freqs, function(f) f$P)                 # SNP x 3 matrix, ALT frequencies
snp <- freqs[[1]][, .(ID, REF, ALT)]

# Cohort genotypes. Column names look like rs123_A, where A is the counted allele.
raw <- fread("work/cohort_geno.raw")
geno_cols <- setdiff(names(raw), c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))
G <- as.matrix(raw[, geno_cols, with = FALSE])
rownames(G) <- raw$IID
counted <- sub("^.*_", "", geno_cols)
ids     <- sub("_[^_]*$", "", geno_cols)

# Line genotypes up with the frequency table, and make both count ALT.
idx <- match(snp$ID, ids)
stopifnot(!anyNA(idx))
G <- G[, idx]; counted <- counted[idx]
flip <- counted == snp$REF
stopifnot(all(counted == snp$ALT | flip))
G[, flip] <- 2 - G[, flip]

# The model, for one sample with ancestry proportions q (they sum to 1):
#   expected ALT frequency at each SNP   f = P %*% q
#   each genotype is Binomial(2, f)
# We look for the q that makes the observed genotypes most likely.
# th holds the first K-1 proportions. The last one is 1 - sum(th).
neg_loglik <- function(th, geno, P) {
  f <- drop(P %*% c(th, 1 - sum(th)))
  -sum(geno * log(f) + (2 - geno) * log(1 - f))
}
neg_gradient <- function(th, geno, P) {
  K <- ncol(P)
  f <- drop(P %*% c(th, 1 - sum(th)))
  w <- geno / f - (2 - geno) / (1 - f)
  -drop(crossprod(P[, -K, drop = FALSE] - P[, K], w))
}
estimate_q <- function(geno, P) {
  ok <- !is.na(geno)
  K <- ncol(P)
  fit <- constrOptim(
    theta = rep(1 / K, K - 1), f = neg_loglik, grad = neg_gradient,
    ui = rbind(diag(K - 1), -1), ci = c(rep(0, K - 1), -1),   # every q >= 0
    geno = geno[ok], P = P[ok, , drop = FALSE])
  c(fit$par, 1 - sum(fit$par))
}

ml <- t(apply(G, 1, estimate_q, P = P))
colnames(ml) <- names(refs)
ml_est <- data.table(IID = rownames(G), ml)

# ---------------------------------------------------------------------------
# One table: labels, ancestry from both methods, in-cohort PCs
# ---------------------------------------------------------------------------
cohort_pcs <- read_sscore("work/cohort_proj.sscore")

out <- Reduce(function(a, b) merge(a, b, by = "IID"),
              list(panel[POP %in% c("ACB", "ASW")], ml_est, pc_est, cohort_pcs))
setorder(out, POP, -AFR)

num <- setdiff(names(out), c("IID", "POP", "UNRELATED"))
out[, (num) := lapply(.SD, signif, 5), .SDcols = num]
fwrite(out, "results/cohort_ancestry_covariates.tsv", sep = "\t")

# ---------------------------------------------------------------------------
# Checks printed to the screen
# ---------------------------------------------------------------------------
cat("\nAgreement between the two methods (Pearson r):\n")
for (a in names(refs))
  cat(sprintf("  %s  r = %.4f\n", a, cor(out[[a]], out[[paste0(a, "_PC")]])))

cat("\nAfrican ancestry by population label:\n")
print(out[, .(n = .N, mean = round(mean(AFR), 3), min = round(min(AFR), 3),
              max = round(max(AFR), 3)), by = POP])

cat("\nSamples where the third component is above 5%:\n")
print(out[EAS > 0.05, .(IID, POP, AFR, EUR, EAS)])

# What does each in-cohort PC track? Compare its eigenvalue with the rest,
# and its correlation with African ancestry.
eig <- scan("work/cohort_pca.eigenval", quiet = TRUE)
cat("\nIn-cohort PCs:\n")
for (i in 1:5)
  cat(sprintf("  PC%d  eigenvalue %.2f   r with African ancestry %6.3f\n",
              i, eig[i], cor(out[[paste0("PC", i)]], out$AFR)))
cat("\nWrote results/cohort_ancestry_covariates.tsv\n")

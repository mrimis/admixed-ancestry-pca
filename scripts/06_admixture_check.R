#!/usr/bin/env Rscript
# 06_admixture_check.R  (optional, called by 06_admixture_check.sh)
# Compare ADMIXTURE's supervised estimates with the estimates from step 4.

suppressPackageStartupMessages(library(data.table))

refs <- c(AFR = "YRI", EUR = "CEU", EAS = "CHB")

fam <- fread("work/admix.fam", header = FALSE)
pop <- fread("work/admix.pop", header = FALSE)[[1]]
Q   <- as.matrix(fread("work/admix.3.Q", header = FALSE))
stopifnot(nrow(Q) == nrow(fam), length(pop) == nrow(fam))

# Which column of Q is which ancestry?
# In supervised mode every reference sample is fixed at 1 in its own column.
col_of <- sapply(refs, function(p) unname(which.max(colMeans(Q[pop == p, , drop = FALSE]))))
stopifnot(!anyDuplicated(col_of))

adm <- data.table(IID = fam[[2]],
                  AFR_ADMIXTURE = Q[, col_of[["AFR"]]],
                  EUR_ADMIXTURE = Q[, col_of[["EUR"]]],
                  EAS_ADMIXTURE = Q[, col_of[["EAS"]]])[pop == "-"]

cov <- fread("results/cohort_ancestry_covariates.tsv")
m <- merge(cov[, .(IID, POP, AFR, EUR, EAS)], adm, by = "IID")
setorder(m, POP, -AFR)
fwrite(m, "results/admixture_check.tsv", sep = "\t")

cat(sprintf("\nADMIXTURE (supervised, K = 3) against step 4, %d unrelated cohort samples:\n", nrow(m)))
for (a in names(refs)) {
  x <- m[[a]]; y <- m[[paste0(a, "_ADMIXTURE")]]
  cat(sprintf("  %s  r = %.4f   mean difference = %.4f   largest difference = %.4f\n",
              a, cor(x, y), mean(abs(x - y)), max(abs(x - y))))
}
cat("\nWrote results/admixture_check.tsv\n")

#!/usr/bin/env Rscript
# 05_plots.R
# Four figures that check the work.
# Run from the repository root after 04_ancestry.R:  Rscript scripts/05_plots.R

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
dir.create("figures", showWarnings = FALSE)

# Colour follows ancestry in every figure. Cohort samples are drawn in ink.
col_anc <- c(African = "#2a78d6", European = "#eb6834", `East Asian` = "#1baf7a")
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e6e5e0"; surface <- "#fcfcfb"

theme_tutorial <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = grid_col, linewidth = 0.3),
        plot.background  = element_rect(fill = surface, colour = NA),
        plot.title = element_text(face = "bold", size = 12, colour = ink),
        plot.subtitle = element_text(colour = ink2, size = 10),
        plot.title.position = "plot",
        axis.title = element_text(colour = ink2), axis.text = element_text(colour = ink2),
        strip.text = element_text(face = "bold", colour = ink, hjust = 0),
        legend.position = "top", legend.justification = "left",
        legend.title = element_blank(), legend.text = element_text(colour = ink2))
save_fig <- function(name, plot, w, h)
  ggsave(file.path("figures", name), plot, width = w, height = h, dpi = 160)

cov   <- fread("results/cohort_ancestry_covariates.tsv")
panel <- fread("data/affy6_samples.panel", header = TRUE, select = 1:2,
               col.names = c("IID", "POP"))
unrel <- fread("work/unrel.king.cutoff.in.id", header = TRUE)[[1]]

# ---------------------------------------------------------------------------
# Figure 1. Reference PCA with the cohort projected onto it
# ---------------------------------------------------------------------------
proj <- fread("work/ref_proj.sscore"); setnames(proj, 1, "IID")
proj <- merge(proj[, .(IID, PC1 = PC1_AVG, PC2 = PC2_AVG)], panel, by = "IID")
ref_name <- c(YRI = "African", CEU = "European", CHB = "East Asian")
refd <- proj[POP %in% names(ref_name) & IID %in% unrel]
refd[, ancestry := factor(ref_name[POP], levels = names(col_anc))]
cohd <- proj[POP %in% c("ACB", "ASW")]
labs_ref <- refd[, .(PC1 = mean(PC1), PC2 = mean(PC2)), by = .(POP, ancestry)]
labs_ref[, label := sprintf("%s reference (%s)", ancestry, POP)]
# Put each label beside its cluster, on the side with free space.
span <- diff(range(proj$PC1))
labs_ref[, hjust := ifelse(POP == "YRI", 1, 0)]
labs_ref[, PC1 := PC1 + ifelse(POP == "YRI", 0, 0.035) * span]
labs_ref[POP == "YRI", PC2 := PC2 - 0.07 * diff(range(proj$PC2))]

fig1 <- ggplot() +
  geom_point(data = refd, aes(PC1, PC2, fill = ancestry), shape = 21, size = 2.2,
             colour = surface, stroke = 0.3) +
  geom_point(data = cohd, aes(PC1, PC2, shape = POP), size = 1.9, colour = ink, alpha = 0.75) +
  geom_text(data = labs_ref, aes(PC1, PC2, label = label, hjust = hjust), colour = ink2, size = 3.3) +
  scale_fill_manual(values = col_anc, guide = "none") +
  scale_shape_manual(values = c(ACB = 16, ASW = 17),
                     labels = c(ACB = "Cohort: ACB", ASW = "Cohort: ASW")) +
  labs(title = "Most cohort samples sit on the African-European edge. Some do not.",
       subtitle = "PCs computed in unrelated reference samples. Cohort samples projected onto them.",
       x = "PC1", y = "PC2") +
  theme_tutorial
save_fig("fig1_reference_pca.png", fig1, 7.5, 5.6)

# ---------------------------------------------------------------------------
# Figure 2. Ancestry of every cohort sample, grouped by population label
# ---------------------------------------------------------------------------
bars <- melt(cov[, .(IID, POP, African = AFR, European = EUR, `East Asian` = EAS)],
             id.vars = c("IID", "POP"), variable.name = "ancestry", value.name = "share")
ord <- cov[order(POP, -AFR), IID]
bars[, IID := factor(IID, levels = ord)]
bars[, ancestry := factor(ancestry, levels = rev(names(col_anc)))]
n_lab <- cov[, .N, by = POP]
strip <- setNames(sprintf("%s  (n = %d)", n_lab$POP, n_lab$N), n_lab$POP)

fig2 <- ggplot(bars, aes(IID, share, fill = ancestry)) +
  geom_col(width = 1) +
  facet_grid(. ~ POP, scales = "free_x", space = "free_x", labeller = as_labeller(strip)) +
  scale_fill_manual(values = col_anc, breaks = names(col_anc),
                    labels = c("African", "European", "East Asian reference (third ancestry check)")) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%"), expand = c(0, 0)) +
  labs(title = "One label, many ancestries",
       subtitle = "Each bar is one person. Bars are sorted by African ancestry within each label.",
       x = NULL, y = "Share of ancestry") +
  theme_tutorial +
  theme(axis.text.x = element_blank(), panel.grid = element_blank(),
        panel.spacing = unit(10, "pt"))
save_fig("fig2_ancestry_by_label.png", fig2, 8.5, 4.2)

# ---------------------------------------------------------------------------
# Figure 3. What relatives do to PC1
# ---------------------------------------------------------------------------
wrong <- fread("work/wrong_relatives.eigenvec"); setnames(wrong, 1, "IID")
both <- rbind(
  merge(wrong[, .(IID, PC1)], cov[, .(IID, AFR)], by = "IID")[, panel := "wrong"],
  cov[, .(IID, PC1, AFR)][, panel := "right"])
r_lab <- both[, .(r = abs(cor(PC1, AFR))), by = panel]
titles <- c(wrong = "Relatives left in the PCA", right = "PCA on unrelated samples, relatives projected")
r_lab[, strip := sprintf("%s\n|r| = %.2f", titles[panel], r)]
both <- merge(both, r_lab[, .(panel, strip)], by = "panel")
both[, strip := factor(strip, levels = r_lab[match(c("wrong", "right"), panel), strip])]

fig3 <- ggplot(both, aes(AFR, PC1)) +
  geom_point(shape = 21, size = 2, fill = col_anc[["African"]], colour = surface, stroke = 0.3) +
  facet_wrap(~ strip, scales = "free_y") +
  scale_x_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(title = "Relatives break the link between PC1 and ancestry",
       subtitle = "Same 229 cohort samples in both panels.",
       x = "African ancestry (maximum likelihood)", y = "In-cohort PC1") +
  theme_tutorial
save_fig("fig3_relatives.png", fig3, 8.5, 4.2)

# ---------------------------------------------------------------------------
# Figure 4. Two independent ancestry estimates agree
# ---------------------------------------------------------------------------
agree <- rbind(cov[, .(ancestry = "African",    ml = AFR, pc = AFR_PC)],
               cov[, .(ancestry = "European",   ml = EUR, pc = EUR_PC)],
               cov[, .(ancestry = "East Asian", ml = EAS, pc = EAS_PC)])
agree[, ancestry := factor(ancestry, levels = names(col_anc))]
pct <- function(x) paste0(100 * x, "%")

fig4 <- ggplot(agree, aes(ml, pc, fill = ancestry)) +
  geom_abline(slope = 1, intercept = 0, colour = ink2, linewidth = 0.3) +
  geom_point(shape = 21, size = 2, colour = surface, stroke = 0.3) +
  facet_wrap(~ ancestry) +
  scale_fill_manual(values = col_anc, guide = "none") +
  scale_x_continuous(labels = pct, limits = c(-0.05, 1)) +
  scale_y_continuous(labels = pct, limits = c(-0.05, 1)) +
  coord_equal() +
  labs(title = "Two independent estimates give the same answer",
       subtitle = "Each point is one cohort sample. The line marks perfect agreement.",
       x = "Maximum likelihood estimate", y = "Estimate from PCA position") +
  theme_tutorial
save_fig("fig4_method_agreement.png", fig4, 8.5, 3.6)

cat("Wrote 4 figures to figures/\n")

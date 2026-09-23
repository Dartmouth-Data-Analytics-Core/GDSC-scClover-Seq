#!/usr/bin/env Rscript
# Figures for the rRNA summary (rule rrna_plots).
#
#   1. Composition: stacked bar per OB, every read of the OB as % of its total,
#      split into mapped rRNA / mapped other / unmapped rRNA /
#      unmapped <= short-max nt / unmapped > short-max nt with no alignment to the rRNA reference (remainder).
#   2. Remainder length profile: read length of the unmapped reads that did
#      not align to the rRNA reference and are longer than short-max, one line per OB,
#      as % of that OB's remainder reads.
#
# Usage:
#   plot_rrna_summary.R <rrna_summary.tsv> <out_composition.png> <out_length.png> \
#       <short_max> <well1.remainder_length_hist.tsv> [<well2...> ...]

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("usage: plot_rrna_summary.R summary.tsv out_comp.png out_len.png short_max hist1.tsv [hist2.tsv ...]")
summary_tsv <- args[1]
out_comp    <- args[2]
out_len     <- args[3]
short_max   <- as.integer(args[4])
hist_files  <- args[5:length(args)]

ob_levels <- function(x) {
  x <- unique(x)
  num <- suppressWarnings(as.integer(sub("^OB", "", x)))
  x[order(is.na(num), num, x)]
}

theme_rrna <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"))

# ---- 1. composition ---------------------------------------------------------
s <- read.delim(summary_tsv, check.names = FALSE, stringsAsFactors = FALSE) %>%
  filter(ob != "ALL")

le_col  <- paste0("unmapped_le", short_max, "nt")
gt_col  <- paste0("remainder_gt", short_max, "nt")
labels  <- c(
  "mapped: rRNA",
  "mapped: other",
  "unmapped: aligns to rRNA reference",
  paste0("unmapped: <= ", short_max, " nt (too short to align)"),
  paste0("unmapped: > ", short_max, " nt, no alignment to rRNA reference")
)

comp <- s %>%
  transmute(
    well, ob, total_reads,
    !!labels[1] := mapped_rRNA + mapped_MtrRNA,
    !!labels[2] := mapped_reads - mapped_rRNA - mapped_MtrRNA,
    !!labels[3] := unmapped_rRNA,
    !!labels[4] := .data[[le_col]],
    !!labels[5] := .data[[gt_col]]
  ) %>%
  pivot_longer(all_of(labels), names_to = "category", values_to = "reads") %>%
  mutate(pct      = 100 * reads / total_reads,
         category = factor(category, levels = rev(labels)),
         ob       = factor(ob, levels = ob_levels(ob)))

pal <- setNames(c("#7f7f7f", "#92c5de", "#2166ac", "#f4a582", "#b2182b"), levels(comp$category))

p1 <- ggplot(comp, aes(ob, pct, fill = category)) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
  facet_wrap(~ well, scales = "free_x") +
  scale_fill_manual(values = pal, name = NULL) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.01)) +
  labs(x = NULL, y = "% of the sample's reads",
       title = "Composition of each sample (OB)") +
  theme_rrna + theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2, reverse = TRUE))

ggsave(out_comp, p1, width = 9, height = 5.5, dpi = 200)

# ---- 2. remainder length profile -------------------------------------------
h <- bind_rows(lapply(hist_files, read.delim, stringsAsFactors = FALSE)) %>%
  filter(ob != "ALL", length > short_max) %>%
  group_by(well, ob) %>%
  mutate(pct = 100 * reads / sum(reads)) %>%
  ungroup() %>%
  mutate(ob = factor(ob, levels = ob_levels(ob)))

p2 <- ggplot(h, aes(length, pct, colour = ob)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ well) +
  labs(x = "read length (nt)", y = "% of the OB's remainder reads",
       colour = NULL,
       title = paste0("Length of unmapped reads > ", short_max, " nt with no alignment to the rRNA reference")) +
  theme_rrna

ggsave(out_len, p2, width = 9, height = 4.5, dpi = 200)

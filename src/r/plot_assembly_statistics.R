#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

assembly_data_raw <- fread("output/assembly_statistics/statistics.txt")

# parse algorithm
assembly_data_raw[grep("assembly.scafSeq", basename(filename)),
                  algorithm := "SOAPdenovo2"]
assembly_data_raw[grep("final.scaffolds.fa", basename(filename)),
                  algorithm := "meraculous"]

# parse data source
assembly_data_raw[grep("bin_reads_by_coverage", filename),
                  data_type := '"PCR-free binned"']
assembly_data_raw[grep("bbduk", filename),
                  data_type := '"PCR-free"']

# parse kmer
assembly_data_raw[, k := as.numeric(gsub(".*run_([[:digit:]]{2})mer.*",
                                         "\\1",
                                         filename))]

# modify units
assembly_data_raw[, length_mb := contig_bp / 1e6]
assembly_data_raw[, coverage := (100 - gap_pct)]
assembly_data_raw[, l50_kb := scaf_L50 / 1e3]
assembly_data_raw[, scaffolds_thousands := n_scaffolds / 1e3]

# column name/order
mv <- c(
    "length_mb" = '"Contig length (Mbp)"',
    "scaffolds_thousands" = '"Scaffolds (K)"',
    "l50_kb" = '"Scaffold"~italic(N)[50]~"(Kbp)"',
    "coverage" = '"Coverage (%)"')


pd <- melt(assembly_data_raw,
           id.vars = c("algorithm", "data_type", "k"),
           measure.vars = names(mv))

# set variable order
pd[, variable_ordered := factor(
    plyr::revalue(variable, replace = mv),
    levels = mv)]

# convert k to factor
pd[, k := factor(k, levels = sort(unique(k)))]

hs <- RColorBrewer::brewer.pal(6, "YlOrRd")
ggplot(pd, aes(x = algorithm, y = value, group = k, fill = k)) +
    theme(
        axis.text.x = element_text(size = 6,
                                   angle = 40,
                                   face = "italic",
                                   hjust = 1,
                                   vjust = 1),
        strip.text = element_text(size = 12),
        panel.spacing = unit(2, "pt")) +
    scale_fill_manual(values = hs[c(2, 4, 6)]) +
    facet_grid(variable_ordered ~ data_type,
               scales = "free_y",
               labeller = label_parsed) +
    xlab(NULL) + ylab(NULL) +
    geom_col(width = 0.6, position = position_dodge(width = 0.8))


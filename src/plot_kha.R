#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(bit64)
library(ggplot2)
library(scales)


###########
# GLOBALS #
###########

hist_before_file <- snakemake@input[["hist"]]
hist_after_file <- snakemake@input[["hist_out"]]
plot_file <- snakemake@output[["plot"]]

# find files (dev)
# hist_files <- list.files("output/020_norm",
#                          pattern = "Male_Bee_1_hist.*txt",
#                          full.names = TRUE)
# 
# names(hist_files) <- ifelse(grepl("-out", hist_files), "Normalised", "Raw")

# dev
# hist_before_file <- "output/020_norm/asw_hist.txt"
# hist_after_file <- "output/020_norm/asw_hist-out.txt"
# plot_file <- "test/asw_kha.pdf"

########
# MAIN #
########

hist_files <- c(Raw = hist_before_file, Normalised = hist_after_file)
hist_data_list <- lapply(hist_files, fread)

CalculateKhaStats <- function(x) {
    setkey(x, `#Depth`)
    x <- x[!last(x)]
    x[, cum_sum := cumsum(as.numeric(Raw_Count))]
    x[, percent_kmers := 100 * cum_sum / sum(Raw_Count)]
    return(x)
}

kha_data_list <- lapply(hist_data_list, CalculateKhaStats)
kha <- rbindlist(kha_data_list, idcol = "type")

gp <- ggplot(kha, aes(x = `#Depth`, y = percent_kmers, colour = type)) +
    theme_minimal(base_size = 8) +
    theme(legend.justification = c("right", "bottom"),
          legend.position = c(0.9, 0.1),
          legend.key.size = unit(0.5, "lines"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_color_viridis_d(guide = guide_legend(title = NULL)) +
    xlab("Depth") + ylab("Cumulative percentage of read 31-mers") +
    ylim(c(0, 100)) +
    scale_x_continuous(trans = log_trans(base = 4),
                       breaks = trans_breaks(function(x) log(x, 4),
                                             function(x) 4^x)) +
    geom_path(alpha = 0.75)


# find the elbows
kha_diff <- kha[`#Depth` <= 256,
                .(diff_pct = diff(diff(Raw_Count) > 0) != 0,
                  `#Depth` = `#Depth`[c(-1, -2)],
                  percent_kmers = percent_kmers[c(-1, -2)]),
                by = type]

kha_diff[, percent_repeat :=
             paste0("~", round(100 - percent_kmers, 0), "% repeats")]

gp_with_elbows <- gp +
    geom_hline(data = kha_diff[diff_pct == TRUE][c(5, 21)],
               mapping = aes(yintercept = percent_kmers,
                             colour = type),
               linetype = 2,
               show.legend = FALSE) +
    geom_text(data = kha_diff[diff_pct == TRUE][c(5, 21)],
               mapping = aes(label = percent_repeat,
                             x = 0),
              hjust = "left",
              vjust = -0.1,
              size = 2,
              show.legend = FALSE)

wo <- grid::convertUnit(grid::unit(483/3, "pt"), "mm", valueOnly = TRUE)
ho <- grid::convertUnit(grid::unit(664/3, "pt"), "mm", valueOnly = TRUE)
ggsave(filename = plot_file,
       plot = gp_with_elbows,
       device = cairo_pdf,
       width = wo,
       height = ho,
       units = "mm")

sessionInfo()

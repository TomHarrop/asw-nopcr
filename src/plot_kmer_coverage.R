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
peak_file <- snakemake@input[["peaks"]]
plot_file <- snakemake@output[["plot"]]

# dev
# hist_before_file <- "output/020_norm/asw_hist.txt"
# hist_after_file <- "output/020_norm/asw_hist-out.txt"
# peak_file <- "output/020_norm/asw_peaks.txt"
# plot_file <- "test/asw_kha.pdf"

########
# MAIN #
########

# read data
peaks <- fread(paste("grep '^[^#]'", peak_file))

hist_files <- c(Raw = hist_before_file, Normalised = hist_after_file)
hist_data_list <- lapply(hist_files, fread)
combined_data <- rbindlist(hist_data_list, idcol = "type")

# arrange plot
combined_data[, type := factor(type, levels = c("Raw", "Normalised"))]

# hlines
mincov <- peaks[2, V1]
p1 <- peaks[2, V2]
maxcov <- peaks[2, V3]
peak_pd <- combined_data[type == "Raw" & between(`#Depth`, mincov, maxcov)]

# plot title
# gt <- paste0(
#     p1, "× 31-mer coverage. ",
#     "Main peak: ", mincov, "×–", maxcov, "×"
# )

# plot
vd <- viridisLite::viridis(3)
line_col <- vd[[1]]
peak_col <- alpha(vd[[2]], 0.5)
kmer_plot <- ggplot(combined_data,
                    aes(x = `#Depth`,
                        y = Unique_Kmers,
                        linetype = type)) +
    theme_minimal(base_size = 8) +
    theme(legend.position = c(5/6, 2/4)) +
    geom_path(alpha = 0.75, colour = line_col) +
    geom_ribbon(data = peak_pd,
                mapping = aes(ymin = 0,
                              ymax = Unique_Kmers,
                              x = `#Depth`), colour = NA, fill = peak_col,
                show.legend = FALSE) +
    scale_linetype_discrete(guide = guide_legend(title = NULL)) +
    scale_y_continuous(
        trans = "log10",
        labels = trans_format("log10", math_format(10^.x)),
        breaks = trans_breaks("log10", function(x) 10^x)) +
    scale_x_continuous(trans = log_trans(base = 4),
                       breaks = trans_breaks(function(x) log(x, 4),
                                             function(x) 4^x)) +
    xlab("31-mer depth") + ylab("Number of unique 31-mers")

# write output
wo <- grid::convertUnit(grid::unit(483, "pt"), "mm", valueOnly = TRUE)
ho <- grid::convertUnit(grid::unit(664/4, "pt"), "mm", valueOnly = TRUE)
ggsave(filename = plot_file,
       plot = kmer_plot,
       device = cairo_pdf,
       width = wo,
       height = ho,
       units = "mm")

# write session info
sessionInfo()

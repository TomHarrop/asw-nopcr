#!/usr/bin/env Rscript

rutils::GenerateMessage("Print kmer distribution plots")

library(data.table)
library(bit64)
library(ggplot2)
library(scales)

# IO files
command.args <- commandArgs(trailingOnly = TRUE)
# command.args <- c(
#     "--fq", "output/bbnorm/ASW_normalised.fastq.gz",
#     "-r", "output/bbnorm/kmer_distribution_plot.pdf")
parsed.args <- argparsR::ParseArguments(
    accepted.switches = list(
        "input_fq" = "--fq",
        "output_pdf" = "-r"),
    command.args)

# data
input_directory <- dirname(parsed.args$input_fq)
output_directory <- dirname(parsed.args$output_pdf)
before_data_file <- list.files(input_directory,
                              pattern = "before",
                              full.names = TRUE)
after_data_file <- list.files(input_directory,
                             pattern = "after",
                             full.names = TRUE)
peaks_data_file <- list.files(input_directory,
           pattern = "peaks",
           full.names = TRUE)

before_data <- fread(before_data_file)
after_data <- fread(after_data_file)
peaks <- fread(peaks_data_file)

# combine
before_data[, type := "Raw"]
after_data[, type := "Normalised"]

combined_data <- rbind(before_data, after_data)

# arrange plot
combined_data[, type := factor(type, levels = c("Raw", "Normalised"))]

# hlines
mincov <- peaks[1, V1]
p1 <- peaks[1, V2]
p2 <- peaks[2, V2]
maxcov <- peaks[2, V3]

# plot
Set1 <- RColorBrewer::brewer.pal(9, "Set1")
kmer_plot <- ggplot(combined_data, aes(x = `#Depth`, y = Unique_Kmers,
                          colour = type)) +
    theme_minimal() +
    theme(legend.position = c(5/6, 2/4)) +
    geom_vline(xintercept = c(mincov, p1, p2, maxcov),
               colour = Set1[9]) +
    geom_path(alpha = 0.75) +
    scale_colour_brewer(palette = "Set1",
                        guide = guide_legend(title = NULL)) +
    scale_y_continuous(
        trans = "log10",
        labels = trans_format("log10", math_format(10^.x)),
        breaks = trans_breaks("log10", function(x) 10^x)) +
    scale_x_continuous(trans = log_trans(base = 4),
                       breaks = trans_breaks(function(x) log(x, 4),
                                             function(x) 4^x)) +
    xlab("31-mer depth") + ylab("Number of unique 31-mers")

# save output
ggsave(filename = parsed.args$output_pdf,
       plot = kmer_plot,
       width = 10,
       height = 7.5)

# logs
log.file <- paste0(output_directory, "/kmer_plot.SessionInfo.txt")
s.inf <- rutils::GitSessionInfo()
writeLines(s.inf, log.file)

quit("no", 0)

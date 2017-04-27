#!/usr/bin/env Rscript

rutils::GenerateMessage("Print base quality plot")

library(data.table)
library(ggplot2)

# IO files
command.args <- commandArgs(trailingOnly = TRUE)
# command.args <- c(
#  "--fq", "output/bbduk/ASW_filtered_trimmed.fastq.gz",
#  "-r", "output/bbduk/quality_histogram_plot.pdf")
parsed.args <- argparsR::ParseArguments(
    accepted.switches = list(
        "input_fq" = "--fq",
        "output_pdf" = "-r"),
    command.args)
input_directory <- dirname(parsed.args$input_fq)
output_directory <- dirname(parsed.args$output_pdf)

# data
qhist_data_file <- list.files(input_directory,
                              pattern = "^qhist.txt$",
                              full.names = TRUE)
qhist <- fread(qhist_data_file)
pd <- melt(qhist, id.vars = "#BaseNum")

# plot
Set1 <- RColorBrewer::brewer.pal(9, "Set1")[-6]
g <- ggplot(pd, aes(x = `#BaseNum`, y = value, colour = variable)) +
    theme_minimal() +
    scale_colour_manual(values = Set1,
                        guide = guide_legend(title = NULL)) +
    xlab("Base number") + ylab("Quality") +
    geom_path()

# save output
ggsave(filename = parsed.args$output_pdf,
       plot = g,
       width = 10,
       height = 7.5)

# logs
log.file <- paste0(output_directory, "/quality_histogram.SessionInfo.txt")
s.inf <- rutils::GitSessionInfo()
writeLines(s.inf, log.file)

quit("no", 0)

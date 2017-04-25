#!/usr/bin/env Rscript

rutils::GenerateMessage("Print read uniqueness histogram")

library(data.table)
library(ggplot2)

# IO files
command.args <- commandArgs(trailingOnly = TRUE)
# command.args <- c(
#     "-y", "output/bbduk/uniqueness_histogram.txt",
#     "-r", "output/bbduk/uniqueness_histogram.pdf")
parsed.args <- argparsR::ParseArguments(
    accepted.switches = list(
        "other_input" = "-y",
        "output_pdf" = "-r"),
    command.args)

# data
output_directory <- dirname(parsed.args$output_pdf)
pd_wide <- fread(parsed.args$other_input)
data_cols <- c("first", "r1_first", "r2_first", "pair", "avg_quality",
               "perfect_prob")
pd <- melt(pd_wide, id.vars = "#count", measure.vars = data_cols)

# plot
Set1 <- RColorBrewer::brewer.pal(9, "Set1")[-6]
g <- ggplot(pd, aes(x = `#count`/1e6,
               y = value,
               colour = variable)) +
    theme_minimal() +
    scale_colour_manual(values = Set1,
                        guide = guide_legend(title = NULL)) +
    xlab("Millions of reads processed") + ylab("Percentage") +
    geom_path()

# save output
ggsave(filename = parsed.args$output_pdf,
       plot = g,
       width = 10,
       height = 7.5)

# logs
log.file <- paste0(output_directory, "/uniqueness_histogram.SessionInfo.txt")
s.inf <- rutils::GitSessionInfo()
writeLines(s.inf, log.file)

quit("no", 0)

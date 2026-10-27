library(ggplot2)
library(dplyr)

read_beast_log <- function(path) {
  lines <- readLines(path)
  header_line <- max(grep("^\\s*Sample\\s+posterior", lines))
  # Read with row.names=NULL so the extra trailing "--" column doesn't confuse R
  df <- read.table(path, header = TRUE, skip = header_line - 1,
                   comment.char = "", fill = TRUE, check.names = FALSE,
                   row.names = NULL)
  # After row.names=NULL, columns shift: first col is the old row.names placeholder
  # rename: drop first unnamed col if it appeared
  if (colnames(df)[1] == "row.names" || !grepl("Sample", colnames(df)[1])) {
    colnames(df)[1] <- "Sample"
  }
  df <- df[!is.na(suppressWarnings(as.numeric(as.character(df$Sample)))), ]
  df$Sample <- as.numeric(as.character(df$Sample))
  df$posterior <- as.numeric(as.character(df$posterior))
  df
}

log_files <- c(
  "1" = "/Users/Flo/Downloads/Boxfish/analyses/beast/run4/run4_beast_1.log",
  "2" = "/Users/Flo/Downloads/Boxfish/analyses/beast/run4/run4_beast_2.log",
  "3" = "/Users/Flo/Downloads/Boxfish/analyses/beast/run4/run4_beast_3.log",
  "4" = "/Users/Flo/Downloads/Boxfish/analyses/beast/run4/run4_beast_4.log"
)

chains <- lapply(names(log_files), function(chain_id) {
  cat("Reading chain", chain_id, "...\n")
  df <- read_beast_log(log_files[[chain_id]])
  burnin_cut <- quantile(df$Sample, 0.25)
  df <- df[df$Sample > burnin_cut, ]
  df$chain <- paste("Chain", chain_id)
  df[, c("Sample", "posterior", "chain")]
})

all_chains <- bind_rows(chains)
all_chains$chain <- factor(all_chains$chain, levels = paste("Chain", 1:4))

p <- ggplot(all_chains, aes(x = Sample, y = posterior)) +
  geom_line(colour = "#2c7bb6", linewidth = 0.3, alpha = 0.8) +
  facet_wrap(~ chain, ncol = 2, scales = "free_x") +
  labs(
    x = "Generation",
    y = "Posterior probability (ln)",
    title = "BEAST2 Run 4 — posterior traces (25% burnin removed)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey90", colour = "grey60"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 11, hjust = 0.5)
  )

out <- "/Users/Flo/Downloads/Boxfish/results/run4-beast-tracer.pdf"
ggsave(out, plot = p, width = 8, height = 6, device = "pdf")
cat("Saved to", out, "\n")

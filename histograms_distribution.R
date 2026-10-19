library(readxl)
library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript plot_distributions.R <input.xlsx> <output_ribosome.png> <output_orthogroup.png>")
}

input_file <- args[1]
output_ribo <- args[2]
output_ortho <- args[3]

df_raw <- read_excel(input_file, col_names = FALSE)

taxons <- as.character(df_raw[1, 3:ncol(df_raw)])
counts_ribo <- as.numeric(df_raw[2, 3:ncol(df_raw)])

df_ribo <- data.frame(Taxon = taxons, Count = counts_ribo) %>%
  filter(!is.na(Taxon)) %>%
  group_by(Taxon) %>%
  summarise(MaxCount = max(Count, na.rm = TRUE)) %>%
  arrange(MaxCount) %>%
  mutate(Taxon = factor(Taxon, levels = Taxon))

percentile_probs <- c(0.01, 0.25, 0.5, 0.75, 1)
quantiles <- quantile(df_ribo$MaxCount, probs = percentile_probs)

get_closest_taxon <- function(q_val, df) {
  closest_row <- df[which.min(abs(df$MaxCount - q_val)), ]
  return(as.character(closest_row$Taxon))
}

quantile_labels <- unique(sapply(quantiles, get_closest_taxon, df = df_ribo))
df_percentiles <- df_ribo %>% filter(Taxon %in% quantile_labels)

plot_ribo <- ggplot(df_ribo, aes(x = Taxon, y = MaxCount)) +
  geom_area(aes(group = 1), fill = "#69b3a2", alpha = 0.6) +
  geom_line(aes(group = 1), color = "#4a8a7a", linewidth = 1) +
  geom_segment(
    data = df_percentiles,
    aes(x = Taxon, xend = Taxon, y = 0, yend = MaxCount),
    linetype = "dashed",
    color = "#e07a5f",
    linewidth = 0.8
  ) +
  scale_x_discrete(breaks = quantile_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", color = "#333333"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Distribution of the number of proteins by ribosome type",
    x = "Taxon and ribosome type",
    y = "Number of proteins in a ribosome"
  )

ortho_counts <- as.numeric(df_raw[[2]][4:nrow(df_raw)])
df_ortho <- data.frame(Count = ortho_counts) %>%
  filter(!is.na(Count))

plot_ortho <- ggplot(df_ortho, aes(x = Count)) +
  geom_histogram(
    fill = "#8d99ae",
    color = "#2b2d42",
    bins = 30,
    alpha = 0.8
  ) +
  scale_x_continuous(n.breaks = 12) +
  scale_y_continuous(n.breaks = 8) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#edf2f4")
  ) +
  labs(
    title = "Distribution of the number of proteins in an orthogroup",
    x = "Number of proteins in an orthogroup",
    y = "Number of orthogroups"
  )

ggsave(output_ribo, plot = plot_ribo, width = 10, height = 6, dpi = 300)
ggsave(output_ortho, plot = plot_ortho, width = 8, height = 5, dpi = 300)
cat("[+] Plots saved to:", output_ribo, "and", output_ortho, "\n")

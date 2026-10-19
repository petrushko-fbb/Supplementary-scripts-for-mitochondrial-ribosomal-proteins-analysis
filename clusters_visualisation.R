library(readxl)
library(writexl)
library(pheatmap)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript cluster_and_plot.R <input.xlsx> <output_clustered.xlsx> <output_binary.pdf> <output_colored.png>")
}

input_file <- args[1]
output_xlsx <- args[2]
output_pdf <- args[3]
output_png <- args[4]

df <- read_excel(input_file, col_names = FALSE)

meta_rows <- 1:3
meta_cols <- 1:2

ortho_names <- as.character(df[-meta_rows, 1, drop = TRUE])
org_names <- as.character(df[3, -meta_cols, drop = TRUE])
raw_taxons <- as.character(df[1, -meta_cols, drop = TRUE])
taxons <- trimws(gsub(",.*", "", raw_taxons))

raw_matrix <- df[-meta_rows, -meta_cols]

bin_matrix <- matrix(0, nrow = nrow(raw_matrix), ncol = ncol(raw_matrix))
bin_matrix[!is.na(raw_matrix) & trimws(raw_matrix) != ""] <- 1
rownames(bin_matrix) <- ortho_names
colnames(bin_matrix) <- org_names

dist_rows <- dist(bin_matrix, method = "binary")
dist_cols <- dist(t(bin_matrix), method = "binary")

hc_rows <- hclust(dist_rows, method = "average")
hc_cols <- hclust(dist_cols, method = "average")

row_order <- hc_rows$order
col_order <- hc_cols$order

new_col_indices <- c(1, 2, col_order + 2)
new_row_indices <- c(1, 2, 3, row_order + 3)
new_df <- df[new_row_indices, new_col_indices]

write_xlsx(new_df, output_xlsx, col_names = FALSE)
cat("[+] Clustered table saved to:", output_xlsx, "\n")

get_compartment <- function(name) {
  if (grepl("_mito$", name, ignore.case = TRUE)) return("Mitochondria")
  if (grepl("_cyto$", name, ignore.case = TRUE)) return("Cytosol")
  if (grepl("_chloro$", name, ignore.case = TRUE)) return("Chloroplast")
  return("Other")
}

get_base_org <- function(name) {
  gsub("_mito$|_cyto$|_chloro$", "", name, ignore.case = TRUE)
}

compartments <- sapply(org_names, get_compartment)
base_orgs <- sapply(org_names, get_base_org)

comp_colors <- c(
  Mitochondria = "#984ea3",
  Chloroplast = "#4daf4a",
  Cytosol = "#ff7f00",
  Other = "#999999"
)

annot_df_org <- data.frame(
  Compartment = compartments,
  Organism = base_orgs,
  row.names = org_names
)

unique_orgs <- unique(base_orgs)
set.seed(42)
org_colors <- sample(rainbow(length(unique_orgs), s = 0.7, v = 0.8))
names(org_colors) <- unique_orgs

ann_colors_org <- list(
  Compartment = comp_colors,
  Organism = org_colors
)

pdf(output_pdf, width = 25, height = 15)
pheatmap(
  bin_matrix,
  color = c("#f0f0f0", "#2ca25f"),
  cluster_rows = hc_rows,
  cluster_cols = hc_cols,
  annotation_col = annot_df_org,
  annotation_colors = ann_colors_org,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_col = 8,
  legend = FALSE
)
dev.off()
cat("[+] Binary heatmap saved to:", output_pdf, "\n")

is_mito <- compartments == "Mitochondria"
plot_matrix <- bin_matrix

for (i in seq_len(nrow(bin_matrix))) {
  row_data <- bin_matrix[i, ]
  if (sum(row_data) == 0) next

  present_outside_mito <- any(row_data[!is_mito] == 1)

  if (!present_outside_mito) {
    taxons_present <- unique(taxons[row_data == 1])

    if (length(taxons_present) > 1) {
      plot_matrix[i, row_data == 1] <- 2
    } else if (length(taxons_present) == 1) {
      t_name <- taxons_present[1]
      if (t_name == "Kinetoplastida") {
        plot_matrix[i, row_data == 1] <- 3
      } else if (t_name == "Mammalia") {
        plot_matrix[i, row_data == 1] <- 4
      } else if (t_name == "Chlorophyta") {
        plot_matrix[i, row_data == 1] <- 5
      } else if (t_name == "Fungi") {
        plot_matrix[i, row_data == 1] <- 6
      }
    }
  }
}

final_palette <- c("#dcdcdc", "#696969", "#8b4513", "#ff7f00", "#377eb8", "#4daf4a", "#984ea3")
my_breaks <- seq(-0.5, 6.5, length.out = 8)

annot_df_taxon <- data.frame(
  Compartment = compartments,
  Taxon = taxons,
  row.names = org_names
)

distinct_colors <- c(
  "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
  "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe",
  "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000",
  "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080"
)

unique_taxons <- unique(taxons)
taxon_colors <- rep(distinct_colors, length.out = length(unique_taxons))
names(taxon_colors) <- unique_taxons

ann_colors_taxon <- list(
  Compartment = comp_colors[names(comp_colors) != "Other"],
  Taxon = taxon_colors
)

png(output_png, width = 34, height = 20, units = "in", res = 600)
pheatmap(
  plot_matrix,
  color = final_palette,
  breaks = my_breaks,
  cluster_rows = hc_rows,
  cluster_cols = hc_cols,
  annotation_col = annot_df_taxon,
  annotation_colors = ann_colors_taxon,
  show_rownames = TRUE,
  fontsize_row = 5,
  show_colnames = TRUE,
  fontsize_col = 14,
  treeheight_row = 150,
  treeheight_col = 100,
  legend = FALSE
)
dev.off()
cat("[+] Taxon-colored heatmap saved to:", output_png, "\n")

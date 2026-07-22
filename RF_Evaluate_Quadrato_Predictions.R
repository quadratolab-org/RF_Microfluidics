# Standalone evaluation script for RF_Classifier_Fetal_onto_Quadrato_2nd_Tri_Only_No_MT_RP
# Run this on already-processed output files — no model retraining needed.
#
# Produces:
#   1. Row-normalized confusion matrix heatmap + CSV
#      (rows = organoid cell type, columns = predicted fetal class, values = fraction)
#   2. Count table: how many organoid cells were assigned each predicted fetal class
#   3. Count table: how many organoid cells belong to each original organoid label
#   4. Max-prediction-probability histogram per organoid cell type
#      (shows whether predictions are confident or borderline)

suppressPackageStartupMessages({
  library(pheatmap)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
  library(dplyr)
})

## ---------- edit this path to match your run --------
out_dir <- "/project2/tuannguy_229/RF_Microfluidics/FetalToQuadrato_RF_AllTri_CC_MT_RP_Excl_3kHVG_Wider_Mtry/Model_from_fetal_Second_Third"
## ----------------------------------------------------


pred_rds  <- file.path(out_dir, "pred_labels_organoid.rds")
probs_rds <- file.path(out_dir, "pred_probs_organoid.rds")
obj_rds   <- file.path(out_dir, "quadrato_with_preds.rds")  # Quadrato output file

stopifnot(file.exists(pred_rds), file.exists(probs_rds), file.exists(obj_rds))

pred  <- readRDS(pred_rds)    # named factor: cell -> predicted fetal class
probs <- readRDS(probs_rds)   # matrix: cells x fetal classes, prediction probabilities
quad  <- readRDS(obj_rds)     # Seurat object with predicted_CellType + RF_pred_from_fetal_Second

# Pull original organoid labels for the predicted cells
# Quadrato cells use predicted_CellType (label-transferred), not CellType
org_labels <- as.character(quad@meta.data[names(pred), "predicted_CellType"])
pred_char  <- as.character(pred)

message("N organoid cells evaluated: ", length(pred_char))

# -------------------------------------------------------
# 1. Row-normalized confusion matrix
#    margin = 1 -> rows (organoid types) sum to 1
#    Answers: "given an organoid cell of type X, what was it predicted as?"
# -------------------------------------------------------
conf_tab   <- table(Organoid = org_labels, Predicted_Fetal = pred_char)
conf_rowN  <- prop.table(conf_tab, margin = 1)   # row-normalized

# CSV
write.csv(
  as.data.frame.matrix(round(conf_rowN, 4)),
  file = file.path(out_dir, "eval_confusion_rowNorm_orgType_vs_predFetal.csv")
)
message("Saved: eval_confusion_rowNorm_orgType_vs_predFetal.csv")

# Heatmap
n_fetal <- ncol(conf_rowN)
n_org   <- nrow(conf_rowN)
pdf(file.path(out_dir, "eval_confusion_rowNorm_heatmap.pdf"),
    width = max(7, n_fetal * 1.2), height = max(5, n_org * 0.6),
    useDingbats = FALSE)
pheatmap(
  conf_rowN,
  color            = colorRampPalette(c("white", "#2171b5"))(100),
  display_numbers  = TRUE,
  number_format    = "%.2f",
  number_color     = "black",
  cluster_rows     = TRUE,
  cluster_cols     = TRUE,
  fontsize_number  = 8,
  main             = "Row-normalized: P(predicted fetal | organoid cell type)"
)
dev.off()
message("Saved: eval_confusion_rowNorm_heatmap.pdf")

# -------------------------------------------------------
# 2. Count of organoid cells assigned to each predicted fetal class
# -------------------------------------------------------
pred_counts <- sort(table(Predicted_Fetal = pred_char), decreasing = TRUE)
pred_count_df <- data.frame(
  Predicted_Fetal = names(pred_counts),
  N_organoid_cells = as.integer(pred_counts),
  Pct_of_total = round(100 * as.numeric(pred_counts) / sum(pred_counts), 2)
)
write.csv(pred_count_df,
          file = file.path(out_dir, "eval_counts_organoid_per_predicted_fetalClass.csv"),
          row.names = FALSE)
message("Organoid cells per predicted fetal class:")
print(pred_count_df)

pdf(file.path(out_dir, "eval_counts_per_predFetal_barplot.pdf"),
    width = 7, height = 5, useDingbats = FALSE)
print(
  ggplot(pred_count_df,
         aes(x = reorder(Predicted_Fetal, -N_organoid_cells), y = N_organoid_cells)) +
    geom_bar(stat = "identity", fill = "#2171b5") +
    geom_text(aes(label = N_organoid_cells), vjust = -0.3, size = 3.5) +
    theme_bw() +
    labs(title = "Organoid cells assigned to each predicted fetal class",
         x = "Predicted fetal class", y = "N organoid cells") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
)
dev.off()
message("Saved: eval_counts_per_predFetal_barplot.pdf")

# -------------------------------------------------------
# 3. Count of organoid cells per original organoid label
# -------------------------------------------------------
org_counts <- sort(table(Organoid_CellType = org_labels), decreasing = TRUE)
org_count_df <- data.frame(
  Organoid_CellType = names(org_counts),
  N_cells = as.integer(org_counts),
  Pct_of_total = round(100 * as.numeric(org_counts) / sum(org_counts), 2)
)
write.csv(org_count_df,
          file = file.path(out_dir, "eval_counts_organoid_original_labels.csv"),
          row.names = FALSE)
message("Organoid cells per original label:")
print(org_count_df)

pdf(file.path(out_dir, "eval_counts_original_orgLabels_barplot.pdf"),
    width = 7, height = 5, useDingbats = FALSE)
print(
  ggplot(org_count_df,
         aes(x = reorder(Organoid_CellType, -N_cells), y = N_cells)) +
    geom_bar(stat = "identity", fill = "#6a51a3") +
    geom_text(aes(label = N_cells), vjust = -0.3, size = 3.5) +
    theme_bw() +
    labs(title = "Organoid cells per original cell type label",
         x = "Organoid cell type", y = "N cells") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
)
dev.off()
message("Saved: eval_counts_original_orgLabels_barplot.pdf")

# -------------------------------------------------------
# 4. Max-probability confidence histogram, faceted by organoid cell type
#    Shows whether predictions are confident (prob > 0.7) or ambiguous (< 0.4)
# -------------------------------------------------------
max_probs <- apply(probs, 1, max)
conf_df   <- data.frame(
  Cell          = rownames(probs),
  MaxProb       = max_probs,
  PredFetal     = pred_char,
  Organoid_Type = org_labels,
  stringsAsFactors = FALSE
)

# Summary: mean and median confidence per organoid type
conf_summary <- conf_df %>%
  group_by(Organoid_Type) %>%
  summarise(
    N              = n(),
    Mean_MaxProb   = round(mean(MaxProb), 3),
    Median_MaxProb = round(median(MaxProb), 3),
    Pct_high_conf  = round(100 * mean(MaxProb >= 0.7), 1),
    Pct_low_conf   = round(100 * mean(MaxProb < 0.4), 1),
    .groups        = "drop"
  )
write.csv(conf_summary,
          file = file.path(out_dir, "eval_confidence_summary_by_orgType.csv"),
          row.names = FALSE)
message("Prediction confidence summary:")
print(conf_summary)

pdf(file.path(out_dir, "eval_confidence_histogram_by_orgType.pdf"),
    width = 11, height = 8, useDingbats = FALSE)
print(
  ggplot(conf_df, aes(x = MaxProb, fill = PredFetal)) +
    geom_histogram(bins = 40, color = "white", linewidth = 0.2) +
    facet_wrap(~ Organoid_Type, scales = "free_y") +
    scale_fill_manual(
      values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(unique(pred_char)))
    ) +
    geom_vline(xintercept = 0.4, linetype = "dashed",
               color = "orange", linewidth = 0.6) +
    geom_vline(xintercept = 0.7, linetype = "dashed",
               color = "red", linewidth = 0.6) +
    theme_bw(base_size = 10) +
    labs(
      title = "Max prediction probability per organoid cell (faceted by original cell type)",
      subtitle = "Orange dashed = 0.4 (ambiguous threshold), Red dashed = 0.7 (high confidence)",
      x = "Max class probability", y = "N cells", fill = "Predicted fetal class"
    )
)
dev.off()
message("Saved: eval_confidence_histogram_by_orgType.pdf")

# -------------------------------------------------------
# 5. Combined table: fraction of cells + mean MaxProb per
#    (Organoid_Type x Predicted_Fetal) bin
#    Answers both "where do cells go?" and "how confident are
#    those assignments?"
#    Each row is one organoid-type/predicted-fetal combination.
#    Frac_of_OrgType sums to 1 within each organoid type (row-normalized).
# -------------------------------------------------------
bin_summary <- conf_df %>%
  group_by(Organoid_Type, PredFetal) %>%
  summarise(
    N_cells        = n(),
    Mean_MaxProb   = round(mean(MaxProb), 3),
    Pct_high_conf  = round(100 * mean(MaxProb >= 0.7), 1),
    .groups        = "drop"
  ) %>%
  group_by(Organoid_Type) %>%
  mutate(
    Frac_of_OrgType = round(N_cells / sum(N_cells), 4)
  ) %>%
  ungroup() %>%
  arrange(Organoid_Type, desc(Frac_of_OrgType))

out_csv <- file.path(out_dir, "eval_combined_fraction_and_confidence.csv")
write.csv(bin_summary, file = out_csv, row.names = FALSE)
message("Combined fraction + confidence table:")
print(bin_summary)

# Bubble plot: size = fraction, color = mean MaxProb
# Lets you see both label assignment and confidence in one figure
n_fetal_classes <- length(unique(bin_summary$PredFetal))
n_org_types     <- length(unique(bin_summary$Organoid_Type))
pdf(file.path(out_dir, "eval_combined_fraction_confidence_bubbleplot.pdf"),
    width     = max(7, n_fetal_classes * 1.4),
    height    = max(5, n_org_types * 0.6),
    useDingbats = FALSE)
print(
  ggplot(bin_summary,
         aes(x = PredFetal, y = Organoid_Type,
             size = Frac_of_OrgType, color = Mean_MaxProb)) +
    geom_point() +
    scale_size_continuous(
      name   = "Fraction of\norganoid type",
      range  = c(1, 12)
    ) +
    scale_color_gradientn(
      name   = "Mean max-class\nprobability",
      colors = c("#d9d9d9", "#fc8d59", "#d73027"),
      limits = c(0, 1)
    ) +
    theme_bw(base_size = 11) +
    labs(
      title    = paste0("Label assignment and prediction confidence",
                        " per organoid cell type"),
      subtitle = paste0("Bubble size = fraction of organoid type assigned",
                        " to fetal class; colour = mean max-class probability"),
      x        = "Predicted fetal class",
      y        = "Organoid cell type"
    ) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
)
dev.off()
message("Saved: eval_combined_fraction_confidence_bubbleplot.pdf")

message("\nAll evaluation outputs saved to: ", out_dir)

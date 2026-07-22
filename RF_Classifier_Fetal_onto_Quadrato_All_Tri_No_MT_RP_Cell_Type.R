# Fetal (Wang 2025) -> Quadrato organoid RF  — All Trimesters
# Train on Wang fetal subclass labels, predict onto Quadrato organoid cells.
# Goal: assign fetal cell type identity to each organoid cell / cluster.
# Five models are trained independently:
#   All           — all trimester cells combined; HVGs from all trimester cells
#   First         — first trimester cells only;   HVGs from first trimester cells
#   Second        — second trimester cells only;  HVGs from second trimester cells
#   Third         — third trimester cells only;   HVGs from third trimester cells
#   Second_Third  — second + third trimester cells combined; HVGs from both
# CC/MT/RP genes are excluded before HVG selection in every model.
# Model selection: mtry tuned by OOB error; final model evaluated on OOB
# predictions (equivalent to leave-one-out; no separate validation split needed).
# Downsampling uses stratified random sampling (per-class cap of 10k cells).
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ranger)          # fast multi-threaded RF
  library(ROCR)
  library(pheatmap)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
})

## ----------------- edit paths -----------------
quadrato_rds  <- "/project2/tuannguy_229/RF_Microfluidics/Microfluidics_Arlotta_Dataset/Final_joined_obj_label_transfer_Arlotta.Robj"
fetal_rds     <- "/project2/tuannguy_229/RF_Microfluidics/Fetal_Dataset/fetal_dataset_Wang_2025_all.Robj"
run_root      <- "/project2/tuannguy_229/RF_Microfluidics/FetalToQuadrato_RF_AllTri_CC_MT_RP_Excl_3kHVG_Wider_Mtry"
fetal_label   <- "subclass"   # column in fetal@meta.data used as training labels
max_per_class <- 10000        # stratified cap: max cells sampled per class
seed          <- 123
n_threads     <- 32    # threads for ranger; set to 1 to disable parallelism
n_trees       <- 500          # trees per candidate mtry during tuning
## ---------------------------------------------------

dir.create(run_root, recursive = TRUE, showWarnings = FALSE)

## --- helpers ---

# Stratified downsampling: sample up to max_n cells per class with a fixed seed.
stratified_sample <- function(seurat_obj, max_n, seed = 123) {
  set.seed(seed)
  cells_keep <- unlist(lapply(levels(Idents(seurat_obj)), function(ct) {
    cc <- WhichCells(seurat_obj, idents = ct)
    if (length(cc) > max_n) sample(cc, max_n) else cc
  }), use.names = FALSE)
  subset(seurat_obj, cells = cells_keep)
}

# Train RF with OOB-based mtry tuning.
# Tuning candidates: sqrt(p), p/10, p/3, p/2 (broader search than standard heuristics).
# Best mtry selected by minimum OOB error; final model refit with
# permutation importance enabled for interpretability.
fit_rf <- function(dat, seed = 123,
                   n_threads = parallel::detectCores(), n_trees = 500) {
  set.seed(seed)
  p          <- ncol(dat) - 1L
  mtry_cands <- unique(pmax(1L, floor(c(sqrt(p), p / 10, p / 3, p / 2))))
  message("Tuning mtry over: ", paste(mtry_cands, collapse = ", "),
          "  |  n_threads = ", n_threads)
  message("Training classes: ",
          paste(levels(factor(dat$CellType)), collapse = ", "))
  print(table(dat$CellType))
  
  # Step 1: tune mtry by OOB error (importance = "none" for speed)
  oob_errors <- sapply(mtry_cands, function(m) {
    ranger(
      CellType ~ ., data = dat,
      num.trees   = n_trees,
      mtry        = m,
      probability = TRUE,
      num.threads = n_threads,
      seed        = seed,
      importance  = "none",
      verbose     = FALSE
    )$prediction.error
  })
  
  best_mtry <- mtry_cands[which.min(oob_errors)]
  message(sprintf("Best mtry: %d  |  OOB error: %.4f",
                  best_mtry, min(oob_errors)))
  
  # Step 2: refit best model with permutation importance
  message("Refitting best model with permutation importance...")
  ranger(
    CellType ~ ., data = dat,
    num.trees   = n_trees,
    mtry        = best_mtry,
    probability = TRUE,
    num.threads = n_threads,
    seed        = seed,
    importance  = "permutation",
    verbose     = FALSE
  )
}

# Compute OOB-based per-class AUC and overall accuracy.
# model$predictions contains OOB probability estimates for every training cell
# (available when probability = TRUE in ranger).
# Compute one-vs-rest AUC for a single class given OOB probs + true labels.
compute_auc_one_class <- function(probs_col, true_labels, class_name) {
  true_bin <- as.integer(true_labels == class_name)
  if (sum(true_bin) == 0) return(NA_real_)
  as.numeric(
    ROCR::performance(
      ROCR::prediction(probs_col, true_bin), "auc"
    )@y.values[[1]]
  )
}

oob_metrics <- function(model, dat, out_prefix,
                        n_boot = 1000, ci_level = 0.95, seed = 123) {
  oob_probs   <- model$predictions          # rows = cells, cols = classes
  classes     <- colnames(oob_probs)
  true_labels <- as.character(dat$CellType)
  n           <- nrow(oob_probs)
  class_n     <- table(true_labels)[classes]   # cells per class (for weighting)
  
  # ---- Per-class one-vs-rest OOB AUC (observed) ----------------------------
  auc_vals <- setNames(
    sapply(classes, function(cl)
      compute_auc_one_class(oob_probs[, cl], true_labels, cl)),
    classes
  )
  
  # ---- Bootstrap CIs on per-class AUC --------------------------------------
  # Resample rows (cells) with replacement; recompute AUC each time.
  set.seed(seed)
  boot_aucs <- matrix(NA_real_, nrow = n_boot, ncol = length(classes),
                      dimnames = list(NULL, classes))
  alpha <- 1 - ci_level
  for (b in seq_len(n_boot)) {
    idx <- sample(n, n, replace = TRUE)
    for (cl in classes) {
      boot_aucs[b, cl] <- compute_auc_one_class(
        oob_probs[idx, cl], true_labels[idx], cl
      )
    }
  }
  auc_ci_lo <- apply(boot_aucs, 2, quantile, probs = alpha / 2,
                     na.rm = TRUE)
  auc_ci_hi <- apply(boot_aucs, 2, quantile, probs = 1 - alpha / 2,
                     na.rm = TRUE)
  
  # ---- Macro and weighted mean AUC -----------------------------------------
  valid      <- !is.na(auc_vals)
  macro_auc  <- mean(auc_vals[valid])
  wts        <- as.numeric(class_n[valid])
  weighted_auc <- sum(auc_vals[valid] * wts) / sum(wts)
  message(sprintf(
    "Macro AUC: %.4f  |  Weighted AUC: %.4f  (n_boot=%d, CI=%.0f%%)",
    macro_auc, weighted_auc, n_boot, ci_level * 100
  ))
  
  # ---- ROC plot ------------------------------------------------------------
  pdf(paste0(out_prefix, "_OOB_ROC.pdf"))
  palette_cols <- colorRampPalette(
    RColorBrewer::brewer.pal(8, "Set2")
  )(length(classes))
  first_plotted <- TRUE
  for (i in seq_along(classes)) {
    cl       <- classes[i]
    true_bin <- as.integer(true_labels == cl)
    if (sum(true_bin) == 0) next
    pred_obj <- ROCR::prediction(oob_probs[, cl], true_bin)
    perf_obj <- ROCR::performance(pred_obj, "tpr", "fpr")
    if (first_plotted) {
      plot(perf_obj, main = "OOB ROC Curves (one-vs-rest)",
           col = palette_cols[i])
      first_plotted <- FALSE
    } else {
      plot(perf_obj, add = TRUE, col = palette_cols[i])
    }
  }
  legend("bottomright", legend = classes, col = palette_cols,
         lty = 1, cex = 0.6, bty = "n")
  dev.off()
  
  # ---- Overall OOB accuracy ------------------------------------------------
  pred_class <- classes[max.col(oob_probs, ties.method = "first")]
  oob_acc    <- mean(pred_class == true_labels)
  message(sprintf("OOB accuracy: %.4f  |  OOB error (ranger): %.4f",
                  oob_acc, model$prediction.error))
  
  # ---- Per-class precision, recall, F1 -------------------------------------
  conf_mat <- table(Predicted = pred_class, True = true_labels)
  per_class_stats <- sapply(classes, function(cl) {
    tp   <- if (cl %in% rownames(conf_mat) && cl %in% colnames(conf_mat))
      conf_mat[cl, cl] else 0L
    fp   <- if (cl %in% rownames(conf_mat))
      sum(conf_mat[cl, ]) - tp else 0L
    fn   <- if (cl %in% colnames(conf_mat))
      sum(conf_mat[, cl]) - tp else 0L
    prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    rec  <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    f1   <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0)
      2 * prec * rec / (prec + rec) else NA_real_
    c(Precision = prec, Recall = rec, F1 = f1)
  })
  stats_df          <- as.data.frame(t(per_class_stats))
  stats_df$Class    <- rownames(stats_df)
  stats_df$AUC      <- auc_vals[stats_df$Class]
  stats_df$AUC_CI_lo <- auc_ci_lo[stats_df$Class]
  stats_df$AUC_CI_hi <- auc_ci_hi[stats_df$Class]
  stats_df$N        <- as.integer(class_n[stats_df$Class])
  stats_df <- stats_df[,
                       c("Class", "N", "Precision", "Recall", "F1",
                         "AUC", "AUC_CI_lo", "AUC_CI_hi")
  ]
  write.csv(stats_df, paste0(out_prefix, "_OOB_classStats.csv"),
            row.names = FALSE)
  message("Per-class OOB stats:")
  print(stats_df)
  
  # ---- Summary AUC table ---------------------------------------------------
  summary_auc <- data.frame(
    Metric = c("Macro_AUC", "Weighted_AUC"),
    Value  = c(macro_auc, weighted_auc)
  )
  write.csv(summary_auc, paste0(out_prefix, "_OOB_summaryAUC.csv"),
            row.names = FALSE)
  
  list(auc_vals     = auc_vals,
       auc_ci_lo    = auc_ci_lo,
       auc_ci_hi    = auc_ci_hi,
       macro_auc    = macro_auc,
       weighted_auc = weighted_auc,
       oob_accuracy = oob_acc,
       oob_error    = model$prediction.error,
       class_stats  = stats_df)
}

# Apply RF to test data and return hard labels + probability matrix.
test_rf <- function(model, testdat, n_pred_threads = 1L) {
  res             <- predict(model, data = testdat, num.threads = n_pred_threads)
  probs           <- res$predictions
  rownames(probs) <- rownames(testdat)    # ranger drops rownames; restore here
  pred            <- factor(
    colnames(probs)[max.col(probs, ties.method = "first")],
    levels = colnames(probs)
  )
  names(pred) <- rownames(testdat)
  list(pred = pred, probs = probs)
}

plot_prob_heatmap <- function(probs, pred, out_prefix) {
  ord       <- order(pred)
  probs_ord <- t(probs[ord, , drop = FALSE])
  
  pred_char <- as.character(pred[ord])
  ann_col   <- data.frame(Predicted        = pred_char,
                          row.names        = names(pred)[ord],
                          stringsAsFactors = FALSE)
  
  unique_preds <- unique(pred_char)
  ann_colors   <- list(
    Predicted = setNames(
      colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(unique_preds)),
      unique_preds
    )
  )
  
  png(paste0(out_prefix, "_prob_heatmap.png"),
      width = 900, height = 900, res = 100)
  pheatmap(probs_ord, scale = "none", show_colnames = FALSE,
           annotation_col    = ann_col,
           annotation_colors = ann_colors,
           cluster_rows      = TRUE,
           cluster_cols      = FALSE)
  dev.off()
}

save_importance <- function(model, out_prefix, top_n = 30) {
  imp_df <- data.frame(
    Gene       = names(model$variable.importance),
    Importance = model$variable.importance
  )
  imp_df <- imp_df[order(imp_df$Importance, decreasing = TRUE), ]
  write.csv(imp_df, paste0(out_prefix, "_feature_importance.csv"),
            row.names = FALSE)
  
  top_df <- head(imp_df, top_n)
  pdf(paste0(out_prefix, "_feature_importance_top", top_n, ".pdf"),
      width = 8, height = 6, useDingbats = FALSE)
  print(
    ggplot(top_df, aes(x = reorder(.data$Gene, .data$Importance),
                       y = .data$Importance)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      theme_bw() +
      labs(title = paste0("Top ", top_n, " Permutation Importances"),
           x = "Gene", y = "Permutation Importance")
  )
  dev.off()
}

save_rds <- function(obj, path) {
  saveRDS(obj, path)
  message("Saved: ", path)
}

## ----------- load -----------------
message("Loading data...")
quadrato <- readRDS(quadrato_rds)
fetal    <- readRDS(fetal_rds)

DefaultAssay(quadrato) <- "RNA"
DefaultAssay(fetal)    <- "RNA"

## ------------ prepare fetal object (labels + CellType column) --------------
fetal$CellType <- as.character(fetal@meta.data[[fetal_label]])
fetal$CellType <- gsub("-", "_", fetal$CellType)   # sanitize: dashes -> underscores
fetal <- subset(fetal, subset = !is.na(CellType))
fetal$CellType <- factor(fetal$CellType)
Idents(fetal)  <- "CellType"

message("Fetal class counts (all trimesters, pre-cap):")
print(table(Idents(fetal)))

## -------- Build organoid test matrix (shared across all models) --------
message("Building organoid test matrix (quadrato cells)...")
quadrato_cells <- rownames(quadrato@meta.data)[
  quadrato@meta.data$data.set == "quadrato"
]
if (length(quadrato_cells) == 0) {
  stop("No cells with data.set == 'quadrato' found in quadrato object.")
}
quadrato_sub <- subset(quadrato, cells = quadrato_cells)
organoid_genes <- rownames(GetAssayData(quadrato_sub, assay = "RNA", layer = "data"))
org_celltypes <- as.character(
  quadrato_sub@meta.data[quadrato_cells, "predicted_CellType"]
)

## -------- Define trimester subsets --------
# "All" uses all fetal cells (no trimester filter).
# HVGs for each model are computed from that model's training cells only.
meta <- fetal@meta.data
trimester_cells <- list(
  All          = rownames(meta),
  First        = rownames(meta)[meta$Group == "First_trimester"],
  Second       = rownames(meta)[meta$Group == "Second_trimester"],
  Third        = rownames(meta)[meta$Group == "Third_trimester"],
  Second_Third = rownames(meta)[meta$Group %in% c("Second_trimester", "Third_trimester")]
)

# Check all subsets have cells
for (tri in names(trimester_cells)) {
  n <- length(trimester_cells[[tri]])
  if (n == 0) stop("No cells found for trimester: ", tri,
                   "; check 'Group' column values.")
  message(tri, " trimester cell count: ", n)
}

## -------- Run one model per trimester --------
cc_list <- Seurat::cc.genes.updated.2019
cc_genes <- unique(c(cc_list$s.genes, cc_list$g2m.genes))

for (tri in names(trimester_cells)) {
  
  message("\n", strrep("=", 60))
  message("===== Trimester: ", tri, " =====")
  message(strrep("=", 60))
  
  # Subset fetal to this trimester
  fetal_tri <- subset(fetal, cells = trimester_cells[[tri]])
  fetal_tri$CellType <- droplevels(fetal_tri$CellType)
  Idents(fetal_tri) <- "CellType"
  message("Class counts (pre-cap):")
  print(table(Idents(fetal_tri)))
  
  # --- CC/MT/RP exclusion then 3k HVG selection for this trimester's cells ---
  all_genes     <- rownames(GetAssayData(fetal_tri, assay = "RNA", layer = "data"))
  mt_genes_all  <- grep("^MT-",    all_genes, value = TRUE, ignore.case = TRUE)
  rp_genes_all  <- grep("^RP[SL]", all_genes, value = TRUE, ignore.case = TRUE)
  genes_allowed <- setdiff(all_genes, c(cc_genes, mt_genes_all, rp_genes_all))
  message("Genes after CC/MT/RP exclusion: ", length(genes_allowed))
  
  fetal_filt <- fetal_tri[genes_allowed, ]
  fetal_filt <- FindVariableFeatures(fetal_filt, selection.method = "vst",
                                     nfeatures = 3000, verbose = FALSE)
  hvg_fetal  <- VariableFeatures(fetal_filt)
  
  genes_use <- intersect(hvg_fetal, organoid_genes)
  rm(fetal_filt); gc()   # free gene-filtered copy; HVGs already extracted
  message("HVGs after organoid intersection (CC/MT/RP excluded): ", length(genes_use))
  if (length(genes_use) < 50) {
    warning("Too few genes for trimester ", tri, "; skipping.")
    next
  }
  
  out_dir <- file.path(run_root, paste0("Model_from_fetal_", tri))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  save_rds(genes_use,
           file.path(out_dir, paste0("genes_use_", tri, "_3kHVG_CC_MT_RP_excl.rds")))
  
  # --- Stratified downsampling ---
  fetal_ds <- stratified_sample(fetal_tri, max_n = max_per_class, seed = seed)
  fetal_ds$CellType <- droplevels(fetal_ds$CellType)
  Idents(fetal_ds)  <- "CellType"
  message("Class counts after stratified cap:")
  print(table(Idents(fetal_ds)))
  save_rds(fetal_ds, file.path(out_dir, paste0("fetal_train_obj_capped_", tri, ".rds")))
  
  # --- Build training matrix ---
  train_mat         <- GetAssayData(fetal_ds, assay = "RNA",
                                    layer = "data")[genes_use, , drop = FALSE]
  train_df          <- data.frame(t(as.matrix(train_mat)))
  train_df$CellType <- factor(fetal_ds@active.ident)
  save_rds(train_df, file.path(out_dir, paste0("train_df_", tri, ".rds")))
  rm(fetal_tri, fetal_ds, train_mat); gc()   # free trimester subset + raw matrix
  
  # --- Train RF ---
  message("Training RF for ", tri, " trimester...")
  rf_model <- fit_rf(train_df, seed = seed,
                     n_threads = n_threads, n_trees = n_trees)
  save_rds(rf_model, file.path(out_dir, paste0("rf_model_fetal_", tri, ".rds")))
  
  # --- OOB metrics ---
  oob_res <- oob_metrics(
    rf_model, train_df,
    out_prefix = file.path(out_dir, paste0("oob_", tri)),
    seed = seed
  )
  save_rds(oob_res, file.path(out_dir, paste0("oob_metrics_", tri, ".rds")))
  
  rm(train_df, oob_res); gc()   # free training data + metrics; no longer needed
  
  # --- Feature importance ---
  save_importance(rf_model, out_prefix = file.path(out_dir, paste0("importance_", tri)))
  
  # --- Build organoid test matrix for this gene set ---
  org_mat  <- GetAssayData(quadrato_sub, assay = "RNA",
                           layer = "data")[genes_use, , drop = FALSE]
  test_dat <- data.frame(t(as.matrix(org_mat)))
  
  # --- Predict onto organoid ---
  message("Predicting onto quadrato organoid cells with ", tri, " trimester model...")
  res   <- test_rf(rf_model, test_dat, n_pred_threads = n_threads)
  pred  <- res$pred
  probs <- res$probs
  
  save_rds(pred,  file.path(out_dir, "pred_labels_organoid.rds"))
  save_rds(probs, file.path(out_dir, "pred_probs_organoid.rds"))
  rm(rf_model, res, org_mat, test_dat); gc()   # free model + test matrices
  
  # --- Save organoid object with predictions ---
  quadrato_out <- quadrato_sub
  col_name     <- paste0("RF_pred_from_fetal_", tri)
  quadrato_out@meta.data[[col_name]] <- NA_character_
  quadrato_out@meta.data[names(pred), col_name] <- as.character(pred)
  save_rds(quadrato_out, file.path(out_dir, "quadrato_with_preds.rds"))
  rm(quadrato_out); gc()   # free organoid copy with predictions
  
  # --- Probability heatmap ---
  plot_prob_heatmap(probs, pred,
                    out_prefix = file.path(out_dir, paste0("probHeat_organoid_", tri)))
  
  # --- Confusion: organoid predicted_CellType vs predicted fetal subclass ---
  conf_tab <- table(Organoid_CellType = org_celltypes, Predicted_Fetal = pred)
  write.csv(as.data.frame(prop.table(conf_tab, margin = 2)),
            file = file.path(out_dir, paste0("confusion_organoidCellType_vs_predFetal_", tri, ".csv")),
            row.names = FALSE)
  
  # --- Dot plot ---
  t_dot <- as.data.frame(prop.table(conf_tab, margin = 2))
  pdf(file.path(out_dir, paste0("dotplot_organoidCellType_vs_predFetal_", tri, ".pdf")),
      width = 8, height = 6, useDingbats = FALSE)
  print(
    ggplot(t_dot, aes(Organoid_CellType, Predicted_Fetal)) +
      geom_point(aes(size = Freq, color = Freq)) +
      theme_bw() +
      ggtitle(paste0("Fetal model: ", tri, " trimester -> Quadrato organoid")) +
      xlab("Organoid Cell Type (predicted_CellType)") +
      ylab("Predicted Fetal Subclass") +
      scale_color_gradient(low = "lightgray", high = "darkmagenta",
                           name = "Fraction") +
      scale_size_continuous(name = "Fraction") +
      guides(color = guide_legend(), size = guide_legend()) +
      theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
  )
  dev.off()
  
  message("Saved outputs -> ", out_dir)
}

message("\nDone. All trimester models saved to: ", run_root)

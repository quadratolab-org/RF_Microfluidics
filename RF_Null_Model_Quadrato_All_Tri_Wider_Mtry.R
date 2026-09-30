# Null-model comparison for Fetal -> Quadrato RF classifiers
# Adapted from Uzquiano et al. (2022, Cell) null-model approach:
#   - Permute CellType labels (preserving original class frequencies)
#   - Train a single null RF with the same hyperparameters
#   - Compare OOB accuracy and per-class AUC between real and null
#
# Loads saved train_df_*.rds and rf_model_fetal_*.rds from each trimester.
# No Seurat objects or HVG recomputation needed.
#
# Outputs per trimester (in Model_from_fetal_*/null_model/):
#   - null RF model (.rds)
#   - null OOB per-class stats CSV
#   - real vs. null comparison CSV
#   - OOB ROC overlay plot (real + null)
suppressPackageStartupMessages({
  library(ranger)
  library(ROCR)
  library(ggplot2)
  library(RColorBrewer)
})

## ----------------- edit paths (must match core script) -----------------
run_root  <- "/project2/tuannguy_229/RF_Microfluidics/FetalToQuadrato_RF_AllTri_CC_MT_RP_Excl_3kHVG_Wider_Mtry"
seed      <- 123
n_threads <- 32
n_trees   <- 500
## -----------------------------------------------------------------------

trimesters <- c("All", "First", "Second", "Third", "Second_Third")

## --- helpers (identical to core script) ---

compute_auc_one_class <- function(probs_col, true_labels, class_name) {
  true_bin <- as.integer(true_labels == class_name)
  if (sum(true_bin) == 0) return(NA_real_)
  as.numeric(
    ROCR::performance(
      ROCR::prediction(probs_col, true_bin), "auc"
    )@y.values[[1]]
  )
}

# Train RF with OOB-based mtry tuning (same as core script).
# Tuning candidates: sqrt(p), p/10, p/3, p/2.
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

## --- main loop ---

all_results <- list()

for (tri in trimesters) {
  
  model_dir  <- file.path(run_root, paste0("Model_from_fetal_", tri))
  train_path <- file.path(model_dir, paste0("train_df_", tri, ".rds"))
  real_path  <- file.path(model_dir, paste0("rf_model_fetal_", tri, ".rds"))
  
  if (!file.exists(train_path)) {
    message("Skipping ", tri, ": train_df not found at ", train_path)
    next
  }
  if (!file.exists(real_path)) {
    message("Skipping ", tri, ": real model not found at ", real_path)
    next
  }
  
  message("\n", strrep("=", 60))
  message("===== Null model: ", tri, " =====")
  message(strrep("=", 60))
  
  train_df   <- readRDS(train_path)
  real_model <- readRDS(real_path)
  
  null_dir <- file.path(model_dir, "null_model")
  dir.create(null_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---- Real-model OOB metrics (recomputed from saved model) ----------------
  real_oob_probs  <- real_model$predictions
  classes         <- colnames(real_oob_probs)
  true_labels     <- as.character(train_df$CellType)
  n_cells         <- nrow(train_df)
  real_pred_class <- classes[max.col(real_oob_probs, ties.method = "first")]
  real_oob_acc    <- mean(real_pred_class == true_labels)
  
  real_auc_vals <- setNames(
    sapply(classes, function(cl)
      compute_auc_one_class(real_oob_probs[, cl], true_labels, cl)),
    classes
  )
  real_macro_auc <- mean(real_auc_vals[!is.na(real_auc_vals)])
  
  message(sprintf("Real model — OOB accuracy: %.4f  |  Macro AUC: %.4f",
                  real_oob_acc, real_macro_auc))
  
  # ---- Build null training data (label permutation) -------------------------
  # Permute existing labels, preserving original class frequencies.
  # This is a stricter null than Uzquiano et al.'s uniform sampling:
  # it holds class balance, gene set, and sample size constant, testing
  # only whether the expression-to-label mapping is informative.
  # Appropriate here because our real model trains on asymmetric class
  # frequencies (stratified cap, not equal allocation).
  set.seed(seed)
  null_train            <- train_df
  null_train$CellType   <- sample(train_df$CellType)  # permute, preserving frequencies
  
  message("Null label distribution (permuted, frequencies preserved):")
  print(table(null_train$CellType))
  
  # ---- Train null model (same fit_rf as real model) ------------------------
  message("Training null RF for ", tri, " trimester...")
  null_model <- fit_rf(null_train, seed = seed,
                       n_threads = n_threads, n_trees = n_trees)
  saveRDS(null_model,
          file.path(null_dir, paste0("null_rf_model_", tri, ".rds")))
  message("Saved: ", file.path(null_dir, paste0("null_rf_model_", tri, ".rds")))
  
  # ---- Null-model OOB metrics ----------------------------------------------
  null_oob_probs  <- null_model$predictions
  null_classes    <- colnames(null_oob_probs)
  null_labels     <- as.character(null_train$CellType)
  null_pred_class <- null_classes[max.col(null_oob_probs, ties.method = "first")]
  null_oob_acc    <- mean(null_pred_class == null_labels)
  
  null_auc_vals <- setNames(
    sapply(null_classes, function(cl)
      compute_auc_one_class(null_oob_probs[, cl], null_labels, cl)),
    null_classes
  )
  null_macro_auc <- mean(null_auc_vals[!is.na(null_auc_vals)])
  
  message(sprintf("Null model — OOB accuracy: %.4f  |  Macro AUC: %.4f",
                  null_oob_acc, null_macro_auc))
  
  # ---- Per-class comparison table ------------------------------------------
  per_class_df <- data.frame(
    Class     = classes,
    Real_AUC  = real_auc_vals[classes],
    Null_AUC  = null_auc_vals[classes],
    Delta_AUC = real_auc_vals[classes] - null_auc_vals[classes]
  )
  write.csv(per_class_df,
            file.path(null_dir, paste0("per_class_real_vs_null_AUC_", tri, ".csv")),
            row.names = FALSE)
  message("Per-class AUC comparison:")
  print(per_class_df)
  
  # ---- Summary comparison --------------------------------------------------
  comp_df <- data.frame(
    Metric     = c("OOB_Accuracy", "Macro_AUC"),
    Real       = c(real_oob_acc, real_macro_auc),
    Null       = c(null_oob_acc, null_macro_auc),
    Delta      = c(real_oob_acc - null_oob_acc,
                   real_macro_auc - null_macro_auc)
  )
  write.csv(comp_df,
            file.path(null_dir, paste0("real_vs_null_comparison_", tri, ".csv")),
            row.names = FALSE)
  message("Summary comparison:")
  print(comp_df)
  
  # ---- ROC overlay plot (real vs. null) ------------------------------------
  palette_cols <- colorRampPalette(
    RColorBrewer::brewer.pal(8, "Set2")
  )(length(classes))
  
  pdf(file.path(null_dir, paste0("ROC_real_vs_null_", tri, ".pdf")),
      width = 10, height = 6, useDingbats = FALSE)
  par(mfrow = c(1, 2))
  
  # Real model ROC
  first_plotted <- TRUE
  for (i in seq_along(classes)) {
    cl       <- classes[i]
    true_bin <- as.integer(true_labels == cl)
    if (sum(true_bin) == 0) next
    pred_obj <- ROCR::prediction(real_oob_probs[, cl], true_bin)
    perf_obj <- ROCR::performance(pred_obj, "tpr", "fpr")
    if (first_plotted) {
      plot(perf_obj,
           main = paste0("Real Model OOB ROC (", tri, ")"),
           col = palette_cols[i])
      first_plotted <- FALSE
    } else {
      plot(perf_obj, add = TRUE, col = palette_cols[i])
    }
  }
  abline(a = 0, b = 1, lty = 2, col = "grey50")
  legend("bottomright", legend = classes, col = palette_cols,
         lty = 1, cex = 0.45, bty = "n")
  
  # Null model ROC
  first_plotted <- TRUE
  for (i in seq_along(null_classes)) {
    cl       <- null_classes[i]
    true_bin <- as.integer(null_labels == cl)
    if (sum(true_bin) == 0) next
    pred_obj <- ROCR::prediction(null_oob_probs[, cl], true_bin)
    perf_obj <- ROCR::performance(pred_obj, "tpr", "fpr")
    if (first_plotted) {
      plot(perf_obj,
           main = paste0("Null Model OOB ROC (", tri, ")"),
           col = palette_cols[i])
      first_plotted <- FALSE
    } else {
      plot(perf_obj, add = TRUE, col = palette_cols[i])
    }
  }
  abline(a = 0, b = 1, lty = 2, col = "grey50")
  legend("bottomright", legend = null_classes, col = palette_cols,
         lty = 1, cex = 0.45, bty = "n")
  
  dev.off()
  
  # ---- Bar plot: real vs. null AUC per class -------------------------------
  bar_df <- data.frame(
    Class = rep(classes, 2),
    AUC   = c(real_auc_vals[classes], null_auc_vals[classes]),
    Model = rep(c("Real", "Null"), each = length(classes))
  )
  bar_df$Model <- factor(bar_df$Model, levels = c("Real", "Null"))
  
  pdf(file.path(null_dir, paste0("barplot_AUC_real_vs_null_", tri, ".pdf")),
      width = max(8, length(classes) * 0.6), height = 6, useDingbats = FALSE)
  print(
    ggplot(bar_df, aes(x = Class, y = AUC, fill = Model)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8),
               width = 0.7) +
      scale_fill_manual(values = c(Real = "steelblue", Null = "grey60")) +
      theme_bw() +
      labs(title = paste0("Per-class OOB AUC: Real vs. Null (", tri, ")"),
           x = "Cell Type", y = "One-vs-Rest AUC") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "red",
                 linewidth = 0.5)
  )
  dev.off()
  
  all_results[[tri]] <- comp_df
  rm(train_df, null_train, real_model, null_model,
     real_oob_probs, null_oob_probs); gc()
  message("Saved null-model outputs -> ", null_dir)
}

## --- Combined summary across all trimesters ---
if (length(all_results) > 0) {
  combined <- do.call(rbind, lapply(names(all_results), function(tri) {
    df <- all_results[[tri]]
    df$Trimester <- tri
    df
  }))
  combined <- combined[, c("Trimester", "Metric", "Real", "Null", "Delta")]
  write.csv(combined,
            file.path(run_root, "null_model_summary_all_trimesters.csv"),
            row.names = FALSE)
  message("\nCombined null-model summary:")
  print(combined)
  message("Saved: ", file.path(run_root, "null_model_summary_all_trimesters.csv"))
}

message("\nDone. Null-model analysis complete.")

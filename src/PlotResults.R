#!/usr/bin/env Rscript

# First figure: compares the "core" (ground+exp+fake) train-test split results, using
# the full ("All") feature set only, across all three modeling pipelines - the six
# classical ML models (./data/ML), the ResNet model (./data/DL), and NetTCR-2.0
# (./data/NetTCR). Each ML/ResNet bar is the mean of its 5-fold CV test scores, with an
# error bar (+-1 SD) and the individual per-fold points overlaid; NetTCR is a single run
# (no CV), so it's a bar with no error bar or points. Saved to ./Results.

# Auto-install any required package that isn't already available, then load it (except
# svglite/systemfonts/fontLiberation, which are only ever used via `::` below, so they
# just need to be installed, not attached).
required_packages <- c(
  "dplyr", "tidyr", "readr", "ggplot2", "forcats",
  "svglite", "systemfonts", "fontLiberation", "factoextra", "patchwork"
)
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(forcats)
library(factoextra)
library(patchwork)

# Arial itself isn't a font this (or most Linux) systems have installed, so without this
# svglite silently substitutes a different font while still labeling text family="Arial"
# requested below - use Liberation Sans (a free, metric-compatible Arial substitute) in
# its place. register_font() makes `family = "Arial"` resolve correctly for layout on
# any graphics device; passing the same font files as `user_fonts` to each svglite
# ggsave() call below additionally makes the saved SVG literally declare
# font-family="Arial" (not "Liberation Sans") - that part only happens if a machine
# opening the SVG for editing has actual Arial installed (typical on Windows/Mac), which
# plain register_font() alone does not achieve since svglite otherwise writes the font
# file's own internal name into the SVG.
liberation_dir <- system.file("fonts", "liberation-fonts", package = "fontLiberation")
arial_faces <- list(
  plain = file.path(liberation_dir, "LiberationSans-Regular.ttf"),
  bold = file.path(liberation_dir, "LiberationSans-Bold.ttf"),
  italic = file.path(liberation_dir, "LiberationSans-Italic.ttf"),
  bolditalic = file.path(liberation_dir, "LiberationSans-BoldItalic.ttf")
)
do.call(systemfonts::register_font, c(list(name = "Arial"), arial_faces))
ARIAL_USER_FONTS <- list(Arial = arial_faces)

ML_DIR <- "./data/ML/Results/Base/test"
ML_FOLDS <- file.path(ML_DIR, "All", "test_fold_metrics.csv")
DL_FOLDS <- "./data/DL/results.csv"
NETTCR_SUMMARY <- "./data/NetTCR/Comparisons/nettcr_summary.csv"
OUT_DIR <- "./Results"

# One subfolder per figure/analysis, so ./Results isn't one flat pile of ~30 files.
DIR_MODEL_COMPARISON <- file.path(OUT_DIR, "01_ModelComparison")
DIR_FEATURE_COMPARISON <- file.path(OUT_DIR, "02_FeatureComparison")
DIR_FEATURE_IMPORTANCE <- file.path(OUT_DIR, "03_FeatureImportance")
DIR_PCA <- file.path(OUT_DIR, "04_PCA")
DIR_CV_BIAS <- file.path(OUT_DIR, "05_CVBias")
DIR_GENERALIZATION <- file.path(OUT_DIR, "06_Generalization")
DIR_IMMREP_BALANCED <- file.path(OUT_DIR, "07_ImmrepBalanced")
DIR_SAMPLE_HEATMAP <- file.path(OUT_DIR, "08_SampleHeatmap")

for (d in c(DIR_MODEL_COMPARISON, DIR_FEATURE_COMPARISON, DIR_FEATURE_IMPORTANCE, DIR_PCA, DIR_CV_BIAS, DIR_GENERALIZATION, DIR_IMMREP_BALANCED, DIR_SAMPLE_HEATMAP)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

metric_labels <- c(
  Accuracy = "Accuracy", ROC_AUC = "ROC-AUC", PR_AUC = "PR-AUC",
  PPV = "PPV (positive set accuracy)", NPV = "NPV (negative set accuracy)"
)
method_colors <- c(ML = "#1f77b4", DL = "#d62728", NetTCR = "#2ca02c", ERGO = "#8c564b", Ensemble = "#333333")

# Abbreviations for the six ML model names, so axis labels are shorter and easier to
# read; ResNet/NetTCR-2.0 are already short enough as-is.
ml_model_abbrev <- c(
  "Logistic Regression" = "LR",
  "Random Forest" = "RF",
  "Gradient Boosting" = "GB",
  "Support Vector Machine" = "SVM",
  "K-Nearest Neighbors" = "KNN",
  "Multilayer Perceptron" = "MLP"
)

# ---------------------------------------------------------------------------
# Per-fold points: the six classical ML models (All feature set) + the ResNet model
# ---------------------------------------------------------------------------
ml_points <- read_csv(ML_FOLDS, show_col_types = FALSE) %>%
  transmute(Method = "ML", Model = unname(ml_model_abbrev[Model]), Fold = as.character(Fold),
            Accuracy, ROC_AUC, PR_AUC, PPV, NPV)

dl_points <- read_csv(DL_FOLDS, show_col_types = FALSE) %>%
  filter(Run == "core_all") %>%
  transmute(Method = "DL", Model = "ResNet", Fold = Fold,
            Accuracy = Test_Accuracy, ROC_AUC = Test_ROC_AUC, PR_AUC = Test_PR_AUC,
            # Test_Precision is PPV under another name; results.csv has no NPV (no raw
            # TN/FN counts to compute it from) so it's left NA - the ResNet bar is
            # simply absent from the NPV plot rather than showing a fabricated value.
            PPV = Test_Precision, NPV = NA_real_)

points <- bind_rows(ml_points, dl_points) %>%
  pivot_longer(cols = c(Accuracy, ROC_AUC, PR_AUC, PPV, NPV), names_to = "Metric", values_to = "Value")

# Bars = mean/SD computed directly from the same per-fold points plotted on top of them
bars <- points %>%
  group_by(Method, Model, Metric) %>%
  summarise(Mean = mean(Value, na.rm = TRUE), SD = sd(Value, na.rm = TRUE), .groups = "drop")

# ---------------------------------------------------------------------------
# NetTCR-2.0: a single run, no CV - bar only, no error bar or points
# ---------------------------------------------------------------------------
nettcr_bar <- read_csv(NETTCR_SUMMARY, show_col_types = FALSE) %>%
  filter(Comparison == "core_80_20") %>%
  transmute(Method = "NetTCR", Model = "NetTCR-2.0", Accuracy, ROC_AUC, PR_AUC, PPV, NPV) %>%
  pivot_longer(cols = c(Accuracy, ROC_AUC, PR_AUC, PPV, NPV), names_to = "Metric", values_to = "Mean") %>%
  mutate(SD = NA_real_)

bars <- bind_rows(bars, nettcr_bar) %>%
  mutate(
    Metric = factor(Metric, levels = names(metric_labels), labels = metric_labels),
    Method = factor(Method, levels = c("ML", "DL", "NetTCR"))
  )
points <- points %>%
  mutate(Metric = factor(Metric, levels = names(metric_labels), labels = metric_labels))

write_csv(bars, file.path(DIR_MODEL_COMPARISON, "core_all_summary.csv"))
write_csv(points, file.path(DIR_MODEL_COMPARISON, "core_all_folds.csv"))

plot_metric <- function(metric_name, out_file) {
  label <- metric_labels[[metric_name]]
  bar_df <- bars %>% filter(Metric == label)
  point_df <- points %>% filter(Metric == label)

  model_order <- bar_df %>% arrange(desc(Mean)) %>% pull(Model)
  bar_df <- bar_df %>% mutate(Model = factor(Model, levels = model_order))
  point_df <- point_df %>% mutate(Model = factor(Model, levels = model_order))

  p <- ggplot(bar_df, aes(x = Model, y = Mean, fill = Method)) +
    geom_col(width = 0.65) +
    geom_errorbar(
      aes(ymin = pmax(Mean - SD, 0), ymax = pmin(Mean + SD, 1)),
      width = 0.2, na.rm = TRUE
    ) +
    geom_jitter(
      data = point_df, aes(x = Model, y = Value),
      inherit.aes = FALSE, width = 0.12, height = 0, size = 1.8, alpha = 0.7, shape = 21,
      fill = "white", color = "black", stroke = 0.4
    ) +
    scale_fill_manual(values = method_colors) +
    scale_y_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
    labs(x = NULL, y = label) +
    theme_minimal(base_size = 13) +
    theme(
      text = element_text(color = "black", size = 20, family = "Arial"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )

  ggsave(
    file.path(DIR_MODEL_COMPARISON, out_file), p, width = 8, height = 6,
    device = svglite::svglite, user_fonts = ARIAL_USER_FONTS
  )
}

plot_metric("Accuracy", "core_accuracy.svg")
plot_metric("ROC_AUC", "core_roc_auc.svg")
plot_metric("PR_AUC", "core_pr_auc.svg")
plot_metric("PPV", "core_ppv.svg")
plot_metric("NPV", "core_npv.svg")

# =============================================================================
# Second figure: feature set comparison on Core - all 6 ML models + ResNet (DL),
# faceted by feature set (All/Seq/Str/Phys, plus DL's feature-free "Pos" baseline,
# which has no ML equivalent). NetTCR isn't included here: it has no feature-set
# variants (it's always raw CDR3 sequence, not one of these engineered sets).
# =============================================================================
feature_set_map <- c(All = "All", Seq = "Seq", Str = "Str", Phys = "Phys")

ml_feature_points <- lapply(names(feature_set_map), function(fs) {
  read_csv(file.path(ML_DIR, fs, "test_fold_metrics.csv"), show_col_types = FALSE) %>%
    transmute(Method = "ML", Model = unname(ml_model_abbrev[Model]), FeatureSet = feature_set_map[[fs]],
              Fold = as.character(Fold), Accuracy, ROC_AUC, PR_AUC, PPV, NPV)
}) %>% bind_rows()

run_to_featureset <- c(
  core_all = "All", core_seq = "Seq", core_str = "Str",
  core_phys = "Phys", core_pos = "Pos"
)
dl_feature_points <- read_csv(DL_FOLDS, show_col_types = FALSE) %>%
  filter(Run %in% names(run_to_featureset)) %>%
  transmute(Method = "DL", Model = "ResNet", FeatureSet = run_to_featureset[Run], Fold = Fold,
            Accuracy = Test_Accuracy, ROC_AUC = Test_ROC_AUC, PR_AUC = Test_PR_AUC,
            PPV = Test_Precision, NPV = NA_real_)

feature_points <- bind_rows(ml_feature_points, dl_feature_points) %>%
  pivot_longer(cols = c(Accuracy, ROC_AUC, PR_AUC, PPV, NPV), names_to = "Metric", values_to = "Value")

feature_bars <- feature_points %>%
  group_by(Method, Model, FeatureSet, Metric) %>%
  summarise(Mean = mean(Value, na.rm = TRUE), SD = sd(Value, na.rm = TRUE), .groups = "drop")

# NetTCR-2.0 is sequence-only (CDR3a/CDR3b/peptide, no structure/physics features), so
# it's dropped into the "Seq" facet alongside the Seq-only ML/ResNet models - a single
# run (no CV), so no error bar or points, same as in the first figure.
nettcr_feature_bar <- read_csv(NETTCR_SUMMARY, show_col_types = FALSE) %>%
  filter(Comparison == "core_80_20") %>%
  transmute(Method = "NetTCR", Model = "NetTCR-2.0", FeatureSet = "Seq",
            Accuracy, ROC_AUC, PR_AUC, PPV, NPV) %>%
  pivot_longer(cols = c(Accuracy, ROC_AUC, PR_AUC, PPV, NPV), names_to = "Metric", values_to = "Mean") %>%
  mutate(SD = NA_real_)

feature_bars <- bind_rows(feature_bars, nettcr_feature_bar) %>%
  mutate(
    Metric = factor(Metric, levels = names(metric_labels), labels = metric_labels),
    FeatureSet = factor(FeatureSet, levels = c("All", "Seq", "Str", "Phys", "Pos")),
    Method = factor(Method, levels = c("ML", "DL", "NetTCR"))
  )
feature_points <- feature_points %>%
  mutate(
    Metric = factor(Metric, levels = names(metric_labels), labels = metric_labels),
    FeatureSet = factor(FeatureSet, levels = c("All", "Seq", "Str", "Phys", "Pos"))
  )

write_csv(feature_bars, file.path(DIR_FEATURE_COMPARISON, "core_features_summary.csv"))
write_csv(feature_points, file.path(DIR_FEATURE_COMPARISON, "core_features_folds.csv"))

# Model order held fixed across every facet (by performance on "All", where every ML/DL
# model has a value) so the same model lines up in the same position panel to panel;
# NetTCR-2.0 only ever appears in the Seq facet, so it's just appended at the end.
feature_model_order <- feature_bars %>%
  filter(FeatureSet == "All", Metric == "Accuracy") %>%
  arrange(desc(Mean)) %>%
  pull(Model) %>%
  c("NetTCR-2.0")

plot_feature_metric <- function(metric_name, out_file) {
  label <- metric_labels[[metric_name]]
  bar_df <- feature_bars %>% filter(Metric == label) %>%
    mutate(Model = factor(Model, levels = feature_model_order))
  point_df <- feature_points %>% filter(Metric == label) %>%
    mutate(Model = factor(Model, levels = feature_model_order))

  p <- ggplot(bar_df, aes(x = Model, y = Mean, fill = Method)) +
    geom_col(width = 0.65) +
    geom_errorbar(
      aes(ymin = pmax(Mean - SD, 0), ymax = pmin(Mean + SD, 1)),
      width = 0.2, na.rm = TRUE
    ) +
    geom_jitter(
      data = point_df, aes(x = Model, y = Value),
      inherit.aes = FALSE, width = 0.12, height = 0, size = 1.5, alpha = 0.7, shape = 21,
      fill = "white", color = "black", stroke = 0.4
    ) +
    facet_wrap(~FeatureSet, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = method_colors) +
    scale_y_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
    labs(x = NULL, y = label) +
    theme_minimal(base_size = 13) +
    theme(
      text = element_text(color = "black", size = 16, family = "Arial"),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 40, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )

  ggsave(
    file.path(DIR_FEATURE_COMPARISON, out_file), p, width = 16, height = 6,
    device = svglite::svglite, user_fonts = ARIAL_USER_FONTS
  )
}

plot_feature_metric("Accuracy", "core_features_accuracy.svg")
plot_feature_metric("ROC_AUC", "core_features_roc_auc.svg")
plot_feature_metric("PR_AUC", "core_features_pr_auc.svg")
plot_feature_metric("PPV", "core_features_ppv.svg")
plot_feature_metric("NPV", "core_features_npv.svg")

# =============================================================================
# Third figure: ML feature importance (Gradient Boosting - the consistent top
# performer throughout - native/impurity-based importance on the "All" feature set),
# and a PCA of the core dataset using all features vs only the top-importance ones.
# =============================================================================
FEATURES_CSV <- "./data/Features/features.csv"
ANNOTATION_CSV <- "./data/ML/Annotation.csv"
TOP_N_FEATURES <- 10
BIPLOT_N_VARS <- 5  # only the top-contributing variables get an arrow drawn in the biplot
IMPORTANCE_MODEL <- "Gradient Boosting"

# --- Feature importance ---
importance <- read_csv(
  file.path(ML_DIR, "All", "validation_feature_importance_summary.csv"), show_col_types = FALSE
) %>%
  filter(Model == IMPORTANCE_MODEL, Method == "native") %>%
  arrange(desc(Importance_Mean))

top_features <- importance %>% slice_head(n = TOP_N_FEATURES)
write_csv(top_features, file.path(DIR_FEATURE_IMPORTANCE, "core_feature_importance_top.csv"))

p_importance <- ggplot(
  top_features %>% mutate(Feature = fct_reorder(Feature, Importance_Mean)),
  aes(x = Feature, y = Importance_Mean)
) +
  geom_col(fill = method_colors[["ML"]], width = 0.7) +
  geom_errorbar(
    aes(ymin = pmax(Importance_Mean - Importance_Std, 0), ymax = Importance_Mean + Importance_Std),
    width = 0.25
  ) +
  coord_flip() +
  labs(x = NULL, y = "Importance") +
  theme_minimal(base_size = 13) +
  theme(
    text = element_text(color = "black", size = 20, family = "Arial"),
    plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 60),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(DIR_FEATURE_IMPORTANCE, "core_feature_importance.svg"), p_importance, width = 8, height = 4,
  device = svglite::svglite, user_fonts = ARIAL_USER_FONTS
)

# --- PCA: all features vs top-importance features, core dataset ---
# Annotation.csv was saved with an implicit leading row-index column that has no header
# entry for it (the same convention pandas' read_csv auto-detects, but readr does not -
# readr::read_csv() would silently shift every column over by one here). Base
# read.csv() handles it correctly, same as pandas.
annotation <- read.csv(ANNOTATION_CSV, stringsAsFactors = FALSE)
feature_matrix <- read_csv(FEATURES_CSV, show_col_types = FALSE)

core_data <- annotation %>%
  filter(source %in% c("ground", "exp", "fake")) %>%
  inner_join(feature_matrix, by = "name")

feature_cols <- setdiff(names(feature_matrix), "name")

# Median-impute each feature column, same strategy evaluate_models() uses for
# training; a column that's entirely NaN (e.g. MHCb-chain measurements for class I
# samples) has no median to impute with, so it's filled with 0 instead - it'll end up
# constant and get dropped just below, same as it effectively is in the models.
impute_column <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(0, length(x)))
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}

X <- as.data.frame(lapply(core_data[feature_cols], impute_column))
names(X) <- feature_cols
# Drop zero-variance columns - scale() divides by SD, undefined/0 for these.
nonconstant <- vapply(X, function(col) sd(col) > 0, logical(1))
X <- X[, nonconstant, drop = FALSE]

meta <- core_data %>%
  transmute(
    name = name,
    Label = factor(response, levels = c(0, 1), labels = c("Non-binder", "Binder")),
    Source = source
  )

pca_all <- prcomp(X, center = TRUE, scale. = TRUE)

top_feature_names <- intersect(top_features$Feature, names(X))
pca_top <- prcomp(X[, top_feature_names, drop = FALSE], center = TRUE, scale. = TRUE)

extract_scores <- function(pca) {
  var_explained <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
  scores <- as.data.frame(pca$x[, 1:2])
  names(scores) <- c("PC1", "PC2")
  scores$VarExplainedPC1 <- round(var_explained[1], 1)
  scores$VarExplainedPC2 <- round(var_explained[2], 1)
  bind_cols(meta, scores)
}

write_csv(extract_scores(pca_all), file.path(DIR_PCA, "core_pca_all_features.csv"))
write_csv(extract_scores(pca_top), file.path(DIR_PCA, "core_pca_top_features.csv"))

# Individuals only (all ~650 features) - a biplot of every feature's loading arrow
# would be unreadable, so just the PCA scores, colored by Label and shaped by Source.
p_pca_all <- fviz_pca_ind(
  pca_all, geom.ind = "point", label = "none",
  col.ind = meta$Label, shape.ind = meta$Source,
  pointsize = 2, alpha.ind = 0.8
) +
  labs(title = NULL, color = "Label", shape = "Source") +
  theme_minimal(base_size = 13) +
  theme(
    text = element_text(color = "black", size = 20, family = "Arial"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )
ggsave(
  file.path(DIR_PCA, "core_pca_all_features.svg"), p_pca_all, width = 6, height = 5,
  device = svglite::svglite, user_fonts = ARIAL_USER_FONTS
)

# Biplot (individuals + variable loading arrows) for the reduced, top-importance-feature
# PCA specifically. The PCA itself is still computed from all TOP_N_FEATURES variables,
# but only the BIPLOT_N_VARS with the highest contribution to these two PCs get an
# arrow drawn, so the biplot stays readable instead of 20 overlapping arrows/labels.
p_pca_top <- fviz_pca_biplot(
  pca_top, geom.ind = "point", label = "var", repel = TRUE, labelsize = 3,
  col.ind = meta$Label, shape.ind = meta$Source, col.var = "black",
  pointsize = 2, alpha.ind = 0.7,
  select.var = list(contrib = BIPLOT_N_VARS)
) +
  labs(title = NULL, color = "Label", shape = "Source") +
  theme_minimal(base_size = 13) +
  theme(
    text = element_text(color = "black", size = 20, family = "Arial"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )
ggsave(
  file.path(DIR_PCA, "core_pca_top_features.svg"), p_pca_top, width = 6, height = 5,
  device = svglite::svglite, user_fonts = ARIAL_USER_FONTS
)

# =============================================================================
# Fourth figure: regular (random) CV on Core vs the non-overlapping, grouped CV runs
# (evaluate_models' bias_column) by MHC type/allele, Peptide, and TRBV gene - each
# keeps every sample sharing a value of that column in a single fold, so the validation
# estimate can't be inflated by near-identical examples leaking between train and
# validation. Faceted by model (the 6 ML models + ResNet), x = which CV this validation
# estimate came from. Uses the *validation* fold metrics (not test) specifically, since
# the whole point is how the within-CV estimate itself shifts once leakage across that
# variable is prevented - the held-out test set is unaffected by how the CV folds were
# built either way. DL also ran a TRAV-grouped CV, but ML never did, so it's left out
# here to keep every model directly comparable on the same set of CV types.
# Effects/Test is excluded: that output directory got overwritten by two different
# bias_column runs sharing it (it carries both a Type and a Peptide fold-composition
# file), so its actual metrics only reflect whichever ran last.
# =============================================================================
cv_bias_sources <- c(
  Core = "./data/ML/Results/Base/test/All/validation_fold_metrics.csv",
  MHCtype = "./data/ML/Results/Effects/MHCtype/validation_fold_metrics.csv",
  Peptide = "./data/ML/Results/Effects/Peptide/validation_fold_metrics.csv",
  TRBV = "./data/ML/Results/Effects/TRBV/validation_fold_metrics.csv"
)
cv_type_colors <- c(
  Core = "#7f7f7f", MHCtype = "#1f77b4", Allele = "#2ca02c",
  Peptide = "#ff7f0e", TRBV = "#9467bd"
)
cv_type_levels <- c("Core", "MHCtype", "Allele", "Peptide", "TRBV")
cv_bias_model_order <- c(unname(ml_model_abbrev), "ResNet")

ml_cv_bias_points <- lapply(names(cv_bias_sources), function(cv_type) {
  read_csv(cv_bias_sources[[cv_type]], show_col_types = FALSE) %>%
    transmute(CVType = cv_type, Model = unname(ml_model_abbrev[Model]), Fold = as.character(Fold),
              Accuracy, ROC_AUC, PR_AUC, PPV, NPV)
}) %>% bind_rows()

dl_run_to_cvtype <- c(
  core_all = "Core", bias_allele = "Allele", bias_peptide = "Peptide", bias_trbv = "TRBV"
)
dl_cv_bias_points <- read_csv(DL_FOLDS, show_col_types = FALSE) %>%
  filter(Run %in% names(dl_run_to_cvtype)) %>%
  transmute(CVType = dl_run_to_cvtype[Run], Model = "ResNet", Fold = Fold,
            Accuracy = Val_Accuracy, ROC_AUC = Val_ROC_AUC, PR_AUC = Val_PR_AUC,
            PPV = Val_Precision, NPV = NA_real_)

cv_bias_points <- bind_rows(ml_cv_bias_points, dl_cv_bias_points) %>%
  pivot_longer(cols = c(Accuracy, ROC_AUC, PR_AUC, PPV, NPV), names_to = "Metric", values_to = "Value")

cv_bias_bars <- cv_bias_points %>%
  group_by(CVType, Model, Metric) %>%
  summarise(Mean = mean(Value, na.rm = TRUE), SD = sd(Value, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Metric = factor(Metric, levels = names(metric_labels), labels = metric_labels),
    CVType = factor(CVType, levels = cv_type_levels),
    Model = factor(Model, levels = cv_bias_model_order)
  )
cv_bias_points <- cv_bias_points %>%
  mutate(
    Metric = factor(Metric, levels = names(metric_labels), labels = metric_labels),
    CVType = factor(CVType, levels = cv_type_levels),
    Model = factor(Model, levels = cv_bias_model_order)
  )

write_csv(cv_bias_bars, file.path(DIR_CV_BIAS, "core_cv_bias_summary.csv"))
write_csv(cv_bias_points, file.path(DIR_CV_BIAS, "core_cv_bias_folds.csv"))

plot_cv_bias_metric <- function(metric_name, out_file) {
  label <- metric_labels[[metric_name]]
  bar_df <- cv_bias_bars %>% filter(Metric == label)
  point_df <- cv_bias_points %>% filter(Metric == label)

  p <- ggplot(bar_df, aes(x = CVType, y = Mean, fill = CVType)) +
    geom_col(width = 0.65) +
    geom_errorbar(
      aes(ymin = pmax(Mean - SD, 0), ymax = pmin(Mean + SD, 1)),
      width = 0.2, na.rm = TRUE
    ) +
    geom_jitter(
      data = point_df, aes(x = CVType, y = Value),
      inherit.aes = FALSE, width = 0.12, height = 0, size = 1.5, alpha = 0.7, shape = 21,
      fill = "white", color = "black", stroke = 0.4
    ) +
    facet_wrap(~Model, nrow = 2, scales = "free_x") +
    scale_fill_manual(values = cv_type_colors) +
    scale_y_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
    labs(x = NULL, y = label, fill = "CV grouped by") +
    theme_minimal(base_size = 13) +
    theme(
      text = element_text(color = "black", size = 14, family = "Arial"),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )

  ggsave(
    file.path(DIR_CV_BIAS, out_file), p, width = 11, height = 8,
    device = svglite::svglite, user_fonts = ARIAL_USER_FONTS
  )
}

plot_cv_bias_metric("Accuracy", "core_cv_bias_accuracy.svg")
plot_cv_bias_metric("ROC_AUC", "core_cv_bias_roc_auc.svg")
plot_cv_bias_metric("PR_AUC", "core_cv_bias_pr_auc.svg")
plot_cv_bias_metric("PPV", "core_cv_bias_ppv.svg")
plot_cv_bias_metric("NPV", "core_cv_bias_npv.svg")

# =============================================================================
# Fifth figure: generalization to ImmRep - two comparisons, each its own plot (not
# faceted). "Core to ImmRep" (train on Core, test on the structured ImmRep subset) puts
# our core-trained models (ML/DL/NetTCR) SIDE BY SIDE with models that never saw our
# Core data at all - NetTCR-2.0 (paper) (retrained only on NetTCR-2.0's own published
# training set) and ERGO-II's two published variants (Auto/LSTM) - both tested on the
# complete ImmRep set rather than just the structured subset, so this is an approximate,
# not exact, like-for-like comparison. "Core+ImmRep (train/test)" (Core+ImmRep combined,
# fixed split) has no such externally-trained counterpart (nothing else was evaluated on
# that exact split), so it stays core-trained models only. Each pipeline/variant labels
# these differently in its own results csv:
#                          Core -> ImmRep              Core+ImmRep (train/test)
#   ML     (Split)        non_immpred_to_immpred       FullRun
#   DL     (Run)            core_to_immrep               core_immrep
#   NetTCR (Comparison)    core_vs_immrep, paper_train_vs_immrep_full   core_immrep_80_20
#   ERGO   (Comparison)    immrep_full                  -
# ML's "FullRun" (a single fixed split matching the DL pipeline's own FullRun
# comparison) and every NetTCR/ERGO run are a single fit with no inner CV; giving them a
# "Fold" of just one value makes sd() return NA for them automatically (same as a
# literal single-element sd() call), so they fall out of the same bar/error-bar/point
# pipeline as the 5-fold CV runs (ML's non_immpred_to_immpred, both DL runs) without
# any special-cased code path, and simply show no error bar/points.
# =============================================================================
ERGO_SUMMARY <- "./data/ERGO/ergo_summary.csv"

generalization_specs <- list(
  "Core to ImmRep" = list(
    ml_file = "./data/ML/Results/Base/non_immpred_to_immpred/All/test_fold_metrics.csv",
    dl_run = "core_to_immrep",
    nettcr_comparisons = c("NetTCR-2.0 (ours)" = "core_vs_immrep", "NetTCR-2.0 (paper)" = "paper_train_vs_immrep_full"),
    ergo_comparison = "immrep_full"
  ),
  "Core+ImmRep (train/test)" = list(
    ml_file = "./data/ML/Results/FullRun/FullRun/All/test_fold_metrics.csv",
    dl_run = "core_immrep",
    nettcr_comparisons = c("NetTCR-2.0 (ours)" = "core_immrep_80_20"),
    ergo_comparison = NULL,
    ensemble_summary = "./data/Results/fullrun_sample_outcomes_summary.csv"
  )
)
generalization_levels <- names(generalization_specs)

generalization_points <- lapply(generalization_levels, function(comparison) {
  spec <- generalization_specs[[comparison]]
  parts <- list()

  if (!is.null(spec$ml_file)) {
    parts$ml <- read_csv(spec$ml_file, show_col_types = FALSE) %>%
      transmute(Method = "ML", Model = unname(ml_model_abbrev[Model]), Fold = as.character(Fold),
                Accuracy, ROC_AUC, PR_AUC, PPV, NPV)
  }

  if (!is.null(spec$dl_run)) {
    parts$dl <- read_csv(DL_FOLDS, show_col_types = FALSE) %>%
      filter(Run == spec$dl_run) %>%
      transmute(Method = "DL", Model = "ResNet", Fold = Fold,
                Accuracy = Test_Accuracy, ROC_AUC = Test_ROC_AUC, PR_AUC = Test_PR_AUC,
                PPV = Test_Precision, NPV = NA_real_)
  }

  if (!is.null(spec$nettcr_comparisons)) {
    parts$nettcr <- lapply(names(spec$nettcr_comparisons), function(model_label) {
      read_csv(NETTCR_SUMMARY, show_col_types = FALSE) %>%
        filter(Comparison == spec$nettcr_comparisons[[model_label]]) %>%
        transmute(Method = "NetTCR", Model = model_label, Fold = "1",
                  Accuracy, ROC_AUC, PR_AUC, PPV, NPV)
    }) %>% bind_rows()
  }

  if (!is.null(spec$ergo_comparison)) {
    parts$ergo <- read_csv(ERGO_SUMMARY, show_col_types = FALSE) %>%
      filter(Comparison == spec$ergo_comparison) %>%
      transmute(Method = "ERGO", Model = Model, Fold = "1",
                Accuracy, ROC_AUC, PR_AUC, PPV, NPV)
  }

  if (!is.null(spec$ensemble_summary)) {
    parts$ensemble <- read_csv(spec$ensemble_summary, show_col_types = FALSE) %>%
      filter(Model == "Ensemble") %>%
      transmute(Method = "Ensemble", Model = "Ensemble", Fold = "1",
                Accuracy, ROC_AUC, PR_AUC, PPV, NPV)
  }

  bind_rows(parts) %>% mutate(Comparison = comparison)
}) %>%
  bind_rows() %>%
  pivot_longer(cols = c(Accuracy, ROC_AUC, PR_AUC, PPV, NPV), names_to = "Metric", values_to = "Value")

generalization_bars <- generalization_points %>%
  group_by(Comparison, Method, Model, Metric) %>%
  summarise(Mean = mean(Value, na.rm = TRUE), SD = sd(Value, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Metric = factor(Metric, levels = names(metric_labels), labels = metric_labels),
    Comparison = factor(Comparison, levels = generalization_levels),
    Method = factor(Method, levels = c("ML", "DL", "NetTCR", "ERGO", "Ensemble"))
  )
generalization_points <- generalization_points %>%
  mutate(
    Metric = factor(Metric, levels = names(metric_labels), labels = metric_labels),
    Comparison = factor(Comparison, levels = generalization_levels)
  )

write_csv(generalization_bars, file.path(DIR_GENERALIZATION, "generalization_summary.csv"))
write_csv(generalization_points, file.path(DIR_GENERALIZATION, "generalization_folds.csv"))

generalization_model_order <- c(
  unname(ml_model_abbrev), "ResNet", "NetTCR-2.0 (ours)", "NetTCR-2.0 (paper)", "ERGO-Auto", "ERGO-LSTM", "Ensemble"
)

# Split into one standalone (non-faceted) plot per comparison, rather than a single
# multi-panel faceted figure - each one is simpler to drop straight into a main figure.
generalization_comparison_slugs <- c(
  "Core to ImmRep" = "core_to_immrep",
  "Core+ImmRep (train/test)" = "core_immrep_split"
)

plot_generalization_metric <- function(metric_name, comparison_name, out_file) {
  label <- metric_labels[[metric_name]]
  bar_df <- generalization_bars %>% filter(Metric == label, Comparison == comparison_name) %>%
    mutate(Model = factor(Model, levels = generalization_model_order))
  point_df <- generalization_points %>% filter(Metric == label, Comparison == comparison_name) %>%
    mutate(Model = factor(Model, levels = generalization_model_order))

  p <- ggplot(bar_df, aes(x = Model, y = Mean, fill = Method)) +
    geom_col(width = 0.65) +
    geom_errorbar(
      aes(ymin = pmax(Mean - SD, 0), ymax = pmin(Mean + SD, 1)),
      width = 0.2, na.rm = TRUE
    ) +
    geom_jitter(
      data = point_df, aes(x = Model, y = Value),
      inherit.aes = FALSE, width = 0.12, height = 0, size = 1.5, alpha = 0.7, shape = 21,
      fill = "white", color = "black", stroke = 0.4
    ) +
    scale_fill_manual(values = method_colors) +
    scale_y_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
    labs(x = NULL, y = label) +
    theme_minimal(base_size = 13) +
    theme(
      text = element_text(color = "black", size = 13, family = "Arial"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )

  ggsave(
    file.path(DIR_GENERALIZATION, out_file), p, width = 4, height = 3,
    device = svglite::svglite, user_fonts = ARIAL_USER_FONTS
  )
}

for (metric_name in c("Accuracy", "ROC_AUC", "PR_AUC", "PPV", "NPV")) {
  for (comparison_name in names(generalization_comparison_slugs)) {
    slug <- generalization_comparison_slugs[[comparison_name]]
    out_file <- paste0("generalization_", slug, "_", tolower(metric_name), ".svg")
    plot_generalization_metric(metric_name, comparison_name, out_file)
  }
}

# =============================================================================
# Sixth figure: core->ImmRep generalization again, but with immrep's heavy class
# imbalance (~598 positive of 3484 overall) removed - immrep is downsampled to ~50/50
# positive/negative many times over (see src/ScoreImmrepBalanced.py) and the same
# metrics recomputed on each balanced resample, to separate how much of the
# generalization gap seen in Figure 6 is simply class imbalance versus a genuine shift
# in peptide/MHC allele/TCR gene distribution between core and immrep. Adds ERGO-II (two
# published model variants, Auto/LSTM) alongside ML/DL/NetTCR - ERGO wasn't trained by
# us at all, so there's no "core" comparison for it beyond this one. NetTCR appears
# twice, explicitly labeled apart: "NetTCR-2.0 (ours)" (trained on our core set, tested
# on the structured immrep subset) vs "NetTCR-2.0 (paper)" (trained purely on NetTCR-2.0's
# own published training set, tested on the complete immrep set - no exposure to our data
# at all, much like ERGO). Every model here (even the single-fit ones, which otherwise
# have no inherent repeated-measurement structure) now gets a real error bar/point
# spread, since the resampling itself is the source of variability.
# =============================================================================
IMMREP_BALANCED_CSV <- "./data/Results/immrep_balanced_resamples.csv"

immrep_balanced_points <- read_csv(IMMREP_BALANCED_CSV, show_col_types = FALSE) %>%
  mutate(Model = if_else(Method == "ML", unname(ml_model_abbrev[Model]), Model)) %>%
  pivot_longer(cols = c(Accuracy, ROC_AUC, PR_AUC, PPV, NPV), names_to = "Metric", values_to = "Value")

immrep_balanced_bars <- immrep_balanced_points %>%
  group_by(Method, Model, Metric) %>%
  summarise(Mean = mean(Value, na.rm = TRUE), SD = sd(Value, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Metric = factor(Metric, levels = names(metric_labels), labels = metric_labels),
    Method = factor(Method, levels = c("ML", "DL", "NetTCR", "ERGO"))
  )
immrep_balanced_points <- immrep_balanced_points %>%
  mutate(Metric = factor(Metric, levels = names(metric_labels), labels = metric_labels))

write_csv(immrep_balanced_bars, file.path(DIR_IMMREP_BALANCED, "immrep_balanced_summary.csv"))
write_csv(immrep_balanced_points, file.path(DIR_IMMREP_BALANCED, "immrep_balanced_resamples.csv"))

immrep_balanced_model_order <- c(
  unname(ml_model_abbrev), "ResNet", "NetTCR-2.0 (ours)", "NetTCR-2.0 (paper)", "ERGO-Auto", "ERGO-LSTM"
)

plot_immrep_balanced_metric <- function(metric_name, out_file) {
  label <- metric_labels[[metric_name]]
  bar_df <- immrep_balanced_bars %>% filter(Metric == label) %>%
    mutate(Model = factor(Model, levels = immrep_balanced_model_order))
  point_df <- immrep_balanced_points %>% filter(Metric == label) %>%
    mutate(Model = factor(Model, levels = immrep_balanced_model_order))

  p <- ggplot(bar_df, aes(x = Model, y = Mean, fill = Method)) +
    geom_col(width = 0.65) +
    geom_errorbar(
      aes(ymin = pmax(Mean - SD, 0), ymax = pmin(Mean + SD, 1)),
      width = 0.2, na.rm = TRUE
    ) +
    geom_jitter(
      data = point_df, aes(x = Model, y = Value),
      inherit.aes = FALSE, width = 0.12, height = 0, size = 0.9, alpha = 0.4, shape = 21,
      fill = "white", color = "black", stroke = 0.3
    ) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
    scale_fill_manual(values = method_colors) +
    scale_y_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
    labs(x = NULL, y = label) +
    theme_minimal(base_size = 13) +
    theme(
      text = element_text(color = "black", size = 13, family = "Arial"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )

  ggsave(
    file.path(DIR_IMMREP_BALANCED, out_file), p, width = 4, height = 3,
    device = svglite::svglite, user_fonts = ARIAL_USER_FONTS
  )
}

plot_immrep_balanced_metric("Accuracy", "immrep_balanced_accuracy.svg")
plot_immrep_balanced_metric("ROC_AUC", "immrep_balanced_roc_auc.svg")
plot_immrep_balanced_metric("PR_AUC", "immrep_balanced_pr_auc.svg")
plot_immrep_balanced_metric("PPV", "immrep_balanced_ppv.svg")
plot_immrep_balanced_metric("NPV", "immrep_balanced_npv.svg")

# =============================================================================
# Seventh figure: per-sample outcomes on the Core+ImmRep FullRun split - the one
# comparison where ML/DL/NetTCR were all scored against the same fixed split (data/DL/
# FullRun/split.csv), so sample-level predictions actually line up across pipelines (see
# src/BuildFullRunSampleOutcomes.py). Each pipeline filtered that shared split down
# differently (ML drops misfolded structures, NetTCR drops non-class-I/too-long
# sequences), so this is restricted to the 161 samples every one of them produced a
# prediction for. One tile per (sample, model): TP/TN/FP/FN against the true label.
# "Ensemble" is a 9th column - the soft-voting mean of every model's predicted
# probability, thresholded at 0.5 - scored identically to the rest in the companion
# summary csv. Rows are faceted by source and, within each facet, sorted by how many of
# the 8 individual models (Ensemble excluded, so it doesn't bias its own sort) got that
# sample right, so easy/hard samples cluster together instead of sitting in file order.
#
# Not built with a dedicated heatmap-annotation library, so the row annotation
# (similarity-to-training, left - one column per attribute) is its own small ggplot,
# sharing the main heatmap's y factor levels and facet_grid(source~., space="free_y")
# structure, then assembled alongside it with patchwork, which aligns panel heights
# across the two since neither plot facets in more than one direction.
# =============================================================================
SAMPLE_OUTCOMES_CSV <- "./data/Results/fullrun_sample_outcomes.csv"

outcomes <- read_csv(SAMPLE_OUTCOMES_CSV, show_col_types = FALSE)

sample_order <- outcomes %>%
  filter(Model != "Ensemble") %>%
  group_by(name, source) %>%
  summarise(FracCorrect = mean(Outcome %in% c("TP", "TN")), .groups = "drop") %>%
  arrange(source, desc(FracCorrect), name)

heatmap_model_order <- c(unname(ml_model_abbrev), "ResNet", "NetTCR-2.0 (ours)", "Ensemble")

outcomes <- outcomes %>%
  mutate(
    name = factor(name, levels = sample_order$name),
    Model = if_else(Model %in% names(ml_model_abbrev), unname(ml_model_abbrev[Model]), Model),
    Model = factor(Model, levels = heatmap_model_order),
    Outcome = factor(Outcome, levels = c("TP", "TN", "FP", "FN"))
  )

write_csv(outcomes, file.path(DIR_SAMPLE_HEATMAP, "fullrun_sample_outcomes.csv"))

outcome_colors <- c(TP = "#1a9850", TN = "#91cf60", FP = "#ffd500", FN = "#d73027")

p_heatmap <- ggplot(outcomes, aes(x = Model, y = name, fill = Outcome)) +
  geom_tile(color = "white", linewidth = 0.1) +
  facet_grid(source ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = outcome_colors) +
  labs(x = NULL, y = NULL, fill = "Outcome") +
  theme_minimal(base_size = 13) +
  theme(
    text = element_text(color = "black", size = 13, family = "Arial"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    strip.text.y = element_blank(),
    legend.position = "top"
  )

# --- Row annotation: per-attribute normalized similarity to training (count of
# training samples sharing that exact Allele/Peptide/TRBV value, divided by that
# attribute's most common training value's count - see
# BuildFullRunSampleOutcomes.py's compute_similarity) ---
similarity_attrs <- c("Allele", "Peptide", "TRBV")
sample_annot <- outcomes %>%
  distinct(name, source, across(all_of(similarity_attrs))) %>%
  pivot_longer(cols = all_of(similarity_attrs), names_to = "Attribute", values_to = "Similarity") %>%
  mutate(Attribute = factor(Attribute, levels = similarity_attrs))

p_row_annot <- ggplot(sample_annot, aes(x = Attribute, y = name, fill = Similarity)) +
  geom_tile(color = "white", linewidth = 0.1) +
  facet_grid(source ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) +
  labs(x = NULL, y = NULL, fill = "Similarity\nto train") +
  theme_minimal(base_size = 13) +
  theme(
    text = element_text(color = "black", size = 13, family = "Arial"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    strip.placement = "outside",
    legend.position = "top"
  )

p_annotated_heatmap <- (p_row_annot + p_heatmap) +
  plot_layout(ncol = 2, widths = c(2.5, 5), guides = "collect") &
  theme(legend.position = "top")

ggsave(
  file.path(DIR_SAMPLE_HEATMAP, "fullrun_sample_heatmap.svg"), p_annotated_heatmap, width = 5, height = 9.5,
  device = svglite::svglite, user_fonts = ARIAL_USER_FONTS
)

cat("Saved plots and tables to", OUT_DIR, "\n")

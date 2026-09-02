# Central configuration for the complete R workflow.
# All formal differential calls are sex-stratified two-group limma results.
# No age-by-sex interaction term is fitted anywhere in this project.

cfg <- list(
  project = list(
    name = "Fuda four-layer age comparison: 28-panel reproducibility package",
    release = "2026-08-28",
    seed = 20260828L,
    strict = TRUE
  ),
  paths = list(
    protein_pairwise_dir = "input/protein_pairwise",
    metabolite_files = c(
      Cardiac_Metabolome = "input/metabolite/heart_data_unfiltered.xlsx",
      Serum_Metabolome = "input/metabolite/serum_data_unfiltered.xlsx",
      Urine_Metabolome = "input/metabolite/urine_data_unfiltered.xlsx"
    ),
    output_dir = "outputs"
  ),
  analysis = list(
    model = "sex-stratified unpaired two-group limma",
    allow_interaction = FALSE,
    paired = FALSE,
    ebayes = list(trend = FALSE, robust = FALSE, proportion = 0.01),
    minimum_non_missing_per_group = 2L,
    formal = list(
      p_column = "adj.P.Val",
      p_cutoff = 0.05,
      minimum_abs_log2FC = log2(1.5)
    )
  ),
  preprocessing = list(
    protein = list(
      input_scale = "log2",
      additional_normalization = "none",
      additional_imputation = "none",
      additional_transform = "none"
    ),
    metabolite = list(
      input_scale = "linear_abundance",
      zero_as_missing = TRUE,
      feature_filter = "at least two positive observations in each compared group",
      normalization = "none; preserve supplied quantitative scale",
      imputation = "none",
      log2_transform = TRUE,
      pseudocount = "none"
    )
  ),
  comparisons = data.frame(
    order = 1:6,
    source_sheet = c("FL-FM", "FM-FH", "FL-FH", "ML-MM", "MM-MH", "ML-MH"),
    comparison_id = c("Female_E_vs_M", "Female_M_vs_L", "Female_E_vs_L",
                      "Male_E_vs_M", "Male_M_vs_L", "Male_E_vs_L"),
    sex = c("Female", "Female", "Female", "Male", "Male", "Male"),
    first_code = c("FL", "FM", "FL", "ML", "MM", "ML"),
    second_code = c("FM", "FH", "FH", "MM", "MH", "MH"),
    first_label = c("E", "M", "E", "E", "M", "E"),
    second_label = c("M", "L", "L", "M", "L", "L"),
    stringsAsFactors = FALSE
  ),
  layers = data.frame(
    order = 1:4,
    layer_id = c("Cardiac_Proteome", "Cardiac_Metabolome",
                 "Serum_Metabolome", "Urine_Metabolome"),
    display_name = c("Cardiac proteome", "Cardiac metabolome",
                     "Serum metabolome", "Urine metabolome"),
    feature_id = c("Accession", "feature_id", "feature_id", "feature_id"),
    label_column = c("PG.Genes", "name", "name", "name"),
    stringsAsFactors = FALSE
  ),
  figures = list(
    formats = c("png", "pdf", "svg"),
    width = 7.2,
    height = 5.8,
    dpi = 600,
    volcano_label_n_each_direction = 8L,
    colors = c(Up = "#C23B33", Down = "#2D6FA3", Not_significant = "#B7B7B7")
  ),
  enrichment = list(
    enabled = TRUE,
    ontology = c("BP", "CC", "MF"),
    pvalue_cutoff_for_export = 1,
    fdr_cutoff_for_formal = 0.05,
    minimum_term_size = 3L,
    maximum_term_size = 500L,
    use_kegg_online_if_cache_missing = TRUE,
    top_n_plot = 15L
  ),
  network = list(
    enabled = TRUE,
    protein_symbols = c("Sirt4", "U2af2"),
    correlation_method = "pearson",
    minimum_complete_pairs = 5L,
    minimum_abs_correlation = 0.60,
    fdr_cutoff = 0.05
  ),
  expected = list(
    number_of_layers = 4L,
    number_of_comparisons_per_layer = 6L,
    number_of_differential_panels = 24L,
    number_of_upset_panels = 4L,
    number_of_main_panels = 28L
  )
)

stopifnot(identical(cfg$analysis$allow_interaction, FALSE))

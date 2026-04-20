#!/usr/bin/env Rscript
# Compute V-measure concordance between subregion-based and full-length genome clusters.
# Compares ClusterPicker and HIV-TRACE results for each subregion against the
# full-length genome reference (ClusterPicker GD=3.0%, BS=90%).
#
# Usage: Rscript compare_cluster_full.R regions.txt

suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(clevr)
  library(readr)
})

# ─── CONFIGURATION ────────────────────────────────────────────────────────────

# Update these paths before running
CP_BASE    <- "path/to/output/ClusterPicker"
TRACE_BASE <- "path/to/output/HIV_trace"
OUTPUT_DIR <- "path/to/output/V_measure"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

GD_THRESHOLDS    <- seq(from = 0.005, to = 0.045, by = 0.005)
BS_TAGS          <- c("bs70", "bs90", "bs99")
BS_NUMBER_FOLDER <- "bs_1000"
CP_LOG_FILE      <- "input_clusterPicks_log.txt"

# Full-length reference: ClusterPicker GD=3.0%, BS=90%
REF_REGION <- "full_length_584"
REF_GD     <- 0.030
REF_GD_TAG <- sprintf("gd%03d", as.integer(round(REF_GD * 1000)))
REF_BS_TAG <- "bs90"

FULL_CP_REF_PATH <- file.path(
  CP_BASE, paste0(REF_REGION, "_CP"), BS_NUMBER_FOLDER,
  paste0(REF_BS_TAG, "_", REF_GD_TAG), CP_LOG_FILE
)

FULL_TRACE_REF_PATH <- file.path(
  TRACE_BASE, REF_REGION,
  paste0(REF_REGION, "_hivtrace_GD", sprintf("%03d", as.integer(round(REF_GD * 1000))), ".json")
)

# ─── FUNCTIONS ────────────────────────────────────────────────────────────────

parse_cp <- function(log_path) {
  if (!file.exists(log_path)) {
    message("ClusterPicker log missing: ", log_path)
    return(NULL)
  }
  lines <- readLines(log_path, warn = FALSE)
  rows  <- grep("^\\d+\\t", lines, value = TRUE)
  if (!length(rows)) {
    message("No data rows in ClusterPicker log: ", log_path)
    return(NULL)
  }
  bind_rows(lapply(rows, function(r) {
    p <- str_split(r, "\\t", n = 6)[[1]]
    tibble(
      SequenceID = str_remove_all(p[4], "\\[|\\]|\\s") |>
        str_split(",") |> unlist() |> str_replace("^Seq_", ""),
      CP_cluster = as.integer(p[1])
    )
  }))
}

parse_trace <- function(json_path) {
  if (!file.exists(json_path)) {
    message("HIV-TRACE JSON missing: ", json_path)
    return(NULL)
  }
  j <- tryCatch(fromJSON(json_path, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(j) || is.null(j$trace_results$Nodes)) {
    message("Invalid HIV-TRACE JSON: ", json_path)
    return(NULL)
  }
  cluster_values <- j$trace_results$Nodes$cluster$values
  cluster_ids    <- j$trace_results$Nodes$id
  if (is.null(cluster_values) || length(cluster_values) == 0)
    cluster_values <- rep(0L, length(cluster_ids))
  tibble(
    SequenceID    = unlist(cluster_ids) |> as.character() |> str_replace("^Seq_", ""),
    TRACE_cluster = unlist(cluster_values) |> as.integer()
  )
}

calc_vm <- function(pred_df, ref_df, pred_col, ref_col) {
  merged <- inner_join(pred_df, ref_df, by = "SequenceID")
  if (nrow(merged) == 0) return(NA_real_)
  if (pred_col == ref_col) {
    pred_col <- paste0(pred_col, ".x")
    ref_col  <- paste0(ref_col,  ".y")
  }
  clevr::v_measure(merged[[pred_col]], merged[[ref_col]])
}

# ─── MAIN ─────────────────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript compare_cluster_full.R regions.txt")

regions <- readLines(args[1], warn = FALSE)
regions <- regions[nzchar(trimws(regions)) & !startsWith(trimws(regions), "#")]

message("Loading full-length reference clusters...")

ref_cp    <- parse_cp(FULL_CP_REF_PATH)
ref_trace <- parse_trace(FULL_TRACE_REF_PATH)

if (is.null(ref_cp))    stop("Full-length CP reference not found: ",       FULL_CP_REF_PATH)
if (is.null(ref_trace)) stop("Full-length HIV-TRACE reference not found: ", FULL_TRACE_REF_PATH)

message(sprintf("CP reference:    %d sequences", nrow(ref_cp)))
message(sprintf("TRACE reference: %d sequences", nrow(ref_trace)))

results <- list()

for (region in regions) {
  region <- trimws(region)
  message(sprintf("\n=== %s ===", region))

  for (gd in GD_THRESHOLDS) {
    gd_tag <- sprintf("gd%03d", as.integer(round(gd * 1000)))

    # ClusterPicker
    for (bs_tag in BS_TAGS) {
      log_path <- file.path(
        CP_BASE, paste0(region, "_CP"), BS_NUMBER_FOLDER,
        paste0(bs_tag, "_", gd_tag), CP_LOG_FILE
      )
      cp_df <- parse_cp(log_path)
      if (!is.null(cp_df)) {
        results[[length(results) + 1]] <- tibble(
          Region               = region,
          Method               = "ClusterPicker",
          GD_percent           = gd * 100,
          Bootstrap            = bs_tag,
          Reference_Method     = "ClusterPicker",
          Reference_GD         = REF_GD * 100,
          Reference_Bootstrap  = REF_BS_TAG,
          V_measure            = round(calc_vm(cp_df, ref_cp, "CP_cluster", "CP_cluster"), 6),
          n_predicted_clusters = length(unique(cp_df$CP_cluster)),
          n_reference_clusters = length(unique(ref_cp$CP_cluster)),
          n_predicted_seqs     = nrow(cp_df),
          n_reference_seqs     = nrow(ref_cp),
          n_shared_seqs        = nrow(inner_join(cp_df, ref_cp, by = "SequenceID"))
        )
      }
    }

    # HIV-TRACE
    json_path <- file.path(
      TRACE_BASE, region,
      paste0(region, "_hivtrace_GD", sprintf("%03d", as.integer(round(gd * 1000))), ".json")
    )
    trace_df <- parse_trace(json_path)
    if (!is.null(trace_df)) {
      results[[length(results) + 1]] <- tibble(
        Region               = region,
        Method               = "HIV-TRACE",
        GD_percent           = gd * 100,
        Bootstrap            = NA_character_,
        Reference_Method     = "HIV-TRACE",
        Reference_GD         = REF_GD * 100,
        Reference_Bootstrap  = NA_character_,
        V_measure            = round(calc_vm(trace_df, ref_trace, "TRACE_cluster", "TRACE_cluster"), 6),
        n_predicted_clusters = length(unique(trace_df$TRACE_cluster)),
        n_reference_clusters = length(unique(ref_trace$TRACE_cluster)),
        n_predicted_seqs     = nrow(trace_df),
        n_reference_seqs     = nrow(ref_trace),
        n_shared_seqs        = nrow(inner_join(trace_df, ref_trace, by = "SequenceID"))
      )
    }
  }
}

# ─── OUTPUT ───────────────────────────────────────────────────────────────────

if (length(results) > 0) {
  final_df <- bind_rows(results) |>
    arrange(Region, GD_percent, Method, Bootstrap)
  out_path <- file.path(OUTPUT_DIR, "vmeasure_results.tsv")
  write_tsv(final_df, out_path)
  message(sprintf("\nWrote results to %s", out_path))
} else {
  message("No V-measure results generated.")
}

message("\n--- Done ---")

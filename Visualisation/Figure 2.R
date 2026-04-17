# =========================================================================================
# Script Name: Figure 2 Generation for Phosphorylation Dynamics
# Description: This script processes simulation output data across different synaptic 
#              plasticity models (Receptor, Slot, Exocytosis) under varying stimulation 
#              protocols. It generates a comprehensive four-panel figure (Figure 2a - 2d) 
#              to visualize the time-series dynamics of various receptor phosphorylation 
#              states and CaMKII subunit metrics.
# =========================================================================================

# Set the working directory to the local GitHub repository clone containing the data folders
setwd("E://Kappa//ICA1") 

# Load necessary packages for data wrangling, string manipulation, functional programming, and visualization
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)
library(purrr)
library(patchwork)

# Define the relative path to the directory containing the simulation CSV results
data_dir <- "./Result"


# Define the target variable for the first analytical panel (Figure 2a) focusing on S845 phosphorylation
y_col <- "S845_phosphorylation_level"

# Custom CSV Reader for Simulation Outputs
# The output files from the simulation engine often contain metadata in the first two rows.
# This function safely bypasses these non-data rows, extracts the correct column headers 
# from the third row, and drops any empty or corrupted columns to return a clean data frame.
read_special_csv <- function(file_path) {
  
  raw_df <- read.csv(
    file_path,
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (nrow(raw_df) < 3) {
    stop(paste("File has fewer than 3 rows:", basename(file_path)))
  }
  
  new_colnames <- as.character(raw_df[3, ])
  dat <- raw_df[-c(1, 2, 3), , drop = FALSE]
  colnames(dat) <- new_colnames
  
  dat <- dat[, colnames(dat) != "" & !is.na(colnames(dat)), drop = FALSE]
  
  return(dat)
}

# Scan the target directory and compile a sorted list of all CSV data files
# Sorting by basename ensures consistent and chronological processing of experimental groups
csv_files <- list.files(
  path = data_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

csv_files <- csv_files[order(basename(csv_files))]

# File Metadata Parser
# This function leverages regular expressions to automatically extract critical experimental 
# background from the standardized filenames. It identifies the two-digit group ID representing 
# the stimulation protocol, the underlying mechanistic model type (Receptor, Slot, or Exocytosis), 
# and a cleaned condition string for downstream grouping.
parse_file_info <- function(file_path) {
  
  fname <- basename(file_path)
  
  group_id <- str_extract(fname, "^[0-9]{2}")
  
  model <- case_when(
    str_detect(fname, regex("Receptor", ignore_case = TRUE)) ~ "Receptor",
    str_detect(fname, regex("Slot", ignore_case = TRUE)) ~ "Slot",
    str_detect(fname, regex("Exocytosis", ignore_case = TRUE)) ~ "Exocytosis",
    TRUE ~ NA_character_
  )
  
  condition_name <- fname %>%
    str_remove("^[0-9]{2}_") %>%
    str_remove("_(Receptor|Slot|Exocytosis)_result\\.csv$") %>%
    str_remove("\\.csv$")
  
  tibble(
    file_path = file_path,
    file_name = fname,
    group_id = group_id,
    model = model,
    condition_name = condition_name
  )
}

# Apply the parser across all discovered files to build a comprehensive metadata registry
file_info <- map_dfr(csv_files, parse_file_info)


# Data Extraction and Transformation Pipeline for S845
# Iterate over the metadata registry to load, validate, and merge all time-series data.
# A derived variable representing the overall phosphorylation level of S845 is calculated 
# by taking the ratio of phosphorylated states (single and double) to the total pool 
# of available S845 sites (unphosphorylated plus phosphorylated).
all_data <- map2_dfr(file_info$file_path, seq_len(nrow(file_info)), function(fp, i) {
  
  info <- file_info[i, ]
  
  # Read file using custom function
  dat <- read_special_csv(fp)
  
  # Check required columns exist to perform the specific calculations
  required_cols <- c("[T]", "S845_P", "S845_PP", "S845_u")
  missing_cols <- setdiff(required_cols, colnames(dat))
  
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Missing required column(s) in file ", info$file_name, ": ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  # Extract time and calculate the derived S845 phosphorylation ratio
  out <- dat %>%
    transmute(
      Time = suppressWarnings(as.numeric(.data[["[T]"]])),
      S845_P_num = suppressWarnings(as.numeric(.data[["S845_P"]])),
      S845_PP_num = suppressWarnings(as.numeric(.data[["S845_PP"]])),
      S845_u_num = suppressWarnings(as.numeric(.data[["S845_u"]]))
    ) %>%
    mutate(
      Value = (S845_P_num + S845_PP_num) / (S845_P_num + S845_PP_num + S845_u_num),
      group_id = info$group_id,
      model = info$model
    ) %>%
    select(Time, Value, group_id, model)
  
  return(out)
})

# Scrub the dataset of any missing or infinite values generated during mathematical transformations
all_data <- all_data %>%
  filter(!is.na(Time), !is.na(Value), is.finite(Value))

# Enforce a specific factor level ordering for the models to ensure consistent color mapping
all_data$model <- factor(
  all_data$model,
  levels = c("Receptor", "Slot", "Exocytosis")
)

# Standardized Subplot Generator
# Given a specific group ID, this function isolates the relevant data and generates a line plot 
# depicting the target variable's temporal dynamics, applying a consistent visual theme.
plot_one_group <- function(gid) {
  
  df <- all_data %>% filter(group_id == gid)
  
  ggplot(df, aes(x = Time, y = Value, color = model)) +
    geom_line(linewidth = 0.4) +
    labs(
      title = title_map[gid],
      x = "Time (s)",
      y = "S845 phosphorylation level",
      color = "Model"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 13),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14)
    )
}

# Generate independent plots for each of the five stimulation protocols
p1 <- plot_one_group("01")
p2 <- plot_one_group("02")
p3 <- plot_one_group("03")
p4 <- plot_one_group("04")
p5 <- plot_one_group("05")

# Assemble the five subplots into a 1x5 horizontal layout using the patchwork package, 
# extracting the legend to a shared position to optimize visual space
p_combined <- (p1 | p2 | p3 | p4 | p5) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# Export the composite Figure 2a as a high-resolution PDF
ggsave(
  filename = file.path(data_dir, paste0("Figure2a_All_groups_", y_col, "_1x5.pdf")),
  plot = p_combined,
  width = 25,
  height = 3
)


# Initialize data processing for the second analytical panel (Figure 2b) focusing on S831 phosphorylation dynamics
y_col <- "S831_phosphorylation_level"
all_data <- map2_dfr(file_info$file_path, seq_len(nrow(file_info)), function(fp, i) {
  
  info <- file_info[i, ]
  
  # Read file using custom function
  dat <- read_special_csv(fp)
  
  # Check required columns exist for S831 computations
  required_cols <- c("[T]", "S831_P", "S831_PP", "S831_u")
  missing_cols <- setdiff(required_cols, colnames(dat))
  
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Missing required column(s) in file ", info$file_name, ": ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  # Extract time and calculate the derived S831 phosphorylation ratio
  out <- dat %>%
    transmute(
      Time = suppressWarnings(as.numeric(.data[["[T]"]])),
      S831_P_num = suppressWarnings(as.numeric(.data[["S831_P"]])),
      S831_PP_num = suppressWarnings(as.numeric(.data[["S831_PP"]])),
      S831_u_num = suppressWarnings(as.numeric(.data[["S831_u"]]))
    ) %>%
    mutate(
      Value = (S831_P_num + S831_PP_num) / (S831_P_num + S831_PP_num + S831_u_num),
      group_id = info$group_id,
      model = info$model
    ) %>%
    select(Time, Value, group_id, model)
  
  return(out)
})

# Remove rows with missing or invalid values
all_data <- all_data %>%
  filter(!is.na(Time), !is.na(Value), is.finite(Value))

# Ensure consistent legend order
all_data$model <- factor(
  all_data$model,
  levels = c("Receptor", "Slot", "Exocytosis")
)

# Redefine the plotting function to update the Y-axis label for S831 metrics
plot_one_group <- function(gid) {
  
  df <- all_data %>% filter(group_id == gid)
  
  ggplot(df, aes(x = Time, y = Value, color = model)) +
    geom_line(linewidth = 0.4) +
    labs(
      title = title_map[gid],
      x = "Time (s)",
      y = "S831 phosphorylation level",
      color = "Model"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 13),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14)
    )
}
p1 <- plot_one_group("01")
p2 <- plot_one_group("02")
p3 <- plot_one_group("03")
p4 <- plot_one_group("04")
p5 <- plot_one_group("05")

p_combined <- (p1 | p2 | p3 | p4 | p5) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(
  filename = file.path(data_dir, paste0("Figure2b_All_groups_", y_col, "_1x5.pdf")),
  plot = p_combined,
  width = 25,
  height = 3
)


# Initialize data processing for the third analytical panel (Figure 2c) focusing on S818 phosphorylation dynamics
y_col <- "S818_phosphorylation_level"
all_data <- map2_dfr(file_info$file_path, seq_len(nrow(file_info)), function(fp, i) {
  
  info <- file_info[i, ]
  
  # Read file using custom function
  dat <- read_special_csv(fp)
  
  # Check required columns exist for S818 computations
  required_cols <- c("[T]", "S818_p", "S818_pp", "S818_u")
  missing_cols <- setdiff(required_cols, colnames(dat))
  
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Missing required column(s) in file ", info$file_name, ": ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  # Extract time and calculate the derived S818 phosphorylation ratio
  out <- dat %>%
    transmute(
      Time = suppressWarnings(as.numeric(.data[["[T]"]])),
      S818_P_num = suppressWarnings(as.numeric(.data[["S818_p"]])),
      S818_PP_num = suppressWarnings(as.numeric(.data[["S818_pp"]])),
      S818_u_num = suppressWarnings(as.numeric(.data[["S818_u"]]))
    ) %>%
    mutate(
      Value = (S818_P_num + S818_PP_num) / (S818_P_num + S818_PP_num + S818_u_num),
      group_id = info$group_id,
      model = info$model
    ) %>%
    select(Time, Value, group_id, model)
  
  return(out)
})

# Remove rows with missing or invalid values
all_data <- all_data %>%
  filter(!is.na(Time), !is.na(Value), is.finite(Value))

# Ensure consistent legend order
all_data$model <- factor(
  all_data$model,
  levels = c("Receptor", "Slot", "Exocytosis")
)

# Redefine the plotting function to update the Y-axis label for S818 metrics
plot_one_group <- function(gid) {
  
  df <- all_data %>% filter(group_id == gid)
  
  ggplot(df, aes(x = Time, y = Value, color = model)) +
    geom_line(linewidth = 0.4) +
    labs(
      title = title_map[gid],
      x = "Time (s)",
      y = "S818 phosphorylation level",
      color = "Model"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 13),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14)
    )
}
p1 <- plot_one_group("01")
p2 <- plot_one_group("02")
p3 <- plot_one_group("03")
p4 <- plot_one_group("04")
p5 <- plot_one_group("05")

p_combined <- (p1 | p2 | p3 | p4 | p5) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(
  filename = file.path(data_dir, paste0("Figure2c_All_groups_", y_col, "_1x5.pdf")),
  plot = p_combined,
  width = 25,
  height = 3
)


# Initialize data processing for the fourth analytical panel (Figure 2d) examining Phosphorylated CaMKII subunits
y_col <- "Phos_CK_subunits"

all_data <- map2_dfr(file_info$file_path, seq_len(nrow(file_info)), function(fp, i) {
  
  info <- file_info[i, ]
  dat <- read_special_csv(fp)
  
  # Check required columns
  if (!"[T]" %in% colnames(dat)) {
    stop(paste0("Missing [T] column in file: ", info$file_name))
  }
  
  if (!y_col %in% colnames(dat)) {
    stop(paste0("Missing y column (", y_col, ") in file: ", info$file_name))
  }
  
  # Direct extraction of values without ratio calculations, unlike previous phosphorylation steps
  out <- dat %>%
    transmute(
      Time = suppressWarnings(as.numeric(.data[["[T]"]])),
      Value = suppressWarnings(as.numeric(.data[[y_col]]))
    ) %>%
    mutate(
      group_id = info$group_id,
      model = info$model
    )
  
  return(out)
})

# Clean NA
all_data <- all_data %>%
  filter(!is.na(Time), !is.na(Value))

# Reset factor order
all_data$model <- factor(
  all_data$model,
  levels = c("Receptor", "Slot", "Exocytosis")
)

# Return the line width back to standard 0.6 and update Y-axis labels for phosphorylated CaMKII subunits
plot_one_group <- function(gid) {
  
  df <- all_data %>% filter(group_id == gid)
  
  ggplot(df, aes(x = Time, y = Value, color = model)) +
    geom_line(linewidth = 0.6) +
    labs(
      title = title_map[gid],
      x = "Time (s)",
      y = "Phos CaMKII subunits",
      color = "Model"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 13),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14)
    )
}

p1 <- plot_one_group("01")
p2 <- plot_one_group("02")
p3 <- plot_one_group("03")
p4 <- plot_one_group("04")
p5 <- plot_one_group("05")

p_combined <- (p1 | p2 | p3 | p4 | p5) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(
  filename = file.path(data_dir, paste0("Figure2d_All_groups_", y_col, "_1x5.pdf")),
  plot = p_combined,
  width = 25,
  height = 3
)
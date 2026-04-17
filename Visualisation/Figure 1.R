# =========================================================================================
# Script Name: Figure 1 Generation for Synaptic Plasticity Models
# Description: This script processes simulation output data across different synaptic 
#              plasticity models (Receptor, Slot, Exocytosis) under varying stimulation 
#              protocols. It generates a four-panel figure (Figure 1a - 1d) visualizing 
#              the time-series dynamics of different synaptic components.
# =========================================================================================

# Set the working directory to the local GitHub repository clone containing the data
setwd("E://Kappa//ICA1") 

# Load necessary packages for data wrangling, functional programming, and visualization
library(ggplot2)    # Core plotting engine
library(dplyr)      # Data manipulation and pipeline operations
library(stringr)    # String parsing and regular expressions for file naming
library(readr)      # Efficient data reading
library(tidyr)      # Tidy data principles
library(purrr)      # Functional programming tools for mapping over files
library(patchwork)  # Combining multiple ggplot objects into comprehensive layouts

# Define the relative path to the directory containing the simulation CSV results
data_dir <- "./Result"


### Section: Data Ingestion and Helper Functions ###

# Define the target variable for the first plot panel (Figure 1a)
y_col <- "Synaptic_receptors"

#' Custom CSV Reader for Simulation Outputs
#' The output files from the simulation engine contain metadata in the first two rows.
#' The actual column names are located in the third row, followed by the data.
#' This function safely extracts the correct headers and filters out empty columns.

read_special_csv <- function(file_path) {
  
  # Read the raw file without treating the first row as headers
  raw_df <- read.csv(
    file_path,
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Validate file structural integrity
  if (nrow(raw_df) < 3) {
    stop(paste("File has fewer than 3 rows:", basename(file_path)))
  }
  
  # Extract the actual column names from the 3rd row and remove the metadata rows
  new_colnames <- as.character(raw_df[3, ])
  dat <- raw_df[-c(1, 2, 3), , drop = FALSE]
  colnames(dat) <- new_colnames
  
  # Drop any ghost columns that might have empty or NA headers
  dat <- dat[, colnames(dat) != "" & !is.na(colnames(dat)), drop = FALSE]
  
  return(dat)
}

# Scan the target directory and compile a sorted list of all CSV data files
csv_files <- list.files(
  path = data_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

# Ensure chronological/alphabetic processing by sorting the basenames
csv_files <- csv_files[order(basename(csv_files))]

#' File Metadata Parser
#' Extracts experimental conditions based on the standardized naming convention 
#' of the simulation output files.

parse_file_info <- function(file_path) {
  
  fname <- basename(file_path)
  
  # The first two digits of the filename represent the stimulation protocol group
  group_id <- str_extract(fname, "^[0-9]{2}")
  
  # Categorize the data based on the underlying mechanistic model used
  model <- case_when(
    str_detect(fname, regex("Receptor", ignore_case = TRUE)) ~ "Receptor",
    str_detect(fname, regex("Slot", ignore_case = TRUE)) ~ "Slot",
    str_detect(fname, regex("Exocytosis", ignore_case = TRUE)) ~ "Exocytosis",
    TRUE ~ NA_character_
  )
  
  # Extract a clean string representing the specific experimental condition
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

# Apply the parser across all found files to create a master metadata table
file_info <- map_dfr(csv_files, parse_file_info)


### Section: Data Processing and Plotting for Figure 1a (Synaptic Receptors) ###

# Iterate over the metadata table to load and merge all time-series data into a single dataframe
all_data <- map2_dfr(file_info$file_path, seq_len(nrow(file_info)), function(fp, i) {
  
  info <- file_info[i, ]
  dat <- read_special_csv(fp)
  
  # Guardrails to ensure necessary simulation tracking columns exist
  if (!"[T]" %in% colnames(dat)) {
    stop(paste0("Missing [T] column in file: ", info$file_name))
  }
  
  if (!y_col %in% colnames(dat)) {
    stop(paste0("Missing y column (", y_col, ") in file: ", info$file_name))
  }
  
  # Extract time and target value, enforcing numeric types and appending metadata
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

# Scrub the dataset of any coercion-induced NAs to prevent plotting artifacts
all_data <- all_data %>%
  filter(!is.na(Time), !is.na(Value))

# Enforce a logical ordering for the models to maintain consistent color mapping in plots
all_data$model <- factor(
  all_data$model,
  levels = c("Receptor", "Slot", "Exocytosis")
)

# Map numeric group IDs to their corresponding biological stimulation protocols
title_map <- c(
  "01" = "Single theta burst train",
  "02" = "Compressed theta burst",
  "03" = "Spaced theta burst",
  "04" = "wLTP burst",
  "05" = "High frequency stimulation"
)

#' Subplot Generation Function
#' Creates a standardized line plot for a specific stimulation protocol (group ID).

plot_one_group <- function(gid) {
  
  df <- all_data %>% filter(group_id == gid)
  
  ggplot(df, aes(x = Time, y = Value, color = model)) +
    geom_line(linewidth = 0.6) +
    labs(
      title = title_map[gid],
      x = "Time (s)",
      y = "Synaptic receptors",
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

# Generate individual panels for each of the 5 stimulation protocols
p1 <- plot_one_group("01")
p2 <- plot_one_group("02")
p3 <- plot_one_group("03")
p4 <- plot_one_group("04")
p5 <- plot_one_group("05")

# Compose the 1x5 layout using patchwork, extracting the legend to avoid redundancy
p_combined <- (p1 | p2 | p3 | p4 | p5) +
  plot_layout(guides = "collect") &   # << key: shared legend
  theme(
    legend.position = "right"
  )

# Export the composite Figure 1a to the designated output directory
ggsave(
  filename = file.path(data_dir, paste0("Figure1a_All_groups_", y_col, "_1x5.pdf")),
  plot = p_combined,
  width = 25,
  height = 3,
  dpi = 300
)


### Section: Data Processing and Plotting for Figure 1b (Slots) ###

# Update the target variable to analyze the total available slot numbers
y_col <- "Slots"

# Re-run the data extraction pipeline targeting the new variable
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

# Redefine the plotting function locally to update the Y-axis label appropriately
plot_one_group <- function(gid) {
  
  df <- all_data %>% filter(group_id == gid)
  
  ggplot(df, aes(x = Time, y = Value, color = model)) +
    geom_line(linewidth = 0.6) +
    labs(
      title = title_map[gid],
      x = "Time (s)",
      y = "Slot numbers",
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
  filename = file.path(data_dir, paste0("Figure1b_All_groups_", y_col, "_1x5.pdf")),
  plot = p_combined,
  width = 25,
  height = 3
)


### Section: Data Processing and Plotting for Figure 1c (Mobile AMPARs) ###

# Update the target variable to analyze the freely diffusing mobile AMPA receptors
y_col <- "AMPA_mobile"

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

# Note: The line width is specifically reduced to 0.25 here for mobile AMPAR visualization
plot_one_group <- function(gid) {
  
  df <- all_data %>% filter(group_id == gid)
  
  ggplot(df, aes(x = Time, y = Value, color = model)) +
    geom_line(linewidth = 0.25) +
    labs(
      title = title_map[gid],
      x = "Time (s)",
      y = "mobile AMPARs",
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
  filename = file.path(data_dir, paste0("Figure1c_All_groups_", y_col, "_1x5.pdf")),
  plot = p_combined,
  width = 25,
  height = 3
)


### Section: Data Processing and Plotting for Figure 1d (Endocytosed AMPARs) ###

# Update the target variable to analyze the internalized AMPA receptor pool
y_col <- "AMPA_endo"

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

# Return the line width back to standard 0.6 and update Y-axis labels for endocytosed receptors
plot_one_group <- function(gid) {
  
  df <- all_data %>% filter(group_id == gid)
  
  ggplot(df, aes(x = Time, y = Value, color = model)) +
    geom_line(linewidth = 0.6) +
    labs(
      title = title_map[gid],
      x = "Time (s)",
      y = "Endo AMPARs",
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
  filename = file.path(data_dir, paste0("Figure1d_All_groups_", y_col, "_1x5.pdf")),
  plot = p_combined,
  width = 25,
  height = 3
)
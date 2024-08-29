library(plyr)
library(dplyr)
library(tidyverse)
library(vegan)
library(zCompositions)
library(phyloseq)
library(metagMisc)
library(data.table)
library(htmlwidgets)
library(geosphere)
library(microViz)
library(reshape2)
library(tidyr)
library(readxl)
library("ggplot2")

# Define files
metadata_file <- "/Users/simplexdna/Desktop/RSDE BAC soil project/RSDE metabarcoding metadata.xlsx"
otu_table_file <- "/Users/simplexdna/Desktop/RSDE BAC soil project/apscale-new/rsde-bac-sediment_rerun_apscale_otu_table_filtered_microdecon-filtered_with_taxonomy.csv"

# Make folder for plots
plot_outdir <- "/Users/simplexdna/Desktop/RSDE BAC soil project/analysis/plots/fancy_plots"

#############################  Metadata processing
# Make sure metadata matches otu table
# otu_df_transposed <- otu_df %>%
#   t() %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "index")
# merged_df <- merge(otu_df_transposed, metadata, by.x = "index", by.y = "Sample ID", all = TRUE)

# Read in metadata
metadata <- as.data.frame(read_excel(metadata_file))

# Cut down metadata and process data types
rownames(metadata) <- metadata[, "Sample ID"]
metadata <- metadata[, -(which(names(metadata) == "Sample ID"))]
metadata <- subset(metadata, select = c("Latitude","Longitude", "Depth"))
metadata$Latitude <- as.numeric(metadata$Latitude)
metadata$Longitude <- as.numeric(metadata$Longitude)

###########################
# Calculate the distance between samples and add the distance as new metadata variables
# Find the most northern sample (highest latitude)
most_northern <- metadata[which.max(metadata$Latitude), ]
# Calculate distances from the most northern sample
metadata$distance_from_north <- apply(metadata, 1, function(row) {
  distHaversine(c(row['Longitude'], row['Latitude']), c(most_northern$Longitude, most_northern$Latitude))
})
metadata$distance_from_north_log <- log(metadata$distance_from_north)
# Find the sample with the smallest depth
smallest_depth <- metadata[which.min(metadata$Depth), ]
# Calculate depth differences from the smallest depth sample
metadata$depth_difference <- apply(metadata, 1, function(row) {
  abs(row['Depth'] - smallest_depth$Depth)
})
# Calculate the combined distance, incorporating both horizontal and vertical components
metadata$combined_distance <- apply(metadata, 1, function(row) {
  horizontal_distance <- distHaversine(c(row['Longitude'], row['Latitude']), c(most_northern$Longitude, most_northern$Latitude))
  depth_difference <- abs(row['Depth'] - smallest_depth$Depth)
  sqrt(horizontal_distance^2 + depth_difference^2)  # Euclidean distance with depth component
})

# Add depth zones
metadata <- metadata %>%
  mutate(
    "Depth category" = case_when(
      Depth > 1000 ~ "Bathypelagic (>1000m)",
      Depth >= 200 & Depth <= 1000 ~ "Mesopelagic (200-1000m)",
      Depth <= 200 ~ "Epipelagic (<200m)",
      TRUE ~ NA_character_  # default condition (optional)
    )
  )
metadata$'Depth category' <- factor(metadata$'Depth category', levels=c("Epipelagic (<200m)", "Mesopelagic (200-1000m)", "Bathypelagic (>1000m)"))

################################
# - Read in OTU table and create the output folder
otu_df <- read.csv(otu_table_file, row.names="ID")
dir.create(plot_outdir, showWarnings = FALSE)

# - Make a phyloseq object
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
otu_df <- otu_df %>%
  mutate_at(vars(ranks), ~replace(., . %in% c("Not available", "Not_available", "Taxonomy unreliable", "No match in database", "Unreliable taxonomy"), NA))
# OTU table: grab all columns starting with LF and turn them into a phyloseq OTU table object
OTU <- otu_table(otu_df[, grepl("^RSDE", names(otu_df))], taxa_are_rows = TRUE)
# Tax table: has to be a matrix for phyloseq to be happy
TAX <- tax_table(as.matrix(otu_df[, ranks]))
# Metadata: row names have to align with sample names in the OTU table
META <- sample_data(metadata)
# Make the phyloseq object
physeq = phyloseq(OTU, TAX, META)
# Threshold rel abundance > 0.00000007 (determined in general_data_exploration.R script)
thresh_relabun = 0.00000007
keepTaxa = (taxa_sums(physeq) / sum(taxa_sums(physeq))) > thresh_relabun
physeq = prune_taxa(keepTaxa, physeq)
# Threshold prevalence >2  (determined in general_data_exploration.R script)
tresh_prev = 2
physeq <- filter_taxa(physeq, function(x) sum(x > 0) > tresh_prev, prune = T)
physeq <- tax_fix(physeq, unknowns = c("Incertae_Sedis", "Incertae_Sedis class", "Gammaproteobacteria_Incertae_Sedis", "Alphaproteobacteria_Incertae_Sedis", "Unknown_Family"))

###############################
# Calculate alpha diversity
alpha_div <- estimate_richness(physeq)

# Merge dfs on row names (indices)
alpha_div_distance <- merge(alpha_div, metadata, by = "row.names")
row.names(alpha_div_distance) <- alpha_div_distance$Row.names
alpha_div_distance$Row.names <- NULL
alpha_div_distance$"Depth_category" <- alpha_div_distance$"Depth category"
#Drop most northern point
alpha_div_distance <- alpha_div_distance[!rownames(alpha_div_distance) %in% 'RSDES271', ]

# Create linear models for each category
models <- alpha_div_distance %>%
  group_by(Depth_category) %>%
  do(model = lm(Observed ~ distance_from_north, data = .))  # Create a linear model for each category
# Add p-values to the data frame
models <- models %>%
  mutate(
    r_squared = summary(model)$r.squared,
    p_value = summary(model)$coefficients[2, 4],  # p-value for the slope
    slope = coef(model)[2],  # Slope
    intercept = coef(model)[1]  # Intercept
  )
plot <- ggplot(alpha_div_distance, aes(x = distance_from_north, y = Observed, color = Depth_category)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, fullrange=TRUE, aes(color=Depth_category, fill=Depth_category)) +
  theme_minimal() +
  labs(
    x = "Distance from most northern sample",
    y = "OTU richness",
    color = "Depth Category"
  )
# Add annotations for each category
for (cat in unique(models$Depth_category)) {
  # Get linear model parameters
  slope <- models %>%
    filter(Depth_category == cat) %>%
    pull(slope)
  
  intercept <- models %>%
    filter(Depth_category == cat) %>%
    pull(intercept)
  
  r_squared <- models %>%
    filter(Depth_category == cat) %>%
    pull(r_squared)
  
  p_value <- models %>%
    filter(Depth_category == cat) %>%
    pull(p_value)
  
  if (cat == "Epipelagic (<200m)") {
    annotation_y = 68000
  } else if (cat == "Mesopelagic (200-1000m)") {
    annotation_y = 65000
  } else if (cat == "Bathypelagic (>1000m)") {
    annotation_y = 62000
  } 
  
  # Annotate the plot with p-values
  plot <- plot +
    annotate(
      "text",
      x = 0,  # Fixed position to the right of the y-axis
      y = annotation_y,  # Corresponds to the regression line
      label = sprintf("p: %.3f, m: %.5f", p_value, slope),
      color = scales::hue_pal()(3)[which(unique(models$Depth_category) == cat)],  # Consistent color
      hjust = 0  # Align text to the right side of the left y-axis
    )
}
# Display the final plot with annotations
plot

ggsave(units="px", width=4000, height=3000, dpi=600, file.path(plot_outdir, "distance_decay.png"), plot = plot)
ggsave(units="px", width=4000, height=3000, dpi=600, file.path(plot_outdir, "distance_decay.svg"), plot = plot)


###############################
# Generate distance matrix for all samples
distance_matrix <- physeq %>%
  tax_transform(rank = "phylum", trans = "clr") %>%
  dist_calc(dist = "euclidean") %>%
  dist_get() %>%
  as.matrix()
# Change format into one variable
dist_df <- melt(distance_matrix, varnames = c("RowSample", "ColSample"), value.name = "Dissimilarity")
# Drop rows for identical samples
dist_df_all <- dist_df[dist_df$RowSample != dist_df$ColSample, ]
# Only keep rows that are against the most northern sample
dist_df <- dist_df[dist_df$RowSample == 'RSDES271', ]
# Drop most northern sample
dist_df <- dist_df[dist_df$ColSample != 'RSDES271', ]

##############################
# Calculate distance and depth
dist_df <- dist_df %>%
  mutate(
    lat_RowSample = metadata[RowSample, "Latitude"],
    lon_RowSample = metadata[RowSample, "Longitude"],
    depth_RowSample = metadata[RowSample, "Depth"]
  )
dist_df <- dist_df %>%
  mutate(
    lat_ColSample = metadata[ColSample, "Latitude"],
    lon_ColSample = metadata[ColSample, "Longitude"],
    depth_ColSample = metadata[ColSample, "Depth"]
  )
# Calculate geographic distance and depth distance
dist_df <- dist_df %>%
  mutate(
    lat_lon_distance = distHaversine(cbind(lon_RowSample, lat_RowSample), cbind(lon_ColSample, lat_ColSample)),
    depth_difference = abs(depth_RowSample - depth_ColSample)  # Absolute difference in depths
  )
dist_df$lat_lon_distance_log <- log(dist_df$lat_lon_distance)
# Add depth zones
dist_df <- dist_df %>%
  mutate(
    "Depth_category" = case_when(
      depth_ColSample > 1000 ~ "Bathypelagic (>1000m)",
      depth_ColSample >= 200 & depth_ColSample <= 1000 ~ "Mesopelagic (200-1000m)",
      depth_ColSample <= 200 ~ "Epipelagic (<200m)",
      TRUE ~ NA_character_  # default condition (optional)
    )
  )
dist_df$'Depth_category' <- factor(dist_df$'Depth_category', levels=c("Epipelagic (<200m)", "Mesopelagic (200-1000m)", "Bathypelagic (>1000m)"))

########################
# Based on dissimilarity to most northern point
# Create linear models for each category
models <- dist_df %>%
  group_by(Depth_category) %>%
  do(model = lm(Dissimilarity ~ lat_lon_distance, data = .))  # Create a linear model for each category
# Add p-values to the data frame
models <- models %>%
  mutate(
    r_squared = summary(model)$r.squared,
    p_value = summary(model)$coefficients[2, 4],  # p-value for the slope
    slope = coef(model)[2],  # Slope
    intercept = coef(model)[1]  # Intercept
  )
plot <- ggplot(dist_df, aes(x = lat_lon_distance, y = Dissimilarity, color = Depth_category)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(
    x = "Distance from North",
    y = "Dissimilarity",
    color = "Depth Category"
  )
# Add annotations for each category
for (cat in unique(models$Depth_category)) {
  # Get linear model parameters
  p_value <- models %>%
    filter(Depth_category == cat) %>%
    pull(p_value)
  
  if (cat == "Epipelagic (<200m)") {
    annotation_y = 59
  } else if (cat == "Mesopelagic (200-1000m)") {
    annotation_y = 57
  } else if (cat == "Bathypelagic (>1000m)") {
    annotation_y = 55
  } 
  # Annotate the plot with R² and p-values
  plot <- plot +
    annotate(
      "text",
      x = 0,  # Fixed position to the right of the y-axis
      y = annotation_y,  # Corresponds to the regression line
      label = sprintf("p: %.3f", p_value),
      color = scales::hue_pal()(3)[which(unique(models$Depth_category) == cat)],  # Consistent color
      hjust = 0  # Align text to the right side of the left y-axis
    )
}
# Display the final plot with annotations
plot

# Plot without categories
diss_model <- lm(Dissimilarity ~ lat_lon_distance, data = dist_df)
p_value_all <- summary(diss_model)$coefficients[2, 4]
plot <- ggplot(dist_df, aes(x = lat_lon_distance, y = Dissimilarity)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(
    x = "Distance from North",
    y = "Dissimilarity",
  ) +
  annotate(
    "text",
    x = 0,  # Fixed position to the right of the y-axis
    y = 55,  # Corresponds to the regression line
    label = sprintf("p: %.3f", p_value_all),
    color = scales::hue_pal()(3)[which(unique(models$Depth_category) == cat)],  # Consistent color
    hjust = 0  # Align text to the right side of the left y-axis
  )

plot

###############################
# Based on distance between all points

# Calculate distance and depth
dist_df_all <- dist_df_all %>%
  mutate(
    lat_RowSample = metadata[RowSample, "Latitude"],
    lon_RowSample = metadata[RowSample, "Longitude"],
    depth_RowSample = metadata[RowSample, "Depth"]
  )
dist_df_all <- dist_df_all %>%
  mutate(
    lat_ColSample = metadata[ColSample, "Latitude"],
    lon_ColSample = metadata[ColSample, "Longitude"],
    depth_ColSample = metadata[ColSample, "Depth"]
  )
# Calculate geographic distance and depth distance
dist_df_all <- dist_df_all %>%
  mutate(
    lat_lon_distance = distHaversine(cbind(lon_RowSample, lat_RowSample), cbind(lon_ColSample, lat_ColSample)),
    depth_difference = abs(depth_RowSample - depth_ColSample)  # Absolute difference in depths
  )
dist_df_all$lat_lon_distance_log <- log(dist_df_all$lat_lon_distance)
########################
# Model and visualize
plot <- ggplot(dist_df_all, aes(x = lat_lon_distance, y = Dissimilarity)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(
    x = "Distance",
    y = "Dissimilarity"
  )
plot
diss_model <- lm(Dissimilarity ~ lat_lon_distance, data = dist_df_all)
summary(diss_model)$coefficients[2, 4]

## Distance
# Visualize the relationship between dissimilarity and distance
plot(dist_df_all$lat_lon_distance, dist_df_all$Dissimilarity, main = "Dissimilarity vs. Distance",
     xlab = "Distance", ylab = "Dissimilarity")
# Fit a linear model
linear_model <- lm(Dissimilarity ~ lat_lon_distance, data = dist_df_all)
summary(linear_model)  # Summary of the linear model
# Fit an exponential model (requires a transformation)
dist_df_all$log_dissimilarity <- log(dist_df_all$Dissimilarity)
exp_model <- lm(log_dissimilarity ~ lat_lon_distance, data = dist_df_all)
summary(exp_model)
# To visualize the exponential fit, use the following:
plot(dist_df_all$lat_lon_distance, dist_df_all$Dissimilarity, main = "Exponential Fit",
     xlab = "Distance", ylab = "Dissimilarity")
lines(dist_df_all$lat_lon_distance, exp(predict(exp_model)), col = "red")  # Red line for the exponential fit

## Depth
# Visualize the relationship between dissimilarity and depth difference
plot(dist_df_all$depth_difference, dist_df_all$Dissimilarity, main = "Dissimilarity vs. Depth difference",
     xlab = "Depth difference", ylab = "Dissimilarity")
# Fit a linear model
linear_model <- lm(Dissimilarity ~ depth_difference, data = dist_df_all)
summary(linear_model)  # Summary of the linear model
lines(dist_df_all$depth_difference, predict(linear_model), col = "red")  # Red line for the exponential fit
# Fit an exponential model (requires a transformation)
dist_df_all$log_dissimilarity <- log(dist_df_all$Dissimilarity)
exp_model <- lm(log_dissimilarity ~ depth_difference, data = dist_df_all)
summary(exp_model)
# To visualize the exponential fit, use the following:
plot(dist_df_all$depth_difference, dist_df_all$Dissimilarity, main = "Exponential Fit",
     xlab = "Depth difference", ylab = "Dissimilarity")
lines(dist_df_all$depth_difference, exp(predict(exp_model)), col = "red")  # Red line for the exponential fit

library(plotly)

fig <- plot_ly(dist_df_all, x = ~Dissimilarity, y = ~lat_lon_distance, z = ~depth_difference, marker = list(size = 1))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Dissimilarity'),
                                   yaxis = list(title = 'lat_lon_distance'),
                                   zaxis = list(title = 'depth_difference')))

fig
library(plyr)
library(dplyr)
library(tidyverse)
library(plotly)
library(readxl)
library(vegan)
library(zCompositions)
library(phyloseq)
library(metagMisc)
library(data.table)
library(leaflet)
library(viridis)
library(mapview)
library(webshot)
library(htmlwidgets)
library(iNEXT)

# Read in data
metadata_file <- "/Users/simplexdna/Desktop/RSDE BAC soil project/RSDE metabarcoding metadata.xlsx"
otu_df_file <- "/Users/simplexdna/Desktop/RSDE BAC soil project/apscale-new/rsde-bac-sediment_rerun_apscale_otu_table_filtered_microdecon-filtered_with_taxonomy.csv"
metadata <- as.data.frame(read_excel(metadata_file))
otu_df <- read.csv(otu_df_file, row.names="ID")

# Make folder for plots
plot_outdir <- "/Users/simplexdna/Desktop/RSDE BAC soil project/analysis/plots/fancy_plots/"
dir.create(plot_outdir, showWarnings = FALSE)

# Define rank for analysis (options: "domain", "phylum", "class", "order", "family", "genus", "species", "OTUs")
rank = "phylum"

#############################  Metadata processing and exploration
# Make sure metadata matches otu table
# otu_df_transposed <- otu_df %>%
#   t() %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "index")
# merged_df <- merge(otu_df_transposed, metadata, by.x = "index", by.y = "Sample ID", all = TRUE)

# Cut down metadata and process data types
rownames(metadata) <- metadata[, "Sample ID"]
metadata$Leg <- as.factor(metadata$Leg)
metadata$Phase <- as.factor(metadata$Phase)
metadata$Latitude <- as.numeric(metadata$Latitude)
metadata$Longitude <- as.numeric(metadata$Longitude)

# Add column with depth categories based on Larissa's recommendations
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
metadata$'Red Sea zone' <- factor(metadata$'Red Sea zone', levels=c("Gulf of Aqaba", "North Red Sea", "Central-north Red Sea", "Central-south Red Sea", "South Red Sea"))

# Show distribution of depth and depth categories
dist_depth_plot <- ggplot(na.omit(metadata), aes(x = Depth)) +
  geom_histogram(bins=60, fill = "blue", color = "black", alpha = 0.7) +
  scale_x_continuous(breaks = seq(0, max(metadata$Depth, na.rm = TRUE), by = 100)) +
  labs(title = "Depth Distribution", x = "Depth", y = "Frequency")
ggsave(file.path(plot_outdir, "distribution_depth.png"), plot = dist_depth_plot)
depth_cat_counts <- as.data.frame(table(metadata$"Depth category"))
depth_cat_counts$Var1 <- factor(depth_cat_counts$Var1, levels=c("Epipelagic (<200m)", "Mesopelagic (200-1000m)", "Bathypelagic (>1000m)"))
dist_depthcat_plot <- ggplot(depth_cat_counts, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "blue", color="black", alpha=0.7) +
  labs(title = "Distribution of depth categories", 
       x = "Depth category", 
       y = "Count")
ggsave(file.path(plot_outdir, "distribution_depth_categories.png"), plot = dist_depthcat_plot)

# Show sampling points
coords = na.omit(metadata[, c("Latitude", "Longitude")])
# Coloured by phases
cols_for_map_phases <- colorFactor(palette = viridis_pal(option = "D")(3), metadata$Phase)
map_phases <- leaflet(coords) %>%
  # Set the initial view (center and zoom level)
  setView(lng = mean(coords$"Longitude"), lat = mean(coords$"Latitude"), zoom = 5) %>%
  addTiles()  %>% # Add default OpenStreetMap tiles
  addCircleMarkers(radius=1,
                   fillOpacity = 0.8,
                   color = ~cols_for_map_phases(metadata$Phase),
                   popup = paste("Latitude: ", coords$Latitude, "<br>Longitude: ", coords$Longitude)) %>% # Add markers for each point
  addLegend('bottomright', pal = cols_for_map_phases, values = metadata$Phase,
            title = 'Phases',
            opacity = 1)
mapshot(map_phases, file = file.path(plot_outdir, "sampling_sites_phases.png"))
# Coloured by legs
cols_for_map_legs <- colorFactor(palette = viridis_pal(option = "D")(5), metadata$Leg)
map_legs <- leaflet(coords) %>%
  # Set the initial view (center and zoom level)
  setView(lng = mean(coords$"Longitude"), lat = mean(coords$"Latitude"), zoom = 5) %>%
  addTiles()  %>% # Add default OpenStreetMap tiles
  addCircleMarkers(radius=1,
                   fillOpacity = 0.8,
                   color = ~cols_for_map_legs(metadata$Leg),
                   popup = paste("Latitude: ", coords$Latitude, "<br>Longitude: ", coords$Longitude)) %>% # Add markers for each point
  addLegend('bottomright', pal = cols_for_map_legs, values = metadata$Leg,
            title = 'Legs',
            opacity = 1)
mapshot(map_legs, file = file.path(plot_outdir, "sampling_sites_legs.png"))
# Coloured by Red Sea Zone
cols_for_map_zone <- colorFactor(palette = c("orange", "blue", "darkgreen", "yellow", "#FF6A65", "#9A9B21", "#00B674", "#00A7F0", "#EB61ED"), metadata$'Zone & Discovery')
map_zone <- leaflet(coords) %>%
  # Set the initial view (center and zoom level)
  setView(lng = mean(coords$"Longitude"), lat = mean(coords$"Latitude"), zoom = 5) %>%
  addTiles() %>% 
  #addProviderTiles("Esri.WorldImagery") %>%
  addProviderTiles("USGS.USImageryTopo") %>%
  addCircleMarkers(radius=1,
               fillOpacity = 1,
               color = ~cols_for_map_zone(metadata$'Zone & Discovery'),
               weight = 6,  # Set the border width
               popup = paste("Latitude: ", coords$Latitude, "<br>Longitude: ", coords$Longitude)) %>% # Add markers for each point
  addLegend('topright', pal = cols_for_map_zone, values = metadata$'Zone & Discovery',
            title = 'Red Sea zone',
            opacity = 1)
mapshot(map_zone, file = file.path(plot_outdir, "map.png"), vwidth = 300, vheight = 550)
# Coloured by depth
cols_for_map_depth <- colorNumeric(palette = "viridis", metadata$Depth)
map_depth <- leaflet(coords) %>%
  # Set the initial view (center and zoom level)
  setView(lng = mean(coords$"Longitude"), lat = mean(coords$"Latitude"), zoom = 5) %>%
  addTiles()  %>% # Add default OpenStreetMap tiles
  addCircleMarkers(radius=1,
                   fillOpacity = 0.8,
                   color = ~cols_for_map_depth(metadata$Depth),
                   popup = paste("Latitude: ", coords$Latitude, "<br>Longitude: ", coords$Longitude)) %>% # Add markers for each point
  addLegend('bottomright', pal = cols_for_map_depth, values = metadata$Depth,
            title = 'Depth [m]',
            opacity = 1)
mapshot(map_depth, file = file.path(plot_outdir, "sampling_sites_depth.png"))
# Coloured by depth category
cols_for_map_depthcat <- colorFactor(palette = viridis_pal(option = "D")(4), metadata$"Depth category")
map_depthcat <- leaflet(coords) %>%
  # Set the initial view (center and zoom level)
  setView(lng = mean(coords$"Longitude"), lat = mean(coords$"Latitude"), zoom = 5) %>%
  addTiles()  %>% # Add default OpenStreetMap tiles
  addCircleMarkers(radius=1,
                   fillOpacity = 0.8,
                   color = ~cols_for_map_depthcat(metadata$"Depth category"),
                   popup = paste("Latitude: ", coords$Latitude, "<br>Longitude: ", coords$Longitude)) %>% # Add markers for each point
  addLegend('bottomright', pal = cols_for_map_depthcat, values = metadata$"Depth category",
            title = 'Depth category',
            opacity = 1)
mapshot(map_depthcat, file = file.path(plot_outdir, "sampling_sites_depth_categories.png"))

# Show Distribution of Legs and Phases
leg_cat_counts <- as.data.frame(table(metadata$"Leg"))
leg_cat_counts$Var1 <- factor(leg_cat_counts$Var1, levels = c("1", "2", "3", "4"))
dist_leg_plot <- ggplot(leg_cat_counts, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "blue", color="black", alpha=0.7) +
  labs(title = "Distribution of leg samples", 
       x = "Leg", 
       y = "Count")
ggsave(file.path(plot_outdir, "distribution_legs.png"), plot = dist_leg_plot)
phase_cat_counts <- as.data.frame(table(metadata$"Phase"))
phase_cat_counts$Var1 <- factor(phase_cat_counts$Var1, levels = c("1", "2", "DSC"))
dist_phase_plot <- ggplot(phase_cat_counts, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "blue", color="black", alpha=0.7) +
  labs(title = "Distribution of phase samples", 
       x = "Phase category", 
       y = "Count")
ggsave(file.path(plot_outdir, "distribution_phases.png"), plot = dist_phase_plot)

############## Make physeq object
# Process otu table
## Drop Sequence and domain info
otu_df <- subset(otu_df, select = -c(Seq))

## Standardize taxonomy for anything not identifiable
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
otu_df <- otu_df %>%
  mutate_at(vars(ranks), ~replace(., . %in% c("Not available", "Not_available", "Unreliable taxonomy"), NA))

# Turn tables into phyloseq object
OTU <- otu_table(otu_df[, grepl("^RSDE", names(otu_df))], taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(otu_df[, ranks]))
META <- sample_data(metadata)
physeq = phyloseq(OTU, TAX, META)

##########################      General data exploration (taken from https://evomics.org/wp-content/uploads/2016/01/phyloseq-Lab-01-Answers.html)
# Overview
microbiome::summarize_phyloseq(physeq)

# Determine the lowest level of taxonomic classification and plot it
# TAKES A LONG TIME, JUST RUN TO GENERATE PLOT ONCE
lowest_levels <- get_max_taxonomic_rank(physeq, return_rank_only = TRUE)
low_tax_plot <- ggplot(as.data.frame(table(lowest_levels)[ranks]), aes(x = lowest_levels, y = Freq)) +
  geom_bar(stat = "identity", fill = "blue", color="black", alpha=0.7) +
  labs(title = "Lowest OTU classification rank frequencies", 
       x = "Rank", 
       y = "Frequency")
ggsave(file.path(plot_outdir, "lowest_otu_classification_rank_frequencies.png"), plot = low_tax_plot)

# Determine seq depth distribution
sdt = data.table(as(sample_data(physeq), "data.frame"),
                 TotalReads = sample_sums(physeq), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) +
  geom_histogram(bins=30, fill = "blue", color="black", alpha=0.7) +
  labs(title = "Histogram of sequencing depth", 
       x = "Number of reads", 
       y = "Sample count")
ggsave(file.path(plot_outdir, "sequencing_depth.png"), plot = pSeqDepth)
pSeqDepth_categories <- pSeqDepth +
  facet_wrap(~Depth.category) +
  ggtitle("Histogram of sequencing depth by sampling depth category")
ggsave(file.path(plot_outdir, "sequencing_depth_by_depth.png"), plot = pSeqDepth_categories)

# Plot total OTU total counts
tdt_tot = data.table(tax_table(physeq),
                     TotalCounts = taxa_sums(physeq),
                     OTU = taxa_names(physeq))
otu_abundances_tot <- ggplot(tdt_tot, aes(TotalCounts)) + 
  geom_histogram(bins=100, fill = "blue", color = "black", alpha = 0.7) + 
  scale_y_continuous(trans='log10') +
  labs(title = "Histogram of total OTU abundances", 
       x = "Number of reads", 
       y = "OTU count")
ggsave(file.path(plot_outdir, "total_otu_abundances.png"), plot = otu_abundances_tot)

# Different format:
abundance_violin <- ggplot(tdt_tot, aes(x=1, y = TotalCounts)) +
  geom_violin(fill = "blue", color = "black", alpha = 0.7) +
  xlab("") +
  scale_y_log10() +
  ylab("Total OTU abundance - log x") +
  theme(axis.text.x = element_blank()) +
  ggtitle("Violin Plot of total OTU abundances")
ggsave(file.path(plot_outdir, "total_otu_abundances_violin.png"), plot = abundance_violin)
# --> what's up around 10k abundance?

# Plot relative OTU total counts
tdt_rel = data.table(tax_table(physeq),
                     RelCounts = taxa_sums(physeq)/sum(taxa_sums(physeq))*100,
                     OTU = taxa_names(physeq))
otu_abundances_rel <- ggplot(tdt_rel, aes(RelCounts)) + 
  geom_histogram(bins=100, fill = "blue", color = "black", alpha = 0.7) + 
  scale_y_continuous(trans='log10') +
  labs(title = "Histogram of relative OTU abundances", 
       x = "Percentage of reads", 
       y = "OTU count")
ggsave(file.path(plot_outdir, "relative_otu_abundances.png"), plot = otu_abundances_rel)

# Number of empty OTUs
cat("Number of empty OTUs:", tdt_tot[(TotalCounts == 0), .N])
# Number of singletons
cat("Number of singletons:", tdt_tot[(TotalCounts == 1), .N])

# Plot cumulative total OTU sum to determine total count filtering threshold
taxcumsum_tot = tdt_tot[, .N, by = TotalCounts]
setkey(taxcumsum_tot, TotalCounts)
taxcumsum_tot[, CumSum := cumsum(N)]
pCumSum_tot <- ggplot(taxcumsum_tot, aes(TotalCounts, CumSum)) + 
  geom_point() +
  xlab("Filtering threshold, minimum total counts") +
  ylab("Number of OTUs filtered") +
  scale_y_continuous(trans='log10') +
  ggtitle("Number of OTUs that would be filtered vs. the minimum total count filtering threshold") +
  xlim(0, 25) # Optional: zoom in on bend
ggsave(file.path(plot_outdir, "total_abundance_filtering_threshold.png"), plot = pCumSum_tot)
# --> 8 used by pipeline, no additional threshold

# Plot cumulative relative OTU sum to determine relative count filtering threshold
taxcumsum_rel = tdt_rel[, .N, by = RelCounts]
setkey(taxcumsum_rel, RelCounts)
taxcumsum_rel[, CumSum := cumsum(N)]
plot_ly(
  data = as.data.frame(taxcumsum_rel),
  x = ~RelCounts,
  y = ~CumSum,
  type = 'scatter',
  mode = 'markers'
) %>%
  layout(
    xaxis = list(
      title = "Filtering threshold in percentage, minimum relative counts",
      range = list(0, 0.0003),
      tickformat = ".9f"
    ),
    yaxis = list(title = "Number of OTUs filtered"),
    title = "Number of OTUs that would be filtered vs. the minimum relative count filtering threshold in percentage"
  )
pCumSum_rel_ggplot <- ggplot(taxcumsum_rel, aes(RelCounts, CumSum)) + 
  geom_point() +
  xlab("Filtering threshold in percentage, minimum relative counts") +
  ylab("Number of OTUs filtered") +
  ggtitle("Number of OTUs that would be filtered vs. the minimum relative count filtering threshold in percentage") +
  xlim(0, 0.0003) # Optional: zoom in on bend
ggsave(file.path(plot_outdir, "relative_abundance_filtering_threshold.png"), plot = pCumSum_rel_ggplot)
# --> 0.00000007 threshold

# Plot OTU prevalence (= in how many samples each OTU appears)
source("/Users/christopherhempel/GDrive private/Private documents/MyScripts/taxa_summary.R", local = TRUE)
mdt = fast_melt(physeq)
prevdt = mdt[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = TaxaID]
prevhist <- ggplot(prevdt, aes(Prevalence)) + 
  geom_histogram(fill = "blue", color = "black", alpha = 0.7, bins=265) + 
  xlim(0, 25) + # Optional: zoom in on peak
  labs(title = "Histogram of OTU prevalence", 
       x = "Prevalences", 
       y = "OTU count")
ggsave(file.path(plot_outdir, "otu_prevalence_histogram.png"), plot = prevhist)
# --> threshold prevalence > 2

# Num taxa appearing in 1 or 2 samples only
cat("Number of OTUs in only 1 sample:", prevdt[(Prevalence == 1), .N])
cat("Number of OTUs in only 2 samples:", prevdt[(Prevalence == 2), .N])

# OTU cumulative sum prevalence
prevcumsum = prevdt[, .N, by = Prevalence]
setkey(prevcumsum, Prevalence)
prevcumsum[, CumSum := cumsum(N)]
plot_ly(
  data = prevcumsum,
  x = ~Prevalence,
  y = ~CumSum,
  type = 'scatter',
  mode = 'markers'
) %>%
  layout(
    xaxis = list(title = "Filtering threshold, prevalence"),
    yaxis = list(title = "Number of OTUs filtered"),
    title = "Number of OTUs that would be filtered vs. the prevalence threshold"
  )
pPrevCumSum_ggplot <- ggplot(prevcumsum, aes(Prevalence, CumSum)) + 
  geom_point() +
  xlab("Filtering threshold, prevalence") +
  ylab("Number of OTUs filtered") +
  ggtitle("Number of OTUs that would be filtered vs. the prevalence threshold")
ggsave(file.path(plot_outdir, "otu_prevalence_cumulativesum.png"), plot = pPrevCumSum_ggplot)
# --> threshold prevalence > 2

# Prevalence vs. total abundance scatter plot
prev_vs_abundance <- phyloseq_prevalence_plot(
  physeq,
  prev.trh = 0.05, #Add horizontal line with prevalence threshold
  taxcolor = "phylum", #Taxonomy rank for coloring the points
  facet = F,
  point_alpha = 0.7,
  showplot = T
)
ggsave(file.path(plot_outdir, "prevalence_vs_abundance.png"), plot = prev_vs_abundance)
plot_ly(
  data = prevdt,
  x = ~Prevalence,
  y = ~TotalCounts,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 1, opacity = 0.6)
) %>%
  layout(
    xaxis = list(title = "Prevalence"),
    yaxis = list(title = "Total OTU abundance", type = "log"),
    title = "Prevalence vs. total abundance, OTUs"
  )

# - Alpha diversity across samples and based on depth category
plot_richness(physeq, measures=c("Shannon"))
p <- plot_richness(physeq, x="Depth.category", measures=c("Shannon"))
p + geom_boxplot(data = p$data, aes(x = Depth.category, y = value, color = NULL), alpha = 0.1)
ggsave(file.path(plot_outdir, "alpha_div_by_depths_cat.png"), plot = p)

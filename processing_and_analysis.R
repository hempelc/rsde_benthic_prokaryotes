library(plyr)
library(dplyr)
library(tidyverse)
library(vegan)
library(zCompositions)
library(phyloseq)
library(metagMisc)
library(data.table)
library(htmlwidgets)
library(iNEXT)
library(microViz)
library(ComplexHeatmap)
library(readxl)
library(GGally)
library(plotly)


# Define files
metadata_file <- "/Users/simplexdna/Desktop/RSDE BAC soil project/RSDE metabarcoding metadata.xlsx"
otu_table_file <- "/Users/simplexdna/Desktop/RSDE BAC soil project/apscale-new/rsde-bac-sediment_rerun_apscale_otu_table_filtered_microdecon-filtered_with_taxonomy.csv"

# Make folder for plots
plot_outdir <- "/Users/simplexdna/Desktop/RSDE BAC soil project/analysis/plots/fancy_plots"

# Define rank for analysis (options: "domain", "phylum", "class", "order", "family", "genus", "species", "OTUs")
rank = "phylum"

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
metadata <- subset(metadata, select = c("Leg", "Phase", "Latitude","Longitude", "Depth", "Discovery site locations", "Red Sea zone", "Zone & Discovery", "Phase & Discovery", "Zone & Phase"))
metadata$Leg <- as.factor(metadata$Leg)
metadata$Phase <- as.factor(metadata$Phase)
metadata$"Red Sea zone" <- as.factor(metadata$"Red Sea zone")
metadata$"Discovery site locations" <- as.factor(metadata$"Discovery site locations")
metadata$"Zone & Discovery" <- as.factor(metadata$"Zone & Discovery")
metadata$Latitude <- as.numeric(metadata$Latitude)
metadata$Longitude <- as.numeric(metadata$Longitude)
metadata$Depth_untransformed <- metadata$Depth

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

# Add new columns for dummy variables for categorical variables
metadata <- metadata %>%
  mutate(
    "Leg 1" = case_when(
      Leg == 1 ~ 1,
      TRUE ~ 0
    )
  )
metadata <- metadata %>%
  mutate(
    "Leg 2" = case_when(
      Leg == 2 ~ 1,
      TRUE ~ 0
    )
  )
metadata <- metadata %>%
  mutate(
    "Leg 3" = case_when(
      Leg == 3 ~ 1,
      TRUE ~ 0
    )
  )
metadata <- metadata %>%
  mutate(
    "Leg 4" = case_when(
      Leg == 4 ~ 1,
      TRUE ~ 0
    )
  )
metadata <- metadata %>%
  mutate(
    "Phase 1" = case_when(
      Phase == 1 ~ 1,
      TRUE ~ 0
    )
  )
metadata <- metadata %>%
  mutate(
    "Phase 2" = case_when(
      Phase == 2 ~ 1,
      TRUE ~ 0
    )
  )
metadata <- metadata %>%
  mutate(
    "Phase DSC" = case_when(
      Phase == "DSC" ~ 1,
      TRUE ~ 0
    )
  )
metadata <- metadata %>%
  mutate(
    "Depth epipelagic" = case_when(
      Depth <= 200 ~ 1,
      TRUE ~ 0
    )
  )
metadata <- metadata %>%
  mutate(
    "Depth mesopelagic" = case_when(
      Depth >= 200 & Depth <= 1000 ~ 1,
      TRUE ~ 0
    )
  )
metadata <- metadata %>%
  mutate(
    "Depth bathypelagic" = case_when(
      Depth > 1000 ~ 1,
      TRUE ~ 0
    )
  )
metadata <- metadata %>%
  mutate(
    "Depth no light" = case_when(
      Depth >= 200 ~ 1,
      TRUE ~ 0
    )
  )
metadata <- metadata %>%
  mutate(
    "Depth light" = case_when(
      Depth < 200 ~ 1,
      TRUE ~ 0
    )
  )
metadata$"Depth no light" <- as.factor(metadata$"Depth no light")
metadata$"Depth light" <- as.factor(metadata$"Depth light")
metadata$"Zone SRS" <- as.integer(metadata$"Red Sea zone" == "South Red Sea")
metadata$"Zone CSRS" <- as.integer(metadata$"Red Sea zone" == "Central-south Red Sea")
metadata$"Zone CNRS" <- as.integer(metadata$"Red Sea zone" == "Central-north Red Sea")
metadata$"Zone NRS" <- as.integer(metadata$"Red Sea zone" == "North Red Sea")
metadata$"Zone GOA" <- as.integer(metadata$"Red Sea zone" == "Gulf of Aqaba")
metadata$"DSC Afifi" <- as.integer(metadata$"Zone & Discovery" == "Al Afifi Brine pool")
metadata$"DSC Atlantis" <- as.integer(metadata$"Zone & Discovery" == "Atlantis 2 Brine pool")
metadata$"DSC Canyon" <- as.integer(metadata$"Zone & Discovery" == "Al Wajh Canyon")
metadata$"DSC Aqaba" <- as.integer(metadata$"Zone & Discovery" == "Aqaba Brine pool")


# - Check if metadata is correlated
correlation_df <- metadata %>%
  select_if(is.numeric) %>%
  cor() %>%
  abs() %>%
  melt() 
categorize <- function(value) {
  if (value == 1) {
    return(1)
  } else if (value >= 0.8 & value < 1) {
    return(0.8)
  } else if (value >= 0.6 & value < 0.8) {
    return(0.6)
  } else if (value >= 0.4 & value < 0.6) {
    return(0.4)
  } else if (value >= 0.2 & value < 0.4) {
    return(0.2)
  } else {
    return(0)
  }
}
correlation_df$"Correlation_category" <- sapply(abs(correlation_df$value), categorize)
correlation_heatmap <- plot_ly(correlation_df, x = ~Var1, y = ~Var2, z = ~Correlation_category, type = "heatmap", colors = colorRamp(colorRampPalette(c("white", "red"))(6)))
correlation_heatmap

heatmap(abs(cor(select_if(metadata, is.numeric))), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topright", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))

metadata %>%
  select_if(is.numeric) %>%
  ggpairs()
################################
# - Read in OTU table and create the output folder
otu_df <- read.csv(otu_table_file, row.names="ID")
dir.create(plot_outdir, showWarnings = FALSE)

# - Make a phyloseq object. Therefore, we preprocess the df, split it into OTU table
# - and taxonomy table, and set up the separate objects for phyloseq, including a metadata object:
# Standardize taxonomy for anything not identifiable. NOTE: This is optional!
# I do this to put everything that could not be identified into the same category
# to delete it later, but the data can also be analyzed non-standardized to get
# a finer distinction between non-assigned OTUs
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


############################
# - Data processing starts here
# - Now that we have an idea of how our data looks, we can process it, for example by
# - using some of the filter tresholds we explored earlier. Note that both the OTU table
# - AND the tax table of the phyloseq object are automatically adapted/changed.
# - This is a major advantage of phyloseq objects:
# Threshold rel abundance > 0.00000007 (determined in general_data_exploration.R script)
thresh_relabun = 0.00000007
keepTaxa = (taxa_sums(physeq) / sum(taxa_sums(physeq))) > thresh_relabun
physeq = prune_taxa(keepTaxa, physeq)
# Threshold prevalence >2  (determined in general_data_exploration.R script)
tresh_prev = 2
physeq <- filter_taxa(physeq, function(x) sum(x > 0) > tresh_prev, prune = T)


########################
# - The fun part - data visualization
# - Another major advantage of phyloseq is that it is very easy to generate nice visualizations,
# - since the data is highly organized and can be easily modified for graphs.
# - On top of that, there are great R packages that build on top of phyloseq and allow
# an even easier generation of graphs. One such package is microViz.

# - For example, check out the documentation of microViz on visualizing
# - compositions: https://david-barnett.github.io/microViz/articles/web-only/compositions.html

# - And here is the documentation to generate very nice heatmaps:
# - https://david-barnett.github.io/microViz/articles/web-only/heatmaps.html

# - With those documentations (and other material you can find on the internet by Google-ing
# - "phyloseq visualizations"), you are well equipped to explore your data yourself! Below,
# - I will showcase 3 plots/functions that I find really useful.

# - First, though, we need to edit our tax table slighly to be able to aggregate the taxonomy
# - for plotting, like so:
physeq <- tax_fix(physeq, unknowns = c("Incertae_Sedis", "Incertae_Sedis class", "Gammaproteobacteria_Incertae_Sedis", "Alphaproteobacteria_Incertae_Sedis", "Unknown_Family"))


###################### BARPLOTS
### Colour stuff
cols = c("#0169fd",
          "#f0da00",
          "#5d3c00",
          "#ffa997",
          "#005869",
          "#f17bff",
          "#204800",
          "#eaffca",
         "#39a2ff",
          "#550033",
          "#62cc00",
          "#00100d",
         "#c3008c",
          "lightgrey")

[
  "#65D9FF",
  "#FFD400",
  "#DD9C00",
  "#FFFDFF",
  "#65C8D9",
  "#FF7BFF",
  "#64C800",
  "#FEFFFA",
  "#99FFFF",
  "#B100A3",
  "#C2FF00",
  "#65A0AD",
  "#D400FF",
  "#FFFFFF"
]
[
  "#0000BD",
  "#D0CA00",
  "#000000",
  "#DF8987",
  "#000000",
  "#D17BFF",
  "#104000",
  "#D8FF00",
  "#000002",
  "#350033",
  "#00A000",
  "#000000",
  "#A3008C",
  "#999999"
]



Copy code
"#FFD400", "#0169FD", "#0000BD"
"#DD9C00", "#f0da00", "#D0CA00"
"#FFFDFF", "#5d3c00", "#000000"
"#65C8D9", "#ffa997", "#DF8987"
"#FF7BFF", "#005869", "#000000"
"#64C800", "#f17bff", "#D17BFF"
"#FEFFFA", "#204800", "#104000"
"#99FFFF", "#eaffca", "#D8FF00"
"#B100A3", "#39a2ff", "#000002"
"#C2FF00", "#550033", "#350033"
"#65A0AD", "#62cc00", "#00A000"
"#D400FF", "#00100d", "#000000"
"#FFFFFF", "#c3008c", "#A3008C"


# Colours for nested hierarchy
hueRank <- "phylum"
hueRankPlural <- "Phyla"
shadeRank <- "class"
nHues <- 13 # "Other" phyla will be shades of grey
nShades <- 3 # "Other" families will be the lightest shade of each hue

### Plots
# Plot barplots merged by zones
merged_physeq_zones <- physeq %>%
  merge_samples("Zone...Discovery") 
sample_data(merged_physeq_zones)$"Zone...Discovery" <- rownames(sample_data(merged_physeq_zones))
## On phylum
bar_chart_merged_zones_phylum <- merged_physeq_zones %>%
  comp_barplot(
    tax_level = "phylum", # define rank to aggregate on
    n_taxa =13, # define the number of taxa to give unique colours
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other phyla", # set custom name for the "other" category
    bar_width = 0.8, # reduce the bar width to 70% of one row
    bar_outline_colour = NA, # is the default (use NA to remove outlines)
    palette = cols
  ) +
  coord_flip()
bar_chart_merged_zones_phylum
ggsave(file.path(plot_outdir, "bar_chart_merged_zones_phylum.svg"), plot = bar_chart_merged_zones_phylum)
ggsave(units="px", width=8000, height=3000, dpi=600, file.path(plot_outdir, "bar_chart_merged_zones_phylum.png"), plot = bar_chart_merged_zones_phylum)
ggsave(units="px", width=8000, height=3000, dpi=600, file.path(plot_outdir, "bar_chart_merged_zones_phylum.svg"), plot = bar_chart_merged_zones_phylum)


## On class - edit palette manually
hierarchicalPalInfo_physeq <- merged_physeq_zones %>%
  tax_sort(by = sum, at = shadeRank) %>%
  tax_sort(by = sum, at = hueRank) %>%
  tax_agg(rank = shadeRank)
  hierarchicalPalInfo <- data.frame(
  hue = as.vector(tt_get(hierarchicalPalInfo_physeq)[, hueRank]),
  shade = as.vector(tt_get(hierarchicalPalInfo_physeq)[, shadeRank]),
  counts = taxa_sums(otu_get(hierarchicalPalInfo_physeq))
)
hierarchicalPalInfo <- hierarchicalPalInfo %>%
  dplyr::mutate(
    hue = forcats::fct_other(
      f = hue, keep = unique(hue)[seq_len(nHues)],
      other_level = paste("Other", hueRankPlural)
    ),
    nChrHue = nchar(as.character(hue)), padHue = max(nChrHue) - nChrHue
  ) %>%
  dplyr::group_by(hue) %>%
  dplyr::mutate(
    shade = forcats::fct_other(
      f = shade, keep = unique(shade)[seq_len(nShades - 1)],
      other_level = "Other"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    nChrShade = nchar(as.character(shade)), padShade = max(nChrShade) - nChrShade,
    Taxa = paste0(hue, ": ", strrep(" ", padHue), shade, strrep(" ", padShade))
  )
hierarchicalPalMatrix <- matrix(
  data = sapply(
    X = seq(from = 30, to = 75, length.out = nShades),
    FUN = function(l) scales::hue_pal(l = l, h.start = 30)(n = nHues)
  ),
  byrow = TRUE, ncol = nHues
)
hierarchicalPalMatrix <- cbind(hierarchicalPalMatrix, grey.colors(n = nShades))
hierarchicalPal <- hierarchicalPalMatrix %>%
  as.vector() %>%
  setNames(unique(hierarchicalPalInfo$Taxa))

bar_chart_merged_zones_class <- hierarchicalPalInfo_physeq %>%
  tax_mutate("Phylum: Class" = hierarchicalPalInfo$Taxa, .keep = "none") %>%
  comp_barplot(
    tax_level = "Phylum: Class", n_taxa = length(hierarchicalPal),
    tax_order = "asis",
    palette = hierarchicalPal,
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other phyla", # set custom name for the "other" category
    bar_width = 0.8, # reduce the bar width to 70% of one row
    bar_outline_colour = NA, # is the default (use NA to remove outlines)
  ) +
  coord_flip()
bar_chart_merged_zones_class
ggsave(units="px", width=8000, height=3000, dpi=600, file.path(plot_outdir, "bar_chart_merged_zones_class.png"), plot = bar_chart_merged_zones_class)
ggsave(units="px", width=8000, height=3000, dpi=600, file.path(plot_outdir, "bar_chart_merged_zones_class.svg"), plot = bar_chart_merged_zones_class)


# Plot barplots merged by depth categories
merged_physeq_depth <- physeq %>%
  merge_samples("Depth.category") 
sample_data(merged_physeq_depth)$"Depth.category" <- rownames(sample_data(merged_physeq_depth))
## On phylum
bar_chart_merged_depth_phylum <- merged_physeq_depth %>%
  comp_barplot(
    tax_level = "phylum", # define rank to aggregate on
    n_taxa =13, # define the number of taxa to give unique colours
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other phyla", # set custom name for the "other" category
    bar_width = 0.8, # reduce the bar width to 70% of one row
    bar_outline_colour = NA, # is the default (use NA to remove outlines)
    palette = cols
  ) +
  coord_flip()
bar_chart_merged_depth_phylum
ggsave(units="px", width=8000, height=3000, dpi=600, file.path(plot_outdir, "bar_chart_merged_depth_phylum.png"), plot = bar_chart_merged_depth_phylum)
ggsave(units="px", width=8000, height=3000, dpi=600, file.path(plot_outdir, "bar_chart_merged_depth_phylum.svg"), plot = bar_chart_merged_depth_phylum)


## On class - edit palette manually
hierarchicalPalInfo_physeq <- merged_physeq_depth %>%
  tax_sort(by = sum, at = shadeRank) %>%
  tax_sort(by = sum, at = hueRank) %>%
  tax_agg(rank = shadeRank)
hierarchicalPalInfo <- data.frame(
  hue = as.vector(tt_get(hierarchicalPalInfo_physeq)[, hueRank]),
  shade = as.vector(tt_get(hierarchicalPalInfo_physeq)[, shadeRank]),
  counts = taxa_sums(otu_get(hierarchicalPalInfo_physeq))
)
hierarchicalPalInfo <- hierarchicalPalInfo %>%
  dplyr::mutate(
    hue = forcats::fct_other(
      f = hue, keep = unique(hue)[seq_len(nHues)],
      other_level = paste("Other", hueRankPlural)
    ),
    nChrHue = nchar(as.character(hue)), padHue = max(nChrHue) - nChrHue
  ) %>%
  dplyr::group_by(hue) %>%
  dplyr::mutate(
    shade = forcats::fct_other(
      f = shade, keep = unique(shade)[seq_len(nShades - 1)],
      other_level = "Other"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    nChrShade = nchar(as.character(shade)), padShade = max(nChrShade) - nChrShade,
    Taxa = paste0(hue, ": ", strrep(" ", padHue), shade, strrep(" ", padShade))
  )
hierarchicalPalMatrix <- matrix(
  data = sapply(
    X = seq(from = 30, to = 75, length.out = nShades),
    FUN = function(l) scales::hue_pal(l = l, h.start = 30)(n = nHues)
  ),
  byrow = TRUE, ncol = nHues
)
hierarchicalPalMatrix <- cbind(hierarchicalPalMatrix, grey.colors(n = nShades))
hierarchicalPal <- hierarchicalPalMatrix %>%
  as.vector() %>%
  setNames(unique(hierarchicalPalInfo$Taxa))

bar_chart_merged_depth_class <- hierarchicalPalInfo_physeq %>%
  tax_mutate("Phylum: Class" = hierarchicalPalInfo$Taxa, .keep = "none") %>%
  comp_barplot(
    tax_level = "Phylum: Class", n_taxa = length(hierarchicalPal),
    tax_order = "asis",
    palette = hierarchicalPal,
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other phyla", # set custom name for the "other" category
    bar_width = 0.8, # reduce the bar width to 70% of one row
    bar_outline_colour = NA, # is the default (use NA to remove outlines)
  ) +
  coord_flip()
bar_chart_merged_depth_class
ggsave(units="px", width=8000, height=3000, dpi=600, file.path(plot_outdir, "bar_chart_merged_depth_class.png"), plot = bar_chart_merged_depth_class)
ggsave(units="px", width=8000, height=3000, dpi=600, file.path(plot_outdir, "bar_chart_merged_depth_class.svg"), plot = bar_chart_merged_depth_class)

#####################
# - Annotated heatmap of the 8 most abundant phyla:
heatmap_top_all_vars <- physeq %>%
  tax_transform("compositional", rank = "phylum") %>% # Turn data into relative abundances
  comp_heatmap(
    grid_col = NA,
    taxa = tax_top(physeq, 8, by = max, rank = "phylum"),
    cluster_rows = FALSE,
#    sample_anno = sampleAnnotation(
#      Depth.category = anno_sample_cat(var = "Depth.category", col = c("#b64773", "darkblue", "darkgreen", "gold"), box_col = NA, border_col = "black"),
#      Leg = anno_sample_cat(var = "Leg", col = c("orange", "lightblue", "pink", "#71a44c"), box_col = NA, border_col = "black"),
#      Phase = anno_sample_cat(var = "Phase", col = c("purple", "yellow", "darkblue"), box_col = NA, border_col = "black"),
#      Lat = anno_sample("Latitude"), 
#      Lon = anno_sample("Longitude"),
#      Depth = anno_sample("Depth"),
#      Leg1 = anno_sample_cat(var = "Leg.1", col = c("orange", "darkblue"), box_col = NA, border_col = "black"),
#      Leg2 = anno_sample_cat(var = "Leg.2", col = c("darkblue", "orange"), box_col = NA, border_col = "black"),
#      Leg3 = anno_sample_cat(var = "Leg.3", col = c("darkblue", "orange"), box_col = NA, border_col = "black"),
#      Leg4 = anno_sample_cat(var = "Leg.4", col = c("darkblue", "orange"), box_col = NA, border_col = "black"),
#      Phase1 = anno_sample_cat(var = "Phase.1", col = c("#007127", "#ff9a6f"), box_col = NA, border_col = "black"),
#      Phase2 = anno_sample_cat(var = "Phase.2", col = c("#ff9a6f", "#007127"), box_col = NA, border_col = "black"),
#      PhaseDSC = anno_sample_cat(var = "Phase.DSC", col = c("#ff9a6f", "#007127"), box_col = NA, border_col = "black"),
#      Depth.upper.photic = anno_sample_cat(var = "Depth.upper.photic", col = c("lightgrey", "black"), box_col = NA, border_col = "black"),
#      Depth.lower.photic = anno_sample_cat(var = "Depth.lower.photic", col = c("lightgrey", "black"), box_col = NA, border_col = "black"),
#      Depth.dysphotic = anno_sample_cat(var = "Depth.dysphotic", col = c("lightgrey", "black"), box_col = NA, border_col = "black"),
#      Depth.aphotic = anno_sample_cat(var = "Depth.aphotic", col = c("black", "lightgrey"), box_col = NA, border_col = "black"),
#      Depth.light = anno_sample_cat(var = "Depth.light", col = c("lightblue", "purple"), box_col = NA, border_col = "black"),
#      Depth.no.light = anno_sample_cat(var = "Depth.no.light", col = c("purple", "lightblue"), box_col = NA, border_col = "black"),
#      Zone = anno_sample_cat(var = "Red.Sea.zone", col = c("orange", "lightblue", "pink", "#71a44c", "gold"), box_col = NA, border_col = "black"),
#      ZoneDisc = anno_sample_cat(var = "Zone...Discovery", col = c("orange", "lightblue", "pink", "#71a44c", "gold", "darkblue", "#b64773", "purple", "yellow"), box_col = NA, border_col = "black"),
#      ZoneSRS = anno_sample_cat(var = "Zone.SRS", col = c("darkblue", "orange"), box_col = NA, border_col = "black"),
#      ZoneCSRS = anno_sample_cat(var = "Zone.CSRS", col = c("darkblue", "orange"), box_col = NA, border_col = "black"),
#      ZoneCNRS = anno_sample_cat(var = "Zone.CNRS", col = c("orange", "darkblue"), box_col = NA, border_col = "black"),
#      ZoneNRS = anno_sample_cat(var = "Zone.NRS", col = c("darkblue", "orange"), box_col = NA, border_col = "black"),
#      ZoneGOA = anno_sample_cat(var = "Zone.GOA", col = c("darkblue", "orange"), box_col = NA, border_col = "black"),
#      DSCAtlantis = anno_sample_cat(var = "DSC.Atlantis", col = c("darkblue", "orange"), box_col = NA, border_col = "black"),
#      DSCAfifi = anno_sample_cat(var = "DSC.Afifi", col = c("darkblue", "orange"), box_col = NA, border_col = "black"),
#      DSCCanyon = anno_sample_cat(var = "DSC.Canyon", col = c("darkblue", "orange"), box_col = NA, border_col = "black"),
#      DSCAqaba = anno_sample_cat(var = "DSC.Aqaba", col = c("darkblue", "orange"), box_col = NA, border_col = "black"),
#      border = TRUE
#    )
  )
heatmap_top_all_vars
# Export manually

#############
# - Correlation heatmap
## Number of taxa to show
num_taxa <- 12

# The following lines are to fix the function because it broke somehow
cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
  grid.text(mat_letters[i, j], x, y)
}
mat_letters = matrix(sample(letters[1:4], 240, replace = TRUE), num_taxa)

# By depth
cor_heatmap_depth <- physeq %>%
  tax_transform("clr", rank = "phylum") %>% # Turn data into clr abundances
  cor_heatmap(taxa = tax_top(physeq, num_taxa, by = max, rank = "phylum"),
              cluster_rows = FALSE,
              tax_anno = NA,
              colors = heat_palette(palette = "PRGn", sym = TRUE, range=c(-1,1)),
              grid_col = NA,
              numbers = heat_numbers(decimals = 2, col = "black", fontface = "plain", fontsize=8),
              seriation_method_col = "Identity",
              vars = c("Depth", "Depth.epipelagic", "Depth.mesopelagic", "Depth.bathypelagic")
  )
cor_heatmap_depth
png(file.path(plot_outdir, "cor_heatmap_depth.png"), width = 350, height = 600, res=100)
cor_heatmap_depth
dev.off()
svglite::svglite(file.path(plot_outdir, "cor_heatmap_depth.svg"), width = 3.5, height = 6)
cor_heatmap_depth
dev.off()

# By zone
cor_heatmap_zones <- physeq %>%
  tax_transform("clr", rank = "phylum") %>% # Turn data into clr abundances
  cor_heatmap(taxa = tax_top(physeq, 12, by = max, rank = "phylum"),
              cluster_rows = FALSE,
              tax_anno = NA,
              colors = heat_palette(palette = "PRGn", sym = TRUE, range=c(-1,1)),
              grid_col = NA,
              numbers = heat_numbers(decimals = 2, col = "black", fontface = "plain", fontsize=8),
              seriation_method_col = "Identity",
              vars = c("Latitude", "Zone.GOA", "Zone.NRS", "Zone.CNRS", "Zone.CSRS", "Zone.SRS")
  )
png(file.path(plot_outdir, "cor_heatmap_zones.png"), width = 420, height = 600, res=100)
cor_heatmap_zones
dev.off()
svglite::svglite(file.path(plot_outdir, "cor_heatmap_zones.svg"), width = 3.5, height = 5)
cor_heatmap_zones
dev.off()


#################
# Determine statistically important variables for RDA
## Function fix
svd <- function (x, nu = min(n, p), nv = min(n, p), LINPACK = TRUE)
{
  print("LINPACK:"); print(LINPACK)  ## added so you can see it's changed
  x <- as.matrix(x)
  if (any(!is.finite(x)))
    stop("infinite or missing values in 'x'")
  dx <- dim(x)
  n <- dx[1L]
  p <- dx[2L]
  if (!n || !p)
    stop("a dimension is zero")
  La.res <- La.svd(x, nu, nv)   ## your problem line
  res <- list(d = La.res$d)
  if (nu)
    res$u <- La.res$u
  if (nv) {
    if (is.complex(x))
      res$v <- Conj(t(La.res$vt))
    else res$v <- t(La.res$vt)
  }
  res
}
assignInNamespace("svd", svd, "base")

## Variable importance tests
physeq_no_dsc <- physeq %>% ps_filter(Phase != "DSC")
rda_env_vars <- metadata %>% select_if(is.numeric) %>% arrange(row.names(metadata))
rda_env_vars_no_dsc <- select_if(data.frame(sample_data(physeq_no_dsc)), is.numeric) %>% arrange(row.names(data.frame(sample_data(physeq_no_dsc))))
# On OTU level
rda_data <- t(as.data.frame(otu_table(tax_transform(physeq, "clr"))))
## Forward selection of variables:
fwd_sel_otu <- ordiR2step(rda(rda_data ~ 1, data = rda_env_vars), # lower model limit (simple!)
                      rda(rda_data ~ ., data = rda_env_vars)) # upper model limit (the "full" model)
# Result BEFORE CLR: Leg 1` + Depth + `Phase 1` + `Leg 2
# Redo selection without DSCs
rda_data <- t(as.data.frame(otu_table(tax_transform(physeq_no_dsc, "clr"))))
## Forward selection of variables:
fwd_sel_otu_noDSC <- ordiR2step(rda(rda_data ~ 1, data = rda_env_vars_no_dsc), # lower model limit (simple!)
                             rda(rda_data ~ ., data = rda_env_vars_no_dsc)) # upper model limit (the "full" model)
# Result: 
# Depth.no.light + Leg.2 + Phase.2 + Depth + 
# Depth.upper.photic + Zone.CSRS + Zone.SRS + Zone.CNRS + Zone.NRS

# On phylum level
rda_data <- t(as.data.frame(otu_table(tax_transform(physeq, "clr", rank = "phylum"))))
## Forward selection of variables:
fwd_sel_phylum <- ordiR2step(rda(rda_data ~ 1, data = rda_env_vars), # lower model limit (simple!)
                      rda(rda_data ~ ., data = rda_env_vars)) # upper model limit (the "full" model)
# Result:
# Longitude + Depth + `Depth no light` + 
#`DSC Canyon` + `Zone GOA` + `Phase 1` + `DSC Atlantis` + 
# `Leg 2` + `DSC Afifi` + `Zone SRS` + `Depth dysphotic` + 
# `Zone NRS`

# Redo selection without DSCs
rda_data <- t(as.data.frame(otu_table(tax_transform(physeq_no_dsc, "clr", rank = "phylum"))))
## Forward selection of variables:
fwd_sel_phylum_noDSC <- ordiR2step(rda(rda_data ~ 1, data = rda_env_vars_no_dsc), # lower model limit (simple!)
                      rda(rda_data ~ ., data = rda_env_vars_no_dsc)) # upper model limit (the "full" model)
# Result:
# Longitude + Depth.light + Depth + Zone.GOA + 
# Phase.1 + Leg.2 + Zone.SRS + Zone.NRS + Depth.dysphotic,

# On class level
rda_data <- t(as.data.frame(otu_table(tax_transform(physeq, "clr", rank = "class"))))
## Forward selection of variables:
fwd_sel_class <- ordiR2step(rda(rda_data ~ 1, data = rda_env_vars), # lower model limit (simple!)
                             rda(rda_data ~ ., data = rda_env_vars)) # upper model limit (the "full" model)
# Result:
# Longitude + `Phase 1` + `Depth no light` + 
#`Zone GOA` + `DSC Canyon` + Depth + `DSC Atlantis` + `Leg 1` + 
#  `Leg 2` + `Zone SRS` + `Depth dysphotic` + `Zone NRS` + `DSC Afifi`

# Redo selection without DSCs
rda_data <- t(as.data.frame(otu_table(tax_transform(physeq_no_dsc, "clr", rank = "class"))))
## Forward selection of variables:
fwd_sel_class_noDSC <- ordiR2step(rda(rda_data ~ 1, data = rda_env_vars_no_dsc), # lower model limit (simple!)
                                   rda(rda_data ~ ., data = rda_env_vars_no_dsc)) # upper model limit (the "full" model)
# Result:
# Longitude + Depth.no.light + Depth + 
# Zone.GOA + Phase.2 + Leg.1 + Leg.2 + Zone.SRS + Zone.NRS


###########################
# RDAs
## Explore
physeq %>%
  tax_transform(rank = "phylum", trans = "clr") %>%
  ord_explore()

# RDA coloured by Red Sea Zone
rda_zone <- physeq %>%
  tax_transform(rank = "phylum", trans = "clr") %>%
  ord_calc(
    constraints = c("Latitude", "Depth"),
    method = "RDA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Red.Sea.zone", fill = "Red.Sea.zone",
    shape = "Red.Sea.zone", alpha = 0.8,
    size = 2
  ) +
  ggside::geom_xsideboxplot(aes(fill = Red.Sea.zone, y = Red.Sea.zone), orientation = "y", varwidth=FALSE, show.legend=FALSE) +
  ggside::geom_ysideboxplot(aes(fill = Red.Sea.zone, x = Red.Sea.zone), orientation = "x", varwidth=FALSE, show.legend=FALSE) +
  ggside::theme_ggside_void()
ggsave(units="px", width=2000, height=1400, file.path(plot_outdir, "rda_zone.png"), plot = rda_zone)
ggsave(units="px", width=2000, height=1400, file.path(plot_outdir, "rda_zone.svg"), plot = rda_zone)

# RDA coloured by Depth 
rda_depth <- physeq %>%
  tax_transform(rank = "phylum", trans = "clr") %>%
  ord_calc(
    constraints = c("Latitude", "Depth"),
    method = "RDA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Depth_untransformed", fill = "Depth_untransformed",
    shape = "Depth.category", alpha = 0.8,
    size = 2
  )  +
  scale_colour_viridis_c(direction=-1) +
  ggside::geom_xsideboxplot(aes(group = Depth.category, y = Depth_untransformed), orientation = "y", varwidth=FALSE) +
  ggside::geom_ysideboxplot(aes(group = Depth.category, x = Depth_untransformed), orientation = "x", varwidth=FALSE) +
  ggside::theme_ggside_void()
ggsave(units="px", width=2000, height=1350, file.path(plot_outdir, "rda_depth.png"), plot = rda_depth)
ggsave(units="px", width=2000, height=1350, file.path(plot_outdir, "rda_depth.svg"), plot = rda_depth)

# RDA highlighting DSC sites
rda_dsc <- physeq %>%
  tax_transform(rank = "phylum", trans = "clr") %>%
  ord_calc(
    constraints = c("Latitude", "Depth"),
    method = "RDA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Discovery.site.locations", fill = "Discovery.site.locations",
    shape = "Discovery.site.locations", alpha = "Phase.DSC",
    size = 2
  ) + 
  scale_shape_girafe_filled()
ggsave(units="px", width=2000, height=1350, file.path(plot_outdir, "rda_dsc.png"), plot = rda_dsc)
ggsave(units="px", width=2000, height=1350, file.path(plot_outdir, "rda_dsc.svg"), plot = rda_dsc)

#### Optional:
# OTUs Red Sea Zone
physeq %>%
  tax_transform(trans = "clr") %>%
  ord_calc(
    constraints = c("Latitude", "Depth"),
    method = "RDA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Red.Sea.zone", fill = "Red.Sea.zone",
    shape = "Red.Sea.zone", alpha = 0.5,
    size = 2
  ) + 
  scale_shape_girafe_filled() 

# OTUs by Depth
physeq %>%
  tax_transform(trans = "clr") %>%
  ord_calc(
    constraints = c("Latitude", "Depth"),
    method = "RDA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Depth", fill = "Depth",
    shape = "Depth.category", alpha = 0.5,
    size = 2
  ) + 
  scale_shape_girafe_filled() 
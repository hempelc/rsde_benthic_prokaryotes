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
metadata_file <- "/Users/simplexdna/Desktop/RSDE BAC soil project/RSDE metabarcoding metadata with env noncor.csv"
otu_table_file <- "/Users/simplexdna/Desktop/RSDE BAC soil project/apscale-new/rsde-bac-sediment_rerun_apscale_otu_table_filtered_microdecon-filtered_with_taxonomy.csv"

# Make folder for plots
plot_outdir <- "/Users/simplexdna/Desktop/RSDE BAC soil project/analysis/plots/fancy_plots"
dir.create(plot_outdir, showWarnings = FALSE)

#############################  Metadata processing
# Make sure metadata matches otu table
# otu_df_transposed <- otu_df %>%
#   t() %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "index")
# merged_df <- merge(otu_df_transposed, metadata, by.x = "index", by.y = "Sample ID", all = TRUE)

# Read in metadata
metadata <- as.data.frame(read_csv(metadata_file))
rownames(metadata) <- metadata[, "Sample ID"]
metadata <- metadata[, -(which(names(metadata) == "Sample ID"))] %>% drop_na()

keep <- rownames(metadata)

# Cut down metadata and process data types
metadata$'Red_Sea_zone' <- factor(metadata$'Red_Sea_zone', levels=c("Gulf of Aqaba", "North Red Sea", "Central-north Red Sea", "Central-south Red Sea", "South Red Sea"))

metadata$"Zone & Discovery" <- as.factor(metadata$"Zone & Discovery")
metadata$Latitude <- as.numeric(metadata$Latitude)
metadata$Depth_untransformed <- metadata$Depth

# Add column with depth categories based on Larissa's recommendations
metadata <- metadata %>%
  mutate(
    "Depth_category" = case_when(
      Depth > 1000 ~ "Bathypelagic (>1000m)",
      Depth >= 200 & Depth <= 1000 ~ "Mesopelagic (200-1000m)",
      Depth <= 200 ~ "Epipelagic (<200m)",
      TRUE ~ NA_character_  # default condition (optional)
    )
  )
metadata$'Depth_category' <- factor(metadata$'Depth_category', levels=c("Epipelagic (<200m)", "Mesopelagic (200-1000m)", "Bathypelagic (>1000m)"))

# Add new columns for dummy variables for categorical variables
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
metadata$"Zone SRS" <- as.integer(metadata$"Red_Sea_zone" == "South Red Sea")
metadata$"Zone CSRS" <- as.integer(metadata$"Red_Sea_zone" == "Central-south Red Sea")
metadata$"Zone CNRS" <- as.integer(metadata$"Red_Sea_zone" == "Central-north Red Sea")
metadata$"Zone NRS" <- as.integer(metadata$"Red_Sea_zone" == "North Red Sea")
metadata$"Zone GOA" <- as.integer(metadata$"Red_Sea_zone" == "Gulf of Aqaba")
#metadata$"DSC Atlantis" <- as.integer(metadata$"Zone & Discovery" == "Atlantis 2 Brine pool")
metadata$"DSC Canyon" <- as.integer(metadata$"Zone & Discovery" == "Al Wajh Canyon")
metadata$"DSC Aqaba" <- as.integer(metadata$"Zone & Discovery" == "Aqaba Brine pool")
metadata$"DSC Afifi" <- as.integer(metadata$"Zone & Discovery" == "Al Afifi Brine pool")
metadata$"DSC" <- as.integer(metadata$"Discovery site locations" != "Non_discovery")

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
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
otu_df <- read_csv(otu_table_file, col_select = c(keep, ranks))

# - Make a phyloseq object. Therefore, we preprocess the df, split it into OTU table
# - and taxonomy table, and set up the separate objects for phyloseq, including a metadata object:
# Standardize taxonomy for anything not identifiable. NOTE: This is optional!
# I do this to put everything that could not be identified into the same category
# to delete it later, but the data can also be analyzed non-standardized to get
# a finer distinction between non-assigned OTUs
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
physeq <- tax_fix(physeq, unknowns = c("Incertae_Sedis", "Incertae_Sedis class", "Gammaproteobacteria_Incertae_Sedis", "Alphaproteobacteria_Incertae_Sedis", "Unknown_Family"))

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

# By other env variables
cor_heatmap_env <- physeq %>%
  tax_transform("clr", rank = "phylum") %>% # Turn data into clr abundances
  cor_heatmap(taxa = tax_top(physeq, num_taxa, by = max, rank = "phylum"),
              cluster_rows = FALSE,
              tax_anno = NA,
              colors = heat_palette(palette = "PRGn", sym = TRUE, range=c(-1,1)),
              grid_col = NA,
              numbers = heat_numbers(decimals = 2, col = "black", fontface = "plain", fontsize=8),
              seriation_method_col = "Identity",
              vars = c("Oxygen_saturation_Mean")
  )
cor_heatmap_env
png(file.path(plot_outdir, "cor_heatmap_env_v2.png"), width = 350, height = 600, res=100)
cor_heatmap_env
dev.off()
svglite::svglite(file.path(plot_outdir, "cor_heatmap_env_v2.svg"), width = 3.5, height = 6)
cor_heatmap_env
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
remove_vars <- c("SpeedofSound_Mean", "Salinity_Mean", "SpecificConductivity_Mean", "Chlorophyll_Mean", "Depth_untransformed", "Depth epipelagic", "Depth mesopelagic", "Depth bathypelagic", "Zone SRS", "Zone CSRS", "Zone CNRS", "Zone NRS", "Zone GOA", "DSC Atlantis", "DSC Canyon", "DSC Aqaba")
rda_env_vars <- metadata %>% select_if(is.numeric) %>% arrange(row.names(metadata))
rda_env_vars[remove_vars] <- NULL
rda_env_vars_no_dsc <- select_if(data.frame(sample_data(physeq_no_dsc)), is.numeric) %>% arrange(row.names(data.frame(sample_data(physeq_no_dsc))))

correlation_df_relevant <- rda_env_vars %>%
  cor() %>%
  abs() %>%
  melt() 
correlation_df_relevant$"Correlation_category" <- sapply(abs(correlation_df_relevant$value), categorize)
correlation_heatmap <- plot_ly(correlation_df_relevant, x = ~Var1, y = ~Var2, z = ~Correlation_category, type = "heatmap", colors = colorRamp(colorRampPalette(c("white", "red"))(6)))
correlation_heatmap

rda_env_vars %>%
    ggpairs()


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
    constraints = c("Latitude", "Depth", "Oxygen_saturation_Mean"),
    method = "RDA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Red_Sea_zone", fill = "Red_Sea_zone",
    shape = "Red_Sea_zone", alpha = 0.8,
    size = 2
  ) +
  ggside::geom_xsideboxplot(aes(fill = Red_Sea_zone, y = Red_Sea_zone), orientation = "y", varwidth=FALSE, show.legend=FALSE) +
  ggside::geom_ysideboxplot(aes(fill = Red_Sea_zone, x = Red_Sea_zone), orientation = "x", varwidth=FALSE, show.legend=FALSE) +
  ggside::theme_ggside_void()
# Somehow the sideboxplots don't work anymore
ggsave(units="px", width=2000, height=1400, file.path(plot_outdir, "rda_zone_v2.png"), plot = rda_zone)
ggsave(units="px", width=2000, height=1400, file.path(plot_outdir, "rda_zone_v2.svg"), plot = rda_zone)


# RDA coloured by Depth 
rda_depth <- physeq %>%
  tax_transform(rank = "phylum", trans = "clr") %>%
  ord_calc(
    constraints = c("Latitude", "Depth", "Oxygen_saturation_Mean"),
    method = "RDA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Depth_untransformed", fill = "Depth_untransformed",
    shape = "Depth_category", alpha = 0.8,
    size = 2
  )  +
  scale_colour_viridis_c(direction=-1) +
  ggside::geom_xsideboxplot(aes(group = Depth_category, y = Depth_untransformed), orientation = "y", varwidth=FALSE) +
  ggside::geom_ysideboxplot(aes(group = Depth_category, x = Depth_untransformed), orientation = "x", varwidth=FALSE) +
  ggside::theme_ggside_void()
# Somehow the sideboxplots don't work anymore
ggsave(units="px", width=2000, height=1350, file.path(plot_outdir, "rda_depth_v2.png"), plot = rda_depth)
ggsave(units="px", width=2000, height=1350, file.path(plot_outdir, "rda_depth_v2.svg"), plot = rda_depth)

# RDA highlighting DSC sites
rda_dsc <- physeq %>%
  tax_transform(rank = "phylum", trans = "clr") %>%
  ord_calc(
    constraints = c("Latitude", "Depth", "Oxygen_saturation_Mean"),
    method = "RDA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Discovery site locations", fill = "Discovery site locations",
    shape = "Discovery site locations", alpha = "DSC",
    size = 2
  ) + 
  scale_shape_girafe_filled()
ggsave(units="px", width=2000, height=1350, file.path(plot_outdir, "rda_dsc_v2.png"), plot = rda_dsc)
ggsave(units="px", width=2000, height=1350, file.path(plot_outdir, "rda_dsc_v2.svg"), plot = rda_dsc)

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
    colour = "Red_Sea_zone", fill = "Red_Sea_zone",
    shape = "Red_Sea_zone", alpha = 0.5,
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

rda_depth <- physeq %>%
  tax_transform(rank = "phylum", trans = "clr") %>%
  ord_calc(
    constraints = c("Latitude", "Depth", "Oxygen_saturation_Mean"),
    method = "RDA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Depth_untransformed", fill = "Depth_untransformed",
    shape = "Depth_category", alpha = 0.8,
    size = 2
  )  +
  scale_colour_viridis_c(direction=-1) +
  ggside::geom_xsideboxplot(aes(group = "Depth_category", y = Depth_untransformed), orientation = "y", varwidth=FALSE) +
  ggside::geom_ysideboxplot(aes(group = "Depth_category", x = Depth_untransformed), orientation = "x", varwidth=FALSE) +
  ggside::theme_ggside_void()

# RDA coloured by Red Sea Zone
rda_zone <- physeq %>%
  tax_transform(rank = "phylum", trans = "clr") %>%
  ord_calc(
    constraints = c("Latitude", "Depth", "Oxygen_saturation_Mean"),
    method = "RDA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Red Sea zone", fill = "Red Sea zone",
    shape = "Red Sea zone", alpha = 0.8,
    size = 2
  ) +
  ggside::geom_xsideboxplot(aes(fill = "Red Sea zone", y = "Red Sea zone"), orientation = "y", varwidth=FALSE, show.legend=FALSE) +
  ggside::geom_ysideboxplot(aes(fill = "Red Sea zone", x = "Red Sea zone"), orientation = "x", varwidth=FALSE, show.legend=FALSE) +
  ggside::theme_ggside_void()

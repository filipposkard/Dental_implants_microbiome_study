# library(lme4)
# library(emmeans)
# library(mia)
library("dplyr")
library(phyloseq)
# library("phyloseqCompanion")
# install.packages("microeco")
library("microeco")
library("magrittr")

project_dir = this.path::here(..=1)
in_dir=file.path(project_dir,"data/phyloseq")
results_dir = file.path(project_dir,"results/compare_color_complexes")
dir.create(results_dir, showWarnings = F)

ps = readRDS(file =file.path(in_dir,"/phyloseq_object_species.rds"))
colnames(tax_table(ps)) =  tools::toTitleCase(phyloseq::rank_names(ps))
ps.rel = ps %>% transform_sample_counts(function(x) {x/sum(x)})


ps_filtered <- filter_taxa(ps, function(x) sum(x > 0) > (0.25 * length(x)), TRUE)

print(paste("Original Taxa:", ntaxa(ps)))
print(paste("Filtered Taxa (Frequency > 0.25):", ntaxa(ps_filtered)))
# 1. Prepare your data ---------------------------------------------------------
# Assuming 'ps' is your phyloseq object with counts

# The methodology states: "OTUs whose frequency is less than 0.6 in all samples are eliminated"
# This typically refers to PREVALENCE (presence in X% of samples).
# A threshold of 0.6 (60%) is quite strict but we will follow the method text exactly.

# Convert phyloseq to microeco object (microtable)
dataset <- microtable$new(sample_table = as.data.frame(sample_data(ps_filtered)),
                          otu_table = as.data.frame(otu_table(ps_filtered)),
                          tax_table = as.data.frame(tax_table(ps_filtered)))

# 2. Filter OTUs (Frequency/Prevalence) ----------------------------------------
# "Eliminate OTUs whose frequency is less than 0.6"
# This means: Keep OTUs present in >= 60% of samples.

# We use the filter_features function in microeco (which works on the microtable)
# Note: microeco calculates 'frequency' as the proportion of samples.

print(paste("Number of taxa remaining after 0.6 frequency filter:", 
            nrow(dataset$tax_table)))

# 3. Calculate Correlations (Spearman) -----------------------------------------
# Initialize the trans_network class
t1 <- trans_network$new(dataset = dataset, 
                        cor_method = "spearman", 
                        taxa_level = "Species") # Or "OTU"/"ASV" if distinct




# CORRECT SEQUENCE (run in this exact order):
t1$cal_network(COR_p_thres = 0.05, COR_optimization = FALSE, COR_cut = 0.6)

# 1. FIRST: Network attributes (before modules)
t1$cal_network_attr()

# 2. THEN: Modules  
t1$cal_module()

# 3. NOW extract tables
head(t1$get_node_table())
head(t1$get_edge_table())


# 3. Now node and edge tables are populated
# After cal_network_attr(), cal_module(), etc.
head(t1$get_node_table())   # ✓ Node attributes (Degree, Betweenness, etc.)
head(t1$get_edge_table())   # ✓ Edge attributes (correlation values, p-values)

# 4. Visualization
t1$plot_network(method = "igraph",
                node_size_scale = TRUE,
                node_color = "Phylum")  # color nodes by Phylum

# 5. Save network for Gephi (optional)
t1$save_network(filepath = "network_spearman_0.6.graphml")

















# 4. Construct Network (Thresholding) ------------------------------------------
# Methodology: "R value greater than 0.6 and P-value less than 0.05"
# 1. Calculate topological properties (Degree, Betweenness, etc.)
# This is mandatory to populate the node tables!
t1$cal_network_attr()

# 2. Get the node table (Check if it's not NULL anymore)
# Note: In some microeco versions, this is stored in res_node_table directly
head(t1$res_node_table)

# 3. Save the network for Gephi (Best for visualization)
# This creates a file you can open in Gephi to make professional network plots
t1$save_network(filepath = "co_occurrence_network.graphml")

# 4. Quick Plot in R
t1$plot_network(
  method = "igraph", 
  node_size_scale = TRUE,    # Scale node size by degree (hub status)
  node_color = "Phylum"      # Color nodes by Phylum
)




#########


# 1. Calculate topological properties (Degree, Betweenness, etc.)
# This is mandatory to populate the node tables!
x <- t1$cal_network_attr()

# 2. Get the node table (Check if it's not NULL anymore)
# Note: In some microeco versions, this is stored in res_node_table directly
head(t1$res_node_table)

# 3. Save the network for Gephi (Best for visualization)
# This creates a file you can open in Gephi to make professional network plots
t1$save_network(filepath = "co_occurrence_network.graphml")

# 4. Quick Plot in R
t1$plot_network(
  method = "igraph", 
  node_size_scale = TRUE,    # Scale node size by degree (hub status)
  node_color = "Phylum"      # Color nodes by Phylum
)



########
t1$cal_network(COR_p_thres = 0.05, 
               COR_optimization = FALSE, # We manually set R threshold, so turn off auto-optimization
               COR_cut = 0.6)            # R > 0.6

# 5. Get Network Statistics ----------------------------------------------------
# Calculate topological properties (degree, modularity, etc.)
t1$cal_network_attr()

# Check the node and edge tables
head(t1$res_node_table) # Attributes of nodes (OTUs)
head(t1$res_edge_table) # The connections

# 6. Visualization -------------------------------------------------------------
# Simple plot in R
t1$plot_network(method = "igraph", 
                node_size_scale = TRUE, 
                node_color = "Phylum") # Color nodes by Phylum

# 7. Save for Gephi (Optional but recommended) ---------------------------------
# Network files (graphml) are best visualized in Gephi
t1$save_network(filepath = "network_spearman_0.6.graphml")


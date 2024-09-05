#Originally written by Matthew Nitshcke
#Modified by Annika Lamb for Lamb et al.

```{r}
# Packages required to analyse this data. If you do not have them, please install them.
library(tidyverse)
library(kmer)
library(devtools)
library(phangorn)
library(GUniFrac)
library(ggtree)
library(patchwork)
library(bioseq)
library(ape)
library(usedist)
library(cowplot)

# Custom functions for loading and handling fasta files
dna_to_DNAbin <- function (dna){
  DNAbin <- as_DNAbin(dna)
  names(DNAbin) <- names(dna)
  return(DNAbin)
}

read_fasta_df <- function (file = "") {
  fasta <- readLines(file)
  ind <- grep(">", fasta)
  s <- data.frame(ind = ind, from = ind + 1, to = c((ind - 
                                                       1)[-1], length(fasta)))
  seqs <- rep(NA, length(ind))
  for (i in 1:length(ind)) {
    seqs[i] <- paste(fasta[s$from[i]:s$to[i]], collapse = "")
  }
  tib <- tibble(label = gsub(">", "", fasta[ind]), sequence = seqs)
  return(tib)
}
```

# Read in and format the sequence tables

```{r}
# Import sequence table
seqs <- read_tsv("20210715_lamb/post_med_seqs/HeatWaveITS2.seqs.absolute.abund_and_meta.txt") %>%
  mutate(sample_name = as.factor(sample_name)) %>%
  filter(!(is.na(sample_name))) # remove the final empty row. ANY SAMPLES NOT TO BE ANALYSED ADD ANOTHER FILTER HERE

# Import sample_name metadata
meta <- read_csv("metadata.csv")

# Convert sequence table to long format
head(seqs)
seqs_long <- seqs %>%
  select(sample_name, `A1dh`:`D17e`) %>% # Select sample_names and the each column contain sequence count data
  pivot_longer(`A1dh`:`D17e`) %>% # make into long dataframe
  filter(value > 0) %>% # Remove zero values
  filter(!(is.na(sample_name))) %>% # Remove non-name samples (tutorial specific)
  left_join(., meta) # add in your metadata

# What is our per sample sequence depth?
sample_totals <- seqs_long %>%
  group_by(sample_name) %>%
  summarise(total = sum(value))

sample_totals

# Create a list of samples to analyse that are above a minimum threshold
keepers <- sample_totals %>%
  filter(total > 1000)
#Remove untracked sample from dataset (from unknown tank)
keepers<-keepers[!(keepers$sample_name=="TT2526"),]

#Remove samples with < 1000 reads
removed <- sample_totals %>%
  filter(total < 1000)
removed

# Filter out the failed samples and convert counts to relative abundance
seqs_long <- seqs_long %>%
  filter(sample_name %in% keepers$sample_name) %>%
  group_by(sample_name) %>%
  mutate(value_rel = value/sum(value)) %>% # Convert to relative abundance
  ungroup() %>%
  mutate(name = as.factor(name)) # Make sequence names a factor

seqs_long$value_rel

#Filter to samples (no controls)
seqs_long<-filter(seqs_long, ControlSample == "Sample")
head(seqs_long)

#Summary of sequencing data
reads_by_species <- seqs_long %>%
  group_by(Species, sample_name) %>%
  summarise(sum = sum(value)) %>%
  ungroup() %>%
  group_by(Species) %>%
  summarise(mean = mean(sum))
reads_by_species

reads_per_sample <- seqs_long %>%
  group_by(sample_name) %>%
  summarise(sum = sum(value)) 
reads_per_sample
min(reads_per_sample$sum)
max(reads_per_sample$sum)
totalreads <- sum(reads_per_sample$sum)
totalreads


# Create a random palette for each sequence
n <- length(levels(seqs_long$name))
seqs_pal = rainbow(n, s=.6, v=.9)[sample(1:n,n, replace = FALSE)]
names(seqs_pal) <- levels(seqs_long$name)

# Read in and format the profile tables

# Read in the profile data
profiles_raw <- read_tsv("20210715_lamb/its2_type_profiles/HeatWaveITS2.profiles.absolute.abund_and_meta.txt", skip = 6)%>%
  select(sample_name = `...2`, `C1-C1b-C1c-C42.2-C1bh-C1br-C1cb-C72k`:`D1-D1u-D2-D4-D2.2`) %>% # rename the second column as "sample_name" and select the profile columns
  filter(!is.na(sample_name))

profiles_raw

profiles_long <- profiles_raw %>%
  pivot_longer(`C1-C1b-C1c-C42.2-C1bh-C1br-C1cb-C72k`:`D1-D1u-D2-D4-D2.2`) %>% # Convert profile table to long format
  mutate(name = paste0("p_", name)) %>% # Add a p_ to the beginning of each profile (Some profiles are single sequence profiles and clash with the Sequence names)
  filter(sample_name %in% seqs_long$sample_name) %>% # Remove samples that do not appear in the Sequence dataframe above (were filtered out)
  mutate(value = as.numeric(value)) %>%
  filter(value > 0) %>% # Remove 0 abundance profiles
  group_by(sample_name) %>%
  mutate(sample_name = as.factor(sample_name),
         name = as.factor(name)) %>%
  ungroup() %>%
  left_join(., meta) # Add in metadata
head(profiles_long)

# What is the total number of profile-related sequences in each sample?
profiles_sum <- profiles_long %>%
  group_by(sample_name) %>%
  summarise(total = sum(value))
profiles_sum

# How many sequences in each sample are not part of a profile?
residual <- left_join(profiles_sum, sample_totals, by = "sample_name") %>%
  mutate(residual = total.y - total.x) %>%
  select(sample_name, value = residual) %>%
  mutate(name = as.factor("non-profile sequences")) %>%
  left_join(., meta)
residual

sample_totals

# Combine the profile and non-profile sequences
profile_data <- rbind(profiles_long, residual) %>%
  group_by(sample_name) %>%
  mutate(value_rel = value/sum(value)) # convert to relative abundance
profile_data

# Create palette for profiles (this is a darker palette)
n <- length(levels(profile_data$name))
profile_pal = rainbow(n, s=.6, v=.6)[sample(1:n,n, replace = FALSE)]
names(profile_pal) <- levels(profile_data$name)

# Read in the sequences and create the unifrac distances

# Import and filter sequence data
fasta <- read_fasta_df("20210715_lamb/post_med_seqs/HeatWaveITS2.seqs.fasta") %>%
  filter(label %in% seqs_long$name) %>% # keep only the sequences present in the sequence table above
deframe() %>%
  as_dna() # convert to bioseq DNA format

# Create a kmer-based pairwise sequence similarity matrix
kdist <- fasta %>%
  dna_to_DNAbin() %>%
  kdistance(k = 7, residues = "DNA", method = "edgar") %>% # test how different k size alters the interpretation
  as.matrix()
#Note - small K = less discriminating power. Large k = more discriminating power but significantly more compute time

# Use heirarchical clustering to create a phylogenetic tree based on the pairwise k-mer distances
tree <- kdist %>% phangorn::upgma()

# Create a wide count table to use in the weighted versions of the unifrac measurements
seqs_wide <- seqs_long %>%
  select(sample_name, name, value) %>% # need the sample names, the sequence names, and the count data
  pivot_wider(names_from = name, values_from = value) %>% # convert from long to wide
  mutate(across(everything(), replace_na, 0)) %>% # NA should = 0
  tibble::column_to_rownames(var = "sample_name") # sample names need to be column names

# Calculate the unifrac distance
unidist <- GUniFrac(seqs_wide, tree)
unifracs <- unidist$unifracs
du <- unifracs[, , "d_1"]		# Weighted UniFrac
#du <- unifracs[, , "d_UW"]		# Unweighted UniFrac	
# du <- unifracs[, , "d_VAW"]		# Variance adjusted weighted UniFrac
# du <- unifracs[, , "d_0"]     	# GUniFrac with alpha 0  
# du <- unifracs[, , "d_0.5"]   	# GUniFrac with alpha 0.5 
# NOTE - try different unifrac measurements and decide which one (information can be found by searching GUniFrac in the help bar) makes sense based upon your study goal!

# Hierarchical clustering of the samples into groups based on their pairwise unifrac distances
hclust_samps <- upgma(du)

# Plot

# Join profiles and sequence data together into single dataframe
all_data <- rbind(seqs_long, profile_data)

# Create a sample order according to their appearance in the clustering tree
ggtree_samps <- ggtree::ggtree(hclust_samps)$data %>%
  filter(isTip == "TRUE") %>%
  arrange(y)

# Generate a plotting dataframe
plot_df <- all_data %>%
  mutate(sample_name = fct_relevel(sample_name, ggtree_samps$label))

# Merge the palettes and replace the non-profile sequences with grey
all_pal <- c(seqs_pal, profile_pal)
all_pal['non-profile sequences'] <- "#808080"

# Generate the bar-plots
theme_set(theme_bw()) # set the theme

p_bar_uni <- ggplot(plot_df, aes(value_rel, sample_name)) +
  geom_bar(stat = "identity", aes(fill = name, color = name), position='fill',width = 1) +
  theme(aspect.ratio = 1, legend.position = "bottom", axis.text.y=element_blank(), axis.ticks.y = element_blank(),
        text = element_text(size=18),
        axis.text.x=element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size = 18))+
  scale_fill_manual(guide="none", values = all_pal, breaks = levels(profile_data$name)) +
  scale_colour_manual(guide="none", values = all_pal, breaks = levels(profile_data$name)) +
  geom_vline(xintercept = 0.5, size = 2) +
  guides(fill=guide_legend(ncol = 1, title = "Profile"))
  
p_bar_uni


# Use ggtree to visualise the clustering of samples
p_tree <- ggtree(hclust_samps, size = 0.2) +
  scale_x_reverse() +
  theme(aspect.ratio = 1.8)

# Add in the sample metadata to the tree
p_tree$data <- left_join(p_tree$data, meta, by = c("label" = "sample_name")) %>%
  mutate(Species = paste0(Species)) # can customise any metadata that is needed here

# Colour the tree tips by the species information
labels=c(expression(italic("A. loripes")), expression(italic("A. kenti")), "KL hybrid")
p_tree_tip <- p_tree + 
  geom_tippoint(aes(color = Species), size = 5) + scale_color_manual(values=c("#d95f02","#7570b3","#1b9e77"),labels = labels) +labs(color = "Offspring group")+theme(text = element_text(size=18),legend.text = element_text(size = 18))
p_tree_tip



# Plot the profile + sequence bar plots and hclust sample tree together

aligned_plots <- plot_grid(
   p_bar_uni, 
   p_tree_tip,
  ncol = 2, 
  align = "h", 
  axis = "lr",
  rel_heights = c(1,1),
  greedy = TRUE,
  labels = c('A','B'),
  label_size = 20
) 


# Display the combined plot
print(aligned_plots)

```

# PERMANOVA of metadata factors and their explaining power of the unifrac distances

```{r}
library(vegan)


# Convert unifrac matrix into dist class
dist <- as.dist(du, diag = FALSE)


# Check the dimensions of the matrix (should equal n x n samples of interest)
dim(du)

# format a metadata dataframe for the permanova
meta
adonis_meta <- meta %>%
  select(sample_name, Species, Tank, Tile, Treatment) %>%
  filter(sample_name %in% rownames(du)) %>% # keep only the samples of interest
  tibble::column_to_rownames(var = "sample_name") # vegan requires sample_names as row names
adonis_meta$Tile<-as.factor(adonis_meta$Tile)
adonis_meta$Tank<-as.factor(adonis_meta$Tank)
adonis_meta 

# call the permanova test and specify the model factors
adonis(dist ~ Species, strata=adonis_meta$Tank, data = adonis_meta)
adonis(dist ~ Treatment, strata=adonis_meta$Tank, data = adonis_meta)


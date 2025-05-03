.libPaths("/scratch/group/kitchen-group/Rlibs")

# install new packages
library("phyloseq")
library("decontam")
library(reshape2)
library(vegan)
library(ggplot2)
library(microshades)

# set your working directory, change "kitchens" to your specific directory name
setwd("/scratch/group/kitchen-group/MARB_689_Molecular_Ecology/class_working_directories/jcj98/project/test")

# import your files
table <- read.table(file = "feature-table.txt", sep = "\t", header = T, row.names = 1, skip = 1, comment.char = "")
OTU = otu_table(table,taxa_are_rows=TRUE)
metadata <- read.table("metadata2.tsv", sep = "\t", header = TRUE, row.names = 1, quote = "", comment.char = "")
SAMPLE <- sample_data(metadata)
taxonomy2 <- read.table(file = "taxonomy.tsv", sep = "\t", header = TRUE ,row.names = 1)
tax_split <- strsplit(as.character(taxonomy2$Taxon), ";")

# Convert to a data frame
tax_df <- do.call(rbind, lapply(tax_split, function(x) {
  # Remove prefix like "D_0__", "D_1__" and keep just the names
  sapply(x, function(y) sub(".*__", "", y))
}))

# Name columns appropriately
colnames(tax_df) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")[1:ncol(tax_df)]

# Add rownames and make it a tax_table
rownames(tax_df) <- rownames(taxonomy2)
TAX2 <- tax_table(as.matrix(tax_df))  # if using phyloseq

TREE <- read_tree("tree.nwk")

# create phyloseq object
physeq <- phyloseq(OTU, TAX2, SAMPLE,TREE)
physeq

# plot alpha diversity
plot_richness(physeq, x = "Description", measures = c("Observed", "Shannon", "Simpson")) +
  geom_boxplot() +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.text.x.bottom = element_text(angle = -90, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    strip.text = element_text(color = "black"),
    plot.title = element_text(color = "black")
  )
# create an ordination plot
ord <- ordinate(physeq,"PCoA","unifrac")
plot_ordination(physeq,ord,type="samples",color="Description")

# taxa bar plot
# Use microshades function prep_mdf to agglomerate, normalize, and melt the phyloseq object
mdf_prep <- prep_mdf(physeq, subgroup_level = "Species")

# Create a color object for the specified data
p<-plot_bar(physeq, fill = "Class", facet_grid="~description")
p=p+facet_grid(~Description, scale="free_x", drop=TRUE)
p

library(RColorBrewer)

# Create a custom palette with as many colors as you need
###Phylum
phylum_list <- unique(mdf_prep$Phylum)
n_phyla <- length(phylum_list)

set.seed(123)  # for reproducibility
phylum_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_phyla)
names(phylum_colors) <- phylum_list
PhylumPlot <- ggplot(mdf_prep, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Description, scales = "free_x") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = -90, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    plot.title = element_text(color = "black"),
    strip.text = element_text(color = "black")
  ) +
  scale_fill_manual(values = phylum_colors)

plot(PhylumPlot)
ggplotly(PhylumPlot)
# Create a species-level label by combining Genus + Species, or using Genus if Species is missing
mdf_prep$SpeciesFull <- ifelse(
  is.na(mdf_prep$Species) | mdf_prep$Species == "", 
  mdf_prep$Genus, 
  paste(mdf_prep$Genus, mdf_prep$Species, sep = " ")
)

# Generate unique species list and assign a custom color to each
species_list <- unique(mdf_prep$SpeciesFull)
n_species <- length(species_list)

set.seed(123)
species_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_species)
names(species_colors) <- species_list

###Genus-Species
install.packages("plotly")
library(plotly)

SpeciesFullPlot <- ggplot(mdf_prep, aes(x = Sample, y = Abundance, fill = SpeciesFull)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Description, scales = "free_x") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = -90, size = 6),
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 8),
    legend.position = "right"
  ) +
  scale_fill_manual(
    name = "Genus–Species",  # ← Custom legend title here
    values = species_colors
  )

ggplotly(SpeciesFullPlot)
###Genus
genus_list <- unique(mdf_prep$Genus)
n_genus <- length(genus_list)

set.seed(123)
genus_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_genus)
names(genus_colors) <- genus_list

GenusPlot <- ggplot(mdf_prep, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Description, scales = "free_x") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = -90, size = 6, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    plot.title = element_text(color = "black"),
    legend.text = element_text(size = 5),
    legend.title = element_text(color = "black"),
    legend.position = "right"
  ) +
  scale_fill_manual(values = genus_colors)

ggplotly(GenusPlot)


####Family

family_list <- unique(mdf_prep$Family)
n_family <- length(family_list)

set.seed(123)
family_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_family)
names(family_colors) <- family_list

FamilyPlot <- ggplot(mdf_prep, aes(x = Sample, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Description, scales = "free_x") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = -90, size = 6, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    plot.title = element_text(color = "black"),
    legend.text = element_text(size = 5),
    legend.title = element_text(color = "black"),
    legend.position = "right"
  ) +
  scale_fill_manual(values = family_colors)

ggplotly(FamilyPlot)

####Order

order_list <- unique(mdf_prep$Order)
n_order <- length(order_list)

set.seed(123)
order_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_order)
names(order_colors) <- order_list

OrderPlot <- ggplot(mdf_prep, aes(x = Sample, y = Abundance, fill = Order)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Description, scales = "free_x") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = -90, size = 6, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    plot.title = element_text(color = "black"),
    legend.text = element_text(size = 5),
    legend.title = element_text(color = "black"),
    legend.position = "right"
  ) +
  scale_fill_manual(values = order_colors)

ggplotly(OrderPlot)

####Class

class_list <- unique(mdf_prep$Class)
n_class <- length(class_list)

set.seed(123)
class_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_class)
names(class_colors) <- class_list

ClassPlot <- ggplot(mdf_prep, aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Description, scales = "free_x") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = -90, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    plot.title = element_text(color = "black"),
    strip.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    legend.text = element_text(size = 5),
    legend.position = "right"
  ) +
  scale_fill_manual(values = class_colors)

plot(ClassPlot)
ggplotly(ClassPlot)

###Top10
top_species <- names(
  sort(tapply(mdf_prep$Abundance, mdf_prep$SpeciesFull, sum), decreasing = TRUE)[1:10]
)

mdf_prep$SpeciesTop10 <- ifelse(
  mdf_prep$SpeciesFull %in% top_species,
  mdf_prep$SpeciesFull,
  "Other"
)

n_species10 <- length(unique(mdf_prep$SpeciesTop10))
species_list10 <- unique(mdf_prep$SpeciesTop10)

set.seed(123)
species_colors10 <- colorRampPalette(brewer.pal(12, "Set3"))(n_species10)
names(species_colors10) <- species_list10

SpeciesFullPlot10 <- ggplot(mdf_prep, aes(x = Sample, y = Abundance, fill = SpeciesTop10)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Description, scales = "free_x") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = -90, size = 6, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    plot.title = element_text(color = "black"),
    strip.text = element_text(color = "black"),
    legend.text = element_text(size = 6),
    legend.title = element_blank(),  # no title
    legend.position = "right"
  ) +
  scale_fill_manual(values = species_colors10)

plot(SpeciesFullPlot10)
ggplotly((SpeciesFullPlot10))
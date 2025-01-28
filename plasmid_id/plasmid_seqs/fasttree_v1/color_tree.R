library(ggtree)
library(dplyr)
print(getwd())
setwd('/Users/mf019/bioinformatics/longread_GWAS/plasmid_id/plasmid_seqs/fasttree_v1/')
print(getwd())

# Read in the newick tree file
tree <- read.tree("/Users/mf019/bioinformatics/longread_GWAS/plasmid_id/plasmid_seqs/fasttree_v1/fasttree_all_pf32_genes_aligned.newick")

# Read in the known types from a txt file
known_types <- read.table("list_of_plasmids.csv", header = TRUE)

# Convert the known types to a regular expression pattern
known_types_pattern <- paste(known_types$type, collapse = "|")

# Add a color column to the known_types table
known_types$color <- rainbow(nrow(known_types))

# Create a dataframe that maps each label to a color
label_colors <- data.frame(
  label = c(known_types$type, "Other"),
  color = c(known_types$color, "black")  # Use black for "Other"
)

# Convert the label column to a factor
label_colors$label <- as.factor(label_colors$label)

# Merge the label_colors dataframe with the tree data
colored_tree$data <- merge(colored_tree$data, label_colors, by = "label", all.x = TRUE)

# Create the base tree
base_tree <- ggtree(tree)

# Add a layer to color the node labels
colored_tree <- base_tree +
  geom_nodelab() +
  geom_tiplab(aes(label = label, color = color), data = label_colors)

# Save the colored tree as a plot
ggsave("colored_tree.png", plot = colored_tree, width = 25, height = 25)

# Display the colored tree
print(colored_tree)


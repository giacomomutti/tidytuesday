library(tidyverse)
library(treeio)
library(ggtree)
library(ggridges)

lemurs <- readr::read_csv('https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2021/2021-08-24/lemur_data.csv')
tax_cod <- read_tsv("data/lemur_codes.txt")
tax_cod$tiplab <- sapply(tax_cod$Latin_name, function(y) paste(unlist(strsplit(y, " "))[1:2], collapse = "_"))
tax_cod$genus <- sapply(tax_cod$Latin_name, function(y) paste(unlist(strsplit(y, " "))[1]))

og_tree <- read.newick("data/Lemur_RAxML_bestTree_total_evidence.tree")
subtree <- ape::drop.tip(og_tree, og_tree$tip.label[!(og_tree$tip.label %in% tax_cod$tiplab)])

lemurs <- lemurs %>%
  left_join(tax_cod, by=c("taxon"="Taxon")) %>% 
  filter(tiplab %in% subtree$tip.label)

tomerge <- lemurs %>% 
  group_by(tiplab) %>% 
  summarise(num=n(), 
            mn_off = mean(n_known_offspring, na.rm=TRUE), 
            age_death_mn = mean(age_at_death_y, na.rm=TRUE))
merged_tree@data$Species=c(gsub(merged_tree@phylo$tip.label,pattern = "_",replacement = " "),rep(NA,times=47-24))

p <- ggtree(subtree) %<+% tomerge


p2 <- p + geom_tiplab(offset = 0.01) + # xlim(0,0.9) +
  geom_tippoint(aes(color = mn_off, size = num)) + 
  scale_color_viridis_c(option = "D") +
  theme(legend.position="bottom")

lemur_age <- lemurs %>% 
  select(tiplab, age_at_death_y, genus) %>% 
  filter(!is.na(age_at_death_y)) %>% 
  as.data.frame()

# Generate distribution of points for each species

facet_plot(p2 + xlim_tree(0.8), panel="AGE", data=lemur_age, ggridges::geom_density_ridges, 
           mapping = aes(x=age_at_death_y, group=label, fill=genus), colour="black", lwd=0.2,
           show.legend=FALSE) +
  scale_fill_manual(values = tableau_color_pal("Classic Cyclic")(13))


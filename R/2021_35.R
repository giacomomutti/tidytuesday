library(tidyverse)
library(treeio)
library(ggtree)
library(ggridges)
library(ggdendro)
library(patchwork)


lemurs <- readr::read_csv('https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2021/2021-08-24/lemur_data.csv')
tax_cod <- read_tsv("data/2021_35/lemur_codes.txt")
tax_cod$tiplab <- sapply(tax_cod$Latin_name, function(y) paste(unlist(strsplit(y, " "))[1:2], collapse = "_"))
tax_cod$genus <- sapply(tax_cod$Latin_name, function(y) paste(unlist(strsplit(y, " "))[1]))

## Reordering tax_cod$genus
tax_cod$genus <- fct_relevel(
  tax_cod$genus,
  "Otolemur", "Galago", "Perodicticus", "Loris", "Nycticebus",
  "Daubentonia", "Cheirogaleus", "Mirza", "Microcebus", "Propithecus",
  "Varecia", "Lemur", "Hapalemur", "Eulemur"
)

og_tree <- read.newick("data/2021_35/Lemur_RAxML_bestTree_total_evidence.tree")
subtree <- ape::drop.tip(og_tree, og_tree$tip.label[!(og_tree$tip.label %in% tax_cod$tiplab)])

lemurs <- lemurs %>%
  left_join(tax_cod, by=c("taxon"="Taxon")) %>% 
  filter(tiplab %in% subtree$tip.label)

tomerge <- lemurs %>% 
  group_by(tiplab) %>% 
  summarise(num=n(), 
            mn_off = mean(n_known_offspring, na.rm=TRUE), 
            age_death_mn = mean(age_at_death_y, na.rm=TRUE))

p <- ggtree(subtree) %<+% tomerge

palette <- ggthemes::tableau_color_pal("Classic Cyclic")(13)
palette <- c(palette[1:8],"#F0183E", palette[9:13])
names(palette) <- levels(tax_cod$genus)

p$data$label <- gsub("_"," ",p$data$label)
p$data$genus_1 <- sapply(p$data$label, function(y) paste(unlist(strsplit(y, " "))[1]))

# p$data$image <- rep("", nrow(p$data))
# 
p + geom_text2(aes(label=node))
# 
# phylopic_info <- data.frame(node = c(39, 38, 30, 36, 31, 48),
#                             phylopic = c("5388a472-994a-48e1-86de-d988c6019e72",
#                                          "d6cfb28f-136e-4a20-a5ac-8eb353c7fc4a",
#                                          "f6aa4027-589b-4bf6-a603-cff6d3eb2771",
#                                          "3674b490-3fc6-4349-9cca-9d64c98e309b",
#                                          "03bbb6de-1690-4117-b459-f6d1b3ed129e",
#                                          "ae169dd2-af0c-4093-a988-21c0588274ae"),
#                             species = c("eulemurs", "lemurs", "nycticebus",
#                                         "varecia", "galagos", "microcebus"))

phylopic_info <- data.frame(node = c(39, 38, 22, 48, 7, 28, 31),
                            rowname = c("Eulemur","Hapalemur","Propithecus","Mirza","Daubentonia","Perodicticus","Otolemur"),
                            nudge = c(0.6,1.4,0.6,1.35,0.58,1.43,0.63),
                            # phylopic = c("5388a472-994a-48e1-86de-d988c6019e72",
                            #              "d6cfb28f-136e-4a20-a5ac-8eb353c7fc4a",
                            #              "f598fb39-facf-43ea-a576-1861304b2fe4",
                            #              "aceb287d-84cf-46f1-868c-4797c4ac54a8",
                            #              "0174801d-15a6-4668-bfe0-4c421fbe51e8",
                            #              "72f2f854-f3cd-4666-887c-35d5c256ab0f"),
                            image_svg = c("data/2021_35/eulemur.svg",
                                          "data/2021_35/lemur.svg",
                                          "data/2021_35/prophitecus.svg",
                                          "data/2021_35/mirza.svg",
                                          "data/2021_35/daubentonia.svg",
                                          "data/2021_35/nycticebus.svg",
                                          "data/2021_35/galago.svg"))


p2 <- p %<+% phylopic_info + geom_nodelab(aes(image=image_svg), geom="image", alpha=.5, nudge_x = -0.01) + 
  geom_tiplab(aes(image=image_svg), geom="image", alpha=.5, nudge_x = -0.06) +
  geom_tiplab(offset = 0.01, family="Lato", fontface = "italic", align=TRUE) + # xlim(0,0.9) +
  geom_tippoint(aes(color = genus_1, size = num)) +
  guides(color="none", size = guide_legend(override.aes = list(color="grey30"))) +
  scale_color_manual(values = palette) +
  # ggimage::geom_phylopic(aes(image=name), size=0.02, alpha=.6, color='steelblue') +
  labs(size = "# samples") + 
  theme(legend.position = "bottom")

lemur_age <- lemurs %>% 
  select(tiplab, age_at_death_y, genus) %>% 
  filter(!is.na(age_at_death_y)) %>% 
  as.data.frame()
lemur_age$tiplab <- gsub("_", " ", lemur_age$tiplab)

# Generate distribution of points for each species
def_plot <- facet_plot(p2 + xlim_tree(0.8),
                       panel="Age distribution",
                       data=lemur_age, 
                       ggridges::geom_density_ridges,
                       mapping = aes(x=age_at_death_y, group=label, fill=genus_1), 
                       colour="black", lwd=0.2, show.legend=FALSE) +
  scale_fill_manual(values = palette)


# hierarchical clustering
# age_levels <- unique(lemurs$age_category)

clust_df <- lemurs %>% 
  filter(age_category=="adult") %>% 
  mutate(age_at_wt_d = replace(age_at_wt_d, age_at_wt_d == 0, 1)) %>% 
  mutate(norm_weight = weight_g/age_at_wt_d) %>% 
  select(dlc_id, hybrid, sex, age_max_live_or_dead_y, 
         norm_weight, avg_daily_wt_change_g) %>% #, preg_status, n_known_offspring,) %>% 
  mutate(age_max_live_or_dead_y = (age_max_live_or_dead_y - mean(age_max_live_or_dead_y)) / sd(age_max_live_or_dead_y)) %>% 
  group_by(dlc_id) %>% 
  mutate(mn_norm_weight = mean(norm_weight, na.rm=T),
         mn_avg_daily = mean(avg_daily_wt_change_g, na.rm=T)) %>%
  filter(!is.na(mn_avg_daily)) %>% 
  select(-c(avg_daily_wt_change_g, norm_weight)) %>% 
  distinct(dlc_id, .keep_all = TRUE) %>% 
  column_to_rownames("dlc_id")

clust_df$hybrid <- as.factor(clust_df$hybrid)
clust_df$sex <- as.factor(clust_df$sex)

df_legend <- palette %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  rename("color"=".") %>% 
  mutate(y_col=seq(1, length(palette)), x_col=1) %>% 
  left_join(phylopic_info)

lemurs_unique <- lemurs %>% 
  select(dlc_id, taxon, tiplab, genus, sex) %>% 
  distinct()

dist_mat <- cluster::daisy(clust_df)
dendrogram <- as.dendrogram(hclust(dist_mat, method = 'average'))
# plot(dendrogram)
divided_dendro <- dendextend::get_subdendrograms(dendrogram, 2)

list_plots <- NULL
num <- 1
for (i in divided_dendro){
  dendrogram_data <- dendro_data(i)
  dendrogram_segments <- dendrogram_data$segments
  dendrogram_ends <- dendrogram_segments %>%
    filter(yend == 0) %>% # filter for terminal dendrogram ends
    left_join(dendrogram_data$labels, by = c("xend"="x")) %>% # .$labels contains the row names from dist_matrix (i.e., sample_name)
    rename(dlc_id = label) %>%
    left_join(lemurs_unique, by = "dlc_id")
  dendroplot <- ggplot() +
    geom_segment(data = dendrogram_segments, 
                 aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_segment(data = dendrogram_ends,
                 aes(x=x, y=y.x, xend=xend, yend=yend, color = genus)) +
    geom_point(data = dendrogram_ends,
                 aes(x=x, y=y.x, color = genus)) +
    scale_color_manual(values = palette) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    coord_flip() +
    theme_void() + 
    theme(legend.position = "none") 
  
  list_plots[[num]] <- dendroplot
  num <- num + 1
  
  max_val <- max(dendrogram_ends$xend)
  min_val <- min(dendrogram_ends$xend)
  
  a <- dendrogram_ends %>% 
    left_join(df_legend, by = c("genus"="rowname")) %>% 
    mutate(new_x = (max_val - min_val) * (.$y_col - 1) / (14 - 1) + min_val) %>% 
    # mutate(new_x = (14-1) * ((dendrogram_ends$xend - min(dendrogram_ends$xend)) / (max(dendrogram_ends$xend) - min(dendrogram_ends$xend))) + 1) %>% 
    ggplot(aes(color=color)) + 
    geom_segment(aes(x = yend, y = xend, yend = jitter(new_x), xend = x_col), alpha=0.4) + 
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) + 
    scale_color_identity() +
    theme_void()
  list_plots[[num]] <- a
  num <- num + 1
}

legend_plot <- df_legend %>% 
  ggplot(aes(x_col, y_col, fill=color)) + 
  geom_tile() + 
  geom_text(aes(label=rowname), family="BebasNeue", fontface="bold", color = "white", size = 9) + 
  geom_image(aes(image=image_svg, x=nudge), color="white",by = "height", size = 0.08) +
  theme_void() + 
  theme(legend.position = "none") + 
  scale_fill_identity() + 
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

a <- list_plots[[1]] + scale_y_reverse(expand=c(0,0)) + ggtitle("Males") + theme(plot.title = element_text(family="BebasNeue", size = 20, hjust = 1)) 
b <- list_plots[[2]] 
c <- legend_plot
d <- list_plots[[4]] + scale_x_reverse(expand=c(0,0)) 
e <- list_plots[[3]] + ggtitle("Females") + theme(plot.title = element_text(family="BebasNeue", size = 20))

def_dendro <- a | b | c | d | e
def_dendro <- def_dendro + patchwork::plot_layout(widths = c(2,1.2,1.2,1.2,2))
# all_plots <- def_plot / def_dendro
ggsave("plots/2021_35_tree_age.png", def_plot, height = 10, width = 14)
ggsave("plots/2021_35_dendro_sex.png", def_dendro, height = 7, width = 14)
# change font add male and female labels and captio


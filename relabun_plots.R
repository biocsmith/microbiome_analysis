### Relative abundance plots using the microshades package: https://github.com/KarstensLab/microshades

library(phyloseq) ### phyloseq object
library(tidyverse)
library(ggplot2) ### plotting
library(microshades) ### relative abundance plots
library(cowplot)

phylo <- readRDS("noncontam2.RDS") ### import phyloseq object

### Order phyla by most abundant

top_phy <- psmelt(phylo %>% microbiome::transform("compositional")) %>% group_by(Phylum) %>% summarise(mean_taxa = mean(Abundance)) %>% arrange(desc(mean_taxa))

### aggregate taxa by genus, rename 'NA' genera to higher taxanomic classification

phylo.gen <- phyloseq::tax_glom(phylo, "Genus", NArm = FALSE)
phylo.gen.rel <- microbiome::transform(phylo.gen, "compositional")
ps <- psmelt(phylo.gen.rel)

ps[ps == ""] <- NA
ps <- ps %>%
    mutate(Phylum = coalesce(Phylum, Kingdom),
           Class = coalesce(Class, Phylum),
           Order = coalesce(Order, Class),
           Family = coalesce(Family, Order),
           Genus = coalesce(Genus, Family))

ps$Phylum <- factor(ps$Phylum, levels = rev(top_phy$Phylum))

### selected_groups is a vector of the top (up to 5) abundant phyla.

color_obs <- create_color_dfs(ps, selected_groups = top_phy$Phylum[1:5], group_level = "Phylum", subgroup_level = "Genus", cvd = TRUE)
mdf <- color_obs$mdf
cdf <- color_obs$cdf
mdf$Phylum <- factor(mdf$Phylum, levels = rev(top_phy$Phylum))
cdf$Phyum <- factor(cdf$Top_Phylum, levels = rev(c(top_phy$Phylum[1:5], "Other")))
# plot <- plot_microshades(mdf, cdf)

GP_legend <- custom_legend(mdf, cdf, legend_key_size = 0.5, legend_text_size = 8)

group_label <- "Phylum Genus" 
x <- "Sample"
y <- "Abundance"

plot <- mdf %>% 
    ggplot(aes_string(x = x, y = y), fill = group_label, width = 1) + 
    scale_fill_manual(name = group_label,
                      values = cdf$hex,
                      breaks = cdf$group) +
    scale_colour_manual(name = group_label,
                        values = cdf$hex,
                        breaks = cdf$group) +
    geom_bar(aes(fill=group), stat="identity", position="stack") + theme_void() + scale_y_continuous(labels = scales::percent, expand = expansion(0)) + 
  theme(legend.position = "none", 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_blank(), 
        panel.spacing = unit(0.2, "lines"), 
        strip.text.x = element_text(size = 10, angle = 90), 
        strip.clip = "off", 
        plot.margin = margin(0.5,0,0.2,0.2, "cm"))

plot_grid(plot, GP_legend, rel_widths = c(1, 0.4), align = "v", axis = "t")


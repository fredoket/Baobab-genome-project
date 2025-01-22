
#Load the libraries
pacman::p_load(tidyverse, janitor,ggpubr, ggtext,
               dplyr, viridis, forcats,ape,
               ggtree, rphylopic, ggimage,
               palmerpenguins, plyr, ggdist,
               colorspace, ragg, gghalves,
               gtsummary)

# Load data
df_dna <- read_csv("data/DNA_extraction.csv") |>
  select(-Sample_ID, -ABR) |>
  filter(PVP != "PST") |>
  filter(Temp != "Cold") |>
  mutate(tissue = case_when(
    Tissue == "Cotyledon leaf" ~ "Cotyledon leaves",
    Tissue == "Mature tree young foliage leaf" ~ "Mature trees leaves",
    Tissue == "Seedling  young foliage leaf" ~ "Seedlings leaves",
    Tissue == "Seed" ~ "Seeds"
  )) |>
  clean_names()


df_dfn <- df_dna |>
  mutate('Purity (260/280)' = case_when(
    purity_260_a280 < 1.78 ~ "Low",
    purity_260_a280 >= 1.78 & purity_260_a280 <= 1.95 ~ "Within",
    purity_260_a280 > 1.95 ~ "High"
  )) |>
  dplyr::rename("Concentration (ng/µL)" = "concentration_ng_m_l")


df_dfn <- df_dfn |>
  mutate('Purity (260/230)' = case_when(
    purity_260_a230 < 2.0 ~ "Low",
    purity_260_a230 >= 2.0 & purity_260_a230 <= 2.2 ~ "Within",
    purity_260_a230 > 2.2 ~ "High"
)) |>
  select(tissue, `Purity (260/280)`, `Purity (260/230)`,
         `Concentration (ng/µL)`, est_qty) |>
  dplyr::rename("Estimated concentration (ng/µL)" = "est_qty")

df_tbl <- df_dfn |>
  tbl_summary(by = tissue) |>
  add_p() |>
  modify_header(label ~ "**Variables**") |>
  modify_spanning_header(c("stat_1", "stat_2", "stat_3", "stat_4") ~ "**Tissue Source**") |>
  bold_labels() |>
  flextable::save_as_docx(path = "table1.docx")
  

df_graph <- df_dna |>
  group_by(abr, category_280) |>
  dplyr::summarise(Count = n()) |>
  dplyr::mutate(Propotion = Count*100/sum(Count))

p1 <- ggplot(df_graph, aes(x = abr, y = Propotion,
                           fill = category_280,
                     width = 0.7)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c('#55e2e2',
                               'grey',
                               '#34495e'))+
  theme_ggdist()+
    theme(axis.title = element_text(size = 13,
      face = "bold"),
axis.text = element_text(size = 10),
axis.text.x = element_text(size = 14),
axis.title.x = element_text(size = 16, face = "bold"),
axis.title.y = element_text(size = 16, face = "bold"),
axis.text.y = element_text(size = 14),

plot.title = element_text(size = 14,
      hjust = 0.5),
plot.subtitle = element_text(size = 12,
         hjust = 0.5))+
xlab("DNA Source") +
  guides(fill=guide_legend(title="Category"))+
  ggtitle("(a)")

df_graph1 <- df_dna |>
  group_by(abr, category_260) |>
  dplyr::summarise(Count = n()) |>
  dplyr::mutate(Propotion = Count*100/sum(Count))

p2 <- ggplot(df_graph1, aes(x = abr, y = Propotion,
                            fill = category_260,
                            width = 0.7)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c('#55e2e2',
                               'grey',
                               '#34495e'))+
  xlab("DNA Source")+
  theme_ggdist() +
  theme(axis.title = element_text(size = 13,
                                  face = "bold"),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 14,
                                  hjust = 0.5),
        plot.subtitle = element_text(size = 12,
                                     hjust = 0.5))+
  xlab("DNA Source") +
  ggtitle("(b)")

ggarrange(p1, p2, ncol = 2, common.legend = TRUE,
          legend = "bottom")


df_raincloud <- ddply(df_dna, c("tissue",
                                "abr"),
                      summarize,
                      Mean = mean(est_qty),
                      SD = sd(est_qty))

#Raincloud plot for the quality values at absorbence A260/A280
df_rcloud <- df_dna |>
  group_by(tissue) |>
  filter(!is.na(est_qty))

#Plot  
ggplot(df_rcloud, aes(x = fct_rev(tissue), y = od_280)) + 
  ggdist::stat_halfeye(
    aes(#fill = tissue
        ),
    adjust = .5, 
    width = .75, 
    .width = 0,
    justification = -.4, 
    point_color = NA) + 
  geom_boxplot(
    aes(#fill = tissue
        ),
    width = .2,
    position = position_nudge(x = 0.09),
    outlier.shape = NA
  )+
  geom_point(position = position_nudge(x = -0.10),
             shape = "|", size = 3)+
  stat_summary(
    geom = "text",
    fun = "median",
    aes(label = round(..y.., 3),
        ),
    family = "Roboto Mono",
    fontface = "bold",
    size = 3.5,
    vjust = -3.9
  )+
  coord_flip(xlim = c(0, NA), clip = "off") +
  scale_y_continuous(
   limits = c(0.8, NA)
  )+
  theme_ggdist()+
  ylab("Absorbance A260/A280")+
  xlab("DNA Source")+
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14))


  
# Mean error barplot for the quantity of DNA per tissue
ggplot(df_raincloud,
       aes(x = abr, y = Mean))+
  geom_bar(stat = "identity", width = 0.5)+
  geom_errorbar(aes(ymin = Mean - SD,
                    ymax = Mean + SD),
                width = 0.2)+
  theme_classic()+
  theme(axis.title = element_text(size = 16,
                                  face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  geom_text(aes(label = round(Mean, 2),
                size = 4,
                vjust = -4))+
  xlab("DNA Source")




df_coty <- df_dna |>
  filter(abr == "CTL")

df_mat <- df_dna |>
  filter(abr == "MTYFL")

df_sap <- df_dna |>
  filter(abr == "SYFL")

df_seed <- df_dna |>
  filter(abr == "SDS")

df_low <- df_dna |>
  filter(category_260 == "Low") |>
  group_by(category_260) |>
  summarise(n = n())

par(mfrow = c(2, 2))

par(bty = "l") 

# Draw box plots for each group in the data frame
boxplot(od_280 ~ category_280,
        data = df_coty,
        main = "A",
        xlab = substitute(paste(bold('Cotyledon Leaves'))),
        ylab = "A230/A280",
        boxwex=.35,
        par(bty = "l"))

boxplot(od_280 ~ category_280,
        data = df_mat,
        main = "B",
        xlab = substitute(paste(bold('Mature Trees Foliage Leaves'))),
        ylab = "A230/A280",
        boxwex=.35,
        par(bty = "l"))

boxplot(od_280 ~ category_280,
        data = df_sap,
        main = "C",
        xlab = substitute(paste(bold('Seedlings Foliage Leaves'))),
        ylab = "A230/A280",
        boxwex=.35,
        par(bty = "l"))

boxplot(od_280 ~ category_280,
        data = df_seed,
        main = "D",
        xlab = substitute(paste(bold('Seeds'))),
        ylab = "A230/A280",
        boxwex=.35,
        par(bty = "l"))

#2. CHLOROPLAST GENOME PHYLOGENETIC TREE
df_tree <- read.tree("data/itol_newik.txt")

outgroup <- "Liquidambar_formosana"

df_rooted <- root(df_tree,
                  outgroup = outgroup,
                  resolve.root = TRUE)

df_rooted$tip.label <- gsub("_", " ", df_rooted$tip.label)

highlight_species <- "Adansonia digitata"

df_rooted |>
  ggtree(branch.length = "none")+
  theme_tree2()+
  xlim(NA, 15)+
  geom_tiplab(geom = "text", aes(label = label),
              offset = 0.2, fontface ="italic",)+
  geom_nodepoint(pch = 1)+
  geom_tippoint(pch = 23)+
  geom_hilight(node=23, fill="forestgreen", alpha=0.3)+
  geom_hilight(node=28, fill="forestgreen", alpha=0.25)+
  geom_hilight(node=26, fill="forestgreen", alpha=0.15)+
  geom_hilight(node=32, fill="forestgreen", alpha=0.1)+
  geom_cladelab(node = 23, label = "Malvaceae",
                color="red", offset= 3.3, align=TRUE)+
  geom_cladelab(node = 28, label = "Caryophyllaceae",
                color="red", offset= 3.3, align=TRUE)+
  geom_cladelab(node = 4, label = "Platanaceae",
                color="red", offset= 3.3, align=TRUE)+
  geom_cladelab(node = 3, label = "Berberidaceae",
                color="red", offset= 3.3, align=TRUE)+
  geom_cladelab(node = 5, label = "Papaveraceae",
                color="red", offset= 3.3, align=TRUE)+
  geom_cladelab(node = 12, label = "Crassulaceae",
                color="red", offset= 3.3, align=TRUE)+
  geom_cladelab(node = 11, label = "Penthoraceae",
                color="red", offset= 3.3, align=TRUE)+
  geom_cladelab(node = 10, label = "Altingiaceae",
                color="red", offset= 3.3, align=TRUE)#+
  #geom_text(aes(label=node), hjust=-.3)


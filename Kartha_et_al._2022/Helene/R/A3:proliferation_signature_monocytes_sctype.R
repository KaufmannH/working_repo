# Are there monocytes with a proliferation signature?
# used sctype

library(Seurat)
library(ggplot2)
library(dplyr)


palette <- c("#70916B", "#A4C089", "#D6EFC3", "#F5F5F5", "#C7EAE5" ,"#5AB4AC", 
"#123230",  "#556967",  "#B8C7C7",  "#E0D5C8",  "#D0BFAD",  "#AE9982", "#BA613E", "#843915", "#471339", "#3B052C")


pbmc <- readRDS("Helene/data/pmbc_seurat_object_sctype.rds")
monocytes <- subset(pbmc, idents = "Non-classical monocytes")
head(monocytes)


#  cell cycle genes
cc.genes <- Seurat::cc.genes 
cc.genes$s.genes

# check proliferation/cell cycle genes
monocytes <- CellCycleScoring(
  object = monocytes,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  set.ident = TRUE  
)

monocytes$CellCyclePhase <- monocytes$Phase
table(monocytes$CellCyclePhase)

# percentages
percentage_data <- monocytes@meta.data %>%
  group_by(Condition, CellCyclePhase) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Condition) %>%
  mutate(
    Percentage = (Count / sum(Count)) * 100,
    ConditionCount = sum(Count)) %>%
  ungroup()
percentage_data

# get cell numbers
numbers_condition <- percentage_data %>%
  distinct(Condition, ConditionCount) 
numbers_condition

# bar plot 
barplot_stim <- ggplot(percentage_data, aes(x = Condition, y = Percentage, fill = CellCyclePhase)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = palette[1:3]) +
  geom_text(
    data = numbers_condition, 
    inherit.aes = FALSE, 
    aes(x = Condition, y = 105, label = ConditionCount), 
    size = 4) +
  theme_classic() +
  labs(
    x = "Condition",
    y = "Percentage of Cells"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7), 
    plot.margin = margin(20, 20, 40, 20)  
  ) 
ggsave("Helene/plots/barplot_proliferation_sctype.png", plot = barplot_stim , width = 8, height = 8, dpi = 300)



# genes for proliferation extra
# Ki-67 = MKI67(nuclear division), PCNA (S), topoisomerase II alpha = TOP2A (replication)


# TODO: do stats




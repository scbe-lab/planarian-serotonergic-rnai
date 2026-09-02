library(readxl)
library(ggplot2)


score <- "/mnt/sda/elena/sero_manuscript/GO_sero/score_GOs_sero.xlsx"
sheet_score <- "all"

colors <- "/mnt/sda/elena/sero_manuscript/GO_sero/leiden_3_colors.xlsx"
sheet_color <- "leiden_3_colors"


scored_sero_go <- readxl::read_excel(score, sheet = sheet_score)
scored_sero_go <- scored_sero_go[rev(order(scored_sero_go$score)),]
scored_sero_go$names <- factor(scored_sero_go$names, levels = unique(scored_sero_go$names))

color_sero <- readxl::read_excel(colors, sheet = sheet_color)


scored_sero_go <- merge(scored_sero_go, color_sero, by = "names")


# Plot using ggplot2
ggplot(scored_sero_go, aes(y = names, x = score)) +
  geom_bar(stat = "identity", fill = scored_sero_go$Color) +
  labs(title = "Bar Plot of Scores of Genes with Annotated Serotonergic GOs",
       x = "Score",
       y = "Names"
       ) +
  theme_minimal()+
  scale_y_discrete(limits=rev)



library(ggplot2)

setwd("re_gecip/machine_learning/AlexBlakes/near_splice/write_up/paper/code/branchpoints/plotting/")

# Get summary stats for splicing variants
df = read.csv("../stats/near_splice_branch_stats_combined.tsv", sep="\t")

# Get summary stats for coding variants from the near-splice MAPS output.
coding = read.csv("~/re_gecip/machine_learning/AlexBlakes/near_splice/write_up/paper/code/near_splice/stats/maps_output.tsv", sep="\t")
coding = subset(coding, xlab %in% c("Synonymous","Missense","Nonsense"))

# Near-splice phyloP, SpliceAI, and MAPS plot.
p = ggplot(df, aes(site, mean, colour = subset)) +
  geom_pointrange(
    aes(ymin = mean - 1.96 * sem,
        ymax = mean + 1.96 *sem
        ),
    size=.6,
    position = position_dodge2(0.5)
    ) +
  scale_x_continuous(breaks = c(-25:10)) +
  scale_color_manual(labels = c("Score >=0.85", "All"),
                     values = c("steelblue2", "royalblue4")) +
  facet_grid(cols = vars(factor(region, levels = c("Branchpoint","Acceptor","Donor"))),
             rows = vars(factor(score, levels = c("phyloP", "SpliceAI","MAPS"))),
             scale = "free",
             space = "free_x",
             switch = "both") +
  theme_bw() +
  theme(
    axis.title = element_blank(), 
    axis.text.x = element_text(size = 14, vjust = 0.5, angle = 90, hjust = 1), 
    axis.text.y = element_text(size = 16),
    strip.text = element_text(size = 18),
    strip.placement = "outside",
    legend.title = element_blank(),
    legend.position = c(0.08,0.90),
    legend.background = element_blank(),
    legend.text = element_text(size = 11),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.spacing.y = unit(1,"lines")
  )
p

# Coding MAPS plot
p1 = ggplot(coding, aes(factor(xlab, levels = c("Synonymous","Missense","Nonsense")),
                        maps)) +
  geom_pointrange(
    aes(ymin = maps - 1.96 * se,
        ymax = maps + 1.96 *se
    ),
    size=1,
    position = position_dodge2(0.5),
    colour = "royalblue4"
  ) +
  scale_y_continuous(limits = c(-0.020,0.152), 
                     position="right") +
  theme_bw() +
  theme(
    axis.title = element_blank(), 
    axis.text.x = element_text(size = 20, vjust = 0.5, angle = 90, hjust = 1), 
    axis.text.y = element_text(size = 20),
    strip.text = element_text(size = 24),
    strip.placement = "outside",
    strip.background = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.08,0.90),
    legend.background = element_blank(),
    legend.text = element_text(size = 16),
    #panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.spacing.y = unit(1,"lines")
  )
p1

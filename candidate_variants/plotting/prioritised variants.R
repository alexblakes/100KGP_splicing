library(ggplot2)
library(stringr)

setwd("~/re_gecip/machine_learning/AlexBlakes/near_splice/write_up/paper/code/candidate_variants/plotting")

# Get summary stats for splicing variants
df = read.csv("../outputs/candidate_dnms.tsv", sep="\t")
df$tier[is.na(df$tier)] <- "Not tiered"
df$tier <- factor(df$tier, levels = c("Not tiered", "3", "1"))
df$eq_outcome <- factor(df$eq_outcome, levels = c("No data",
                                                  "Case not solved",
                                                  "Case solved, variant unknown",
                                                  "Case solved, different variant",
                                                  "Case solved, same variant"
))

# Get summary stats for coding variants from the near-splice MAPS output.

# SpliceAI scores
spliceai_plot = ggplot(df, aes(site)) +
  geom_dotplot(aes(fill = DS_any, group = DS_any),
               binwidth = 1, 
               stackgroups = TRUE,
               method = "histodot",
               binpositions = "all",
               dotsize = 1,
               stackratio = 1) +
  scale_x_continuous(breaks = seq(-25,10,5)) +
  scale_y_continuous(NULL, breaks = NULL) +
  scale_fill_gradient2(low = "white", 
                       high = "firebrick1",
                       na.value="grey80",
                       limits = c(0,1),
                       breaks = c(0,1),
                       "SpliceAI"
  ) +
  facet_grid(cols = vars(factor(str_to_title(region), levels = c("Branch","Acceptor","Donor"))),
             scale = "free",
             space = "free",
             switch = "both") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.line.x = element_line(),
    axis.text.x = element_text(size = 25, vjust = 0.5, angle = 90, hjust = 1), 
    axis.text.y = element_blank(),
    strip.text = element_text(size = 30),
    strip.placement = "outside",
    strip.background = element_blank(),
    legend.title = element_text(size = 25),
    legend.position = c(0.1, 0.6),
    legend.background = element_blank(),
    legend.text = element_text(size = 25),
    panel.grid = element_blank(),
    panel.spacing.y = unit(1,"lines"),
    panel.border = element_blank()
  )

spliceai_plot

# Tiering scores
tiering_plot = ggplot(df, aes(site)) +
  geom_dotplot(aes(fill = tier, group = tier),
               binwidth = 1, 
               stackgroups = TRUE,
               method = "histodot",
               binpositions = "all",
               dotsize = 1,
               stackratio = 1
  ) +
  scale_x_continuous(breaks = seq(-25,10,5)) +
  scale_y_continuous(NULL, breaks = NULL) +
  scale_fill_manual(values = c("Not tiered" = "gray80",
                               "3" = "olivedrab3",
                               "1" = "forestgreen"),
                    "Tier") +
  facet_grid(cols = vars(factor(str_to_title(region), levels = c("Branch","Acceptor","Donor"))),
             scale = "free",
             space = "free",
             switch = "both") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.line.x = element_line(),
    axis.text.x = element_text(size = 25, vjust = 0.5, angle = 90, hjust = 1), 
    axis.text.y = element_blank(),
    strip.text = element_text(size = 30),
    strip.placement = "outside",
    strip.background = element_blank(),
    legend.title = element_text(size = 25),
    legend.position = c(0.1, 0.6),
    legend.background = element_blank(),
    legend.text = element_text(size = 25),
    panel.grid = element_blank(),
    panel.spacing.y = unit(1,"lines"),
    panel.border = element_blank()
  ) +
  guides(fill = guide_legend(reverse=TRUE,
                             override.aes = list(size = 15))
  )
tiering_plot

# EQ outcomes
outcome_plot = ggplot(df, aes(site)) +
  geom_dotplot(aes(fill = eq_outcome, group = eq_outcome),
               binwidth = 1, 
               stackgroups = TRUE,
               method = "histodot",
               binpositions = "all",
               dotsize = 1,
               stackratio = 1
  ) +
  scale_x_continuous(breaks = seq(-25,10,5)) +
  scale_y_continuous(NULL, breaks = NULL) +
  scale_fill_manual(values = c("Case solved, same variant" = "forestgreen",
                               "Case solved, different variant" = "royalblue",
                               "Case solved, variant unknown" = "gray40",
                               "Case not solved" = "gray80",
                               "No data" = "white"),
                    "Outcome") +
  facet_grid(cols = vars(factor(str_to_title(region), levels = c("Branch","Acceptor","Donor"))),
             scale = "free",
             space = "free",
             switch = "both") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.line.x = element_line(),
    axis.text.x = element_text(size = 25, vjust = 0.5, angle = 90, hjust = 1), 
    axis.text.y = element_blank(),
    strip.text = element_text(size = 30),
    strip.placement = "outside",
    strip.background = element_blank(),
    legend.title = element_text(size = 25),
    legend.position = c(0.2, 0.6),
    legend.background = element_blank(),
    legend.text = element_text(size = 25),
    panel.grid = element_blank(),
    panel.spacing.y = unit(1,"lines"),
    panel.border = element_blank()
  ) +
  guides(fill = guide_legend(reverse=TRUE,
                             override.aes = list(size = 15))
  )
outcome_plot

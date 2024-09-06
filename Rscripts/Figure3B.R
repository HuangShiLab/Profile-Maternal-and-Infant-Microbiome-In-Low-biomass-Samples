#######################################
# Monther-infant Study
#
# Figure 3B
#
# Author: HOU Shuwen
#######################################


# milk: cornbeans
{pcoa_results <- ordinate(milk, method = "PCoA", distance = distance_milk)
pcoa_plot <- plot_ordination(milk, pcoa_results, color = 'cornbeans') + 
  geom_point(size = 2, alpha = 0.6) +
  theme_minimal() + stat_ellipse() +
  ggtitle("Milk-CornBeans")
pcoa_plot}

# meconium: delivery way PCoA
pcoa_results <- ordinate(meconium, method = "PCoA", distance = distance_meconium)
pcoa_plot <- plot_ordination(meconium, pcoa_results, color = 'delivery') + 
  geom_point(size = 2, alpha = 0.8) +
  stat_ellipse() +
  scale_color_manual(
    values = c("cesarean" = "skyblue", "eutocia" = "palegreen3", "forceps" = "grey"),
    labels = c("cesarean" = "Cesarean", "eutocia" = "Vaginal", "forceps" = "Forceps")
  ) +
  labs(color = "Delivery Mode") + 
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, color = "black"))
pcoa_plot
ggsave("./plots/3B.png", pcoa_plot, width = 4.2, height = 2.5)

# meconium: other covariates PCoA
{pcoa_plot <- plot_ordination(meconium, pcoa_results, color = 'pregnancy') + 
  geom_point(size = 3) +
  theme_minimal() + stat_ellipse() +
  ggtitle("Meconium-pregnancy")
pcoa_plot
ggsave("./plots/3B.png", pcoa_plot, width = 4, height = 3)

pcoa_results <- ordinate(meconium, method = "PCoA", distance = distance_meconium)
pcoa_plot <- plot_ordination(meconium, pcoa_results, color = 'education') + 
  geom_point(size = 2, alpha = 0.6) +
  theme_minimal() + stat_ellipse() +
  ggtitle("Meconium-education")
pcoa_plot

pcoa_results <- ordinate(meconium, method = "PCoA", distance = distance_meconium)
pcoa_plot <- plot_ordination(meconium, pcoa_results, color = 'IB1B_quan') + 
  geom_point(size = 2, alpha = 0.6) + stat_ellipse() +
  theme_minimal() +
  ggtitle("Meconium-IL1B")
pcoa_plot
}





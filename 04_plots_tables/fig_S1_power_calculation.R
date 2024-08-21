###############################################################
# Power calculation sourced from https://towardsdatascience.com/how-to-calculate-statistical-power-for-your-meta-analysis-e108ee586ae8
# The above is based on https://journals.sagepub.com/doi/full/10.3102/1076998609346961
source("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/code/plot_themes.R")
###############################################################


###############################################################
# Heterogeniety (".33" for small, "1" for moderate, & "3" for large)
# n1 = average number of samples in discovery for cases, n2 for controls
calculate_power <- function(eff_size, heterogeneity, n1 = 23, n2 = 23, num_datasets = 7, p = 0.05)
{
  sigma <- (1 + heterogeneity) * (((n1 + n2) / (n1 * n2)) + ((eff_size^2) / (2 * (n1 + n2)))) / num_datasets
  lambda <- eff_size / sqrt(sigma)
  zval <- qnorm(1 - p / 2, 0, 1)
  power <- 1 - pnorm(zval - lambda) + pnorm(-zval - lambda)
  return(power)
}
###############################################################


###############################################################
eff_size_list <- c(0:200) / 200
combined_eff_power_frame <- data.frame(eff_size = rep(eff_size_list, 4),
  heterogeneity = c(rep(0, length(eff_size_list)), rep(0.33, length(eff_size_list)), rep(1, length(eff_size_list)), rep(3, length(eff_size_list))))
computed_power <- apply(combined_eff_power_frame, 1, function(x) calculate_power(x[1], x[2]))
combined_eff_power_frame$power <- computed_power
combined_eff_power_frame$heterogeneity_label <- c(rep("Fixed effect (0.315)", length(eff_size_list)), rep("RE low heterogeneity (0.36)", length(eff_size_list)), rep("RE moderate heterogeneity (0.445)", length(eff_size_list)), rep("RE high heterogeneity (0.635)", length(eff_size_list)))
###############################################################


###############################################################
vertical_lines <- data.frame(x1 = c(0.3148, 0.3602, 0.44539, 0.63549),
  y1 = c(0, 0, 0, 0),
  x2 = c(0.3148, 0.3602, 0.44539, 0.63549),
  y2 = c(0.8, 0.8, 0.8, 0.8),
  heterogeneity_label = c("Fixed effect (0.315)", "RE low heterogeneity (0.36)", "RE moderate heterogeneity (0.445)", "RE high heterogeneity (0.635)"))
###############################################################


###############################################################
combined_eff_power_frame$heterogeneity_label <- ordered(combined_eff_power_frame$heterogeneity_label, levels = c("Fixed effect (0.315)", "RE low heterogeneity (0.36)", "RE moderate heterogeneity (0.445)", "RE high heterogeneity (0.635)"))
vertical_lines$heterogeneity_label <- ordered(vertical_lines$heterogeneity_label, levels = c("Fixed effect (0.315)", "RE low heterogeneity (0.36)", "RE moderate heterogeneity (0.445)", "RE high heterogeneity (0.635)"))
power_curve <- ggplot(combined_eff_power_frame, aes(x = eff_size, y = power, colour = heterogeneity_label)) +
  xlab("Effect size") + ylab("Power") +
  geom_point(size = point_size) +
  geom_line(linewidth = line_size_text_repel_mm) +
  geom_segment(x = 0, xend = 0.636, y = 0.8, yend = 0.8, linewidth = line_size_text_repel_mm, colour = black_text_colour) +
  geom_segment(data = vertical_lines, aes(x = x1, y = y1, xend = x2, yend = y2, colour = heterogeneity_label), linewidth = line_size_text_repel_mm) +
  theme(panel.grid.major = element_line(linewidth = line_size_text_repel_mm), legend.position = c(0.65, 0.25), legend.background = element_rect(colour = "transparent", fill = "#FEFEFEDD")) +
  labs(colour = "Effect size with 80% power\nat p−value = 0.05") + guides(colour = guide_legend(ncol = 1, override.aes = list(size = text_size_labels * 6 / 9))) +
  scale_colour_manual(values = c("#1965B0", "#4EB265", "#F6C141", "#DC050C"))
ggsave("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/figures/fig_S1_power_calculation.pdf", power_curve, width = 4, height = 4, device = cairo_pdf)
###############################################################

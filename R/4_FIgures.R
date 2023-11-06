produceFigure <- function(dataTib, path, cycle) {
  # Set the storage folder
  if (!dir.exists(path)) {
    dir.create(path)
  }
  # Define aesthetics
  orange <- "#F14124"
  black <- "#000000"
  grey <- "#525252"
  blue <- "#124AAD"
  textSize <- 25
  themeElement <- theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = textSize),
    axis.text.y = element_text(size = textSize),
    axis.text.x = element_text(size = textSize, angle = 0),
    axis.title = element_text(size = textSize, vjust = 1),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.ticks.length = unit(.4, "cm")
  )
  colors <- c("bHP" = orange, "BNP" = blue, "log GDP" = black)

  # Produce the figure
  y_name <- ifelse(cycle, "log GDP", "%-Deviation")
  figure_1 <- dataTib %>%
    ggplot(aes(x = Date)) +
    geom_line(aes(y = bHP, color = "bHP"), linewidth = .75) +
    geom_line(aes(y = BND, color = "BNP"), linewidth = .75) +
    scale_x_date(date_breaks = "5 year", date_labels = c("%Y"), name = "") +
    scale_y_continuous(name = y_name) +
    scale_color_manual(values = colors) +
    theme_bw() +
    themeElement

  if (cycle) {
    figure_2 <- figure_1 +
      geom_hline(yintercept = 0, linewidth = .75)
    fileName <- "/Cycle"
  } else {
    figure_2 <- figure_1 +
      geom_line(aes(y = logGDP, color = "log GDP"), linewidth = .75)
    fileName <- "/Trend"
  }

  ggsave(figure_2, filename = paste0(path, fileName, ".pdf"), height = 8, width = 14)
}

produceGDPFigure <- function(dataTib, path){
  # Define aesthetics
  orange <- "#F14124"
  black <- "#000000"
  grey <- "#525252"
  blue <- "#124AAD"
  themeElement <- produceThemeElement(textSize = 25)
  colors <- c("GDP" = orange, "log GDP" = black)
  
  # Produce the figure
  
  max_first  <- max(dataTib$logGDP)   # Specify max of first y axis
  max_second <- max(dataTib$GDP) # Specify max of second y axis
  min_first  <- min(dataTib$logGDP)   # Specify min of first y axis
  min_second <- min(dataTib$GDP) # Specify min of second y axis
  
  # scale and shift variables calculated based on desired mins and maxes
  scale1 = (max_second - min_second)/(max_first - min_first)

  # Function to scale secondary axis
  scale_function <- function(x, scale, shift){
    return ((x)*scale - shift)
  }
  
  # Function to scale secondary variable values
  inv_scale_function <- function(x, scale, shift){
    return ((x + shift)/scale)
  }
  
  scale <- scale1 * .8
  
  figure <- dataTib %>%
    ggplot(aes(x = Date)) +
    geom_line(aes(y = inv_scale_function(GDP / 1000, scale / 1000, 0) + 1090, color = "GDP"), linewidth = .75) +
    geom_line(aes(y = logGDP, color = "log GDP"), linewidth = .75) +
    scale_x_date(date_breaks = "5 year", date_labels = c("%Y"), name = "") +
    scale_y_continuous(name = "log GDP",
                       sec.axis = sec_axis(~ scale_function(., scale / 1000, 0) + 1090, name = "GDP")) +
    scale_color_manual(values = colors) +
    theme_bw() +
    themeElement
  
  ggsave(figure, filename = paste0(path, "/GDP_figure.jpg"), height = 8, width = 14)
}


produceThemeElement <- function(textSize){
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
  return(themeElement)
}


produceTCFigure <- function(dataTib, path, cycle) {
  # Set the storage folder
  if (!dir.exists(path)) {
    dir.create(path)
  }
  # Define aesthetics
  orange <- "#F14124"
  black <- "#000000"
  grey <- "#525252"
  blue <- "#124AAD"
  themeElement <- produceThemeElement(textSize = 25)
  colors <- c("bHP" = orange, "BNP" = blue, "log GDP" = black)

  # Prepare the data
  fileName <- ifelse(cycle, "cylce", "trend")
  figureTib <- mutate(dataTib, 
                           "bHP" = bHP_result[,fileName],
                           BND = BND_results[,fileName])
  
  # Produce the figure
  y_name <- ifelse(cycle, "log GDP", "%-Deviation")
  figure_1 <- figureTib %>%
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
  } else {
    figure_2 <- figure_1 +
      geom_line(aes(y = logGDP, color = "log GDP"), linewidth = .75)
  }

  ggsave(figure_2, filename = paste0(path, "/", fileName, ".pdf"), height = 8, width = 14)
}

#' Custom ggplot2 Theme and Colors for Tidy-Style Plots
#'
#' Provides a consistent theme and color palette for publication-ready plots
#' 
#' @param base_size Base font size (default: 10)
#' @param grid_lines Add horizontal grid lines? ("none", "major", "minor", "both")
#' @param grid_color Color for grid lines (default: "gray90")
#' @return A ggplot2 theme object
#' @export
#' 
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, wt)) +
#'   geom_point() +
#'   maxwell_theme()
#' 
#' @import ggplot2

maxwell_theme <- function(base_size = 10, 
                          grid_lines = "none", 
                          grid_color = "gray90") {
  
  # Validate grid_lines input
  grid_options <- c("none", "major", "minor", "both")
  if (!grid_lines %in% grid_options) {
    stop("grid_lines must be one of: ", paste(grid_options, collapse = ", "))
  }
  
  # Base theme
  theme_obj <- theme_minimal(base_size = base_size) +
    theme(
      text = element_text(family = "Arial", color = "black"),
      plot.title = element_text(size = rel(1.6), face = "bold", hjust = 0.5, 
                               margin = margin(b = 5)),
      axis.title.x = element_text(margin = margin(t = 10), hjust = 0.5),  
      axis.title.y = element_text(margin = margin(r = 10), hjust = 0.5), 
      axis.title = element_text(size = rel(1.4), face = "plain"),
      axis.text = element_text(size = rel(1.2), color = "black"),
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks.length = unit(0.25, "cm"),  
      axis.ticks = element_line(linewidth = 0.5, color = "black"),
      panel.border = element_blank(),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      legend.position = "top",
      legend.title = element_text(size = rel(1.2)),
      legend.text = element_text(size = rel(1.0)),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # Handle grid lines
  if (grid_lines != "none") {
    theme_obj <- theme_obj + theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = if (grid_lines %in% c("major", "both")) 
        element_line(color = grid_color, linewidth = 0.3) else element_blank(),
      panel.grid.minor.y = if (grid_lines %in% c("minor", "both")) 
        element_line(color = grid_color, linewidth = 0.15) else element_blank()
    )
  } else {
    theme_obj <- theme_obj + theme(panel.grid = element_blank())
  }
  
  return(theme_obj)
}

#' Tidy Color Palette
#' 
#' Default color palette for tidyplot graphics
#' 
#' @export
#' @return A named vector of hex colors
maxwell_colors <- function() {
  c(
    "gray" = "#929292",
    "black" = "#000000",
    "blue" = "#0C58CA",
    "orange" = "#FF8C27"
  )
}

#' Apply  Color Scale
#' 
#' Applies the default color scale to ggplot2 plots
#' 
#' @param ... Arguments passed to scale_color_manual/scale_fill_manual
#' @export
#' @return A ggplot2 scale function
scale_fill <- function(...) {
  scale_fill_manual(values = maxwell_colors(), ...)
}

#' @rdname scale_tidyplot
#' @export
scale_color <- function(...) {
  scale_color_manual(values = maxwell_colors(), ...)
}

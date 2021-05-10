library(ggrepel)
library(scales)
library(grid)


extract_ggrepel <- function(p) {
  # source: https://stackoverflow.com/questions/45065567/getting-coordinates-for-the-label-locations-from-ggrepel
  # Get x and y plot ranges 
  ggp <- ggplot_build(p)
  xrg <- ggp$layout$panel_params[[1]]$x.range
  yrg <- ggp$layout$panel_params[[1]]$y.range
  print(p)
  grid.force()
  nodes <- grid.get("textrepeltree", grep = TRUE)
  kids <- childNames(nodes)

  # fix for some upstream change BS
  #browser() 
  kids <- kids[regexpr("textrepelgrob", kids) != -1]

  # Function: get the x and y positions of a single ggrepel label
  get.xy.pos.labs <- function(n) {
    grid.force()
    grb <- grid.get(n)
    data.frame(
      x = xrg[1]+diff(xrg)*convertX(grb$x, "native", valueOnly = TRUE),
      y = yrg[1]+diff(yrg)*convertY(grb$y, "native", valueOnly = TRUE),
      label = grb$label
    )
  }
  # Get positions of all ggrepel labels
  dts <- do.call(rbind, lapply(kids, get.xy.pos.labs))
  dts
}



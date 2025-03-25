############################################
### Helper function for plotting futures ###
############################################

# 2024-07-15

# helper function for ensemble plots of potential futures

plot_futures <- function(future_tibble, output_name){
  
  p1 <- ggplot() + 
    ggtitle("SSP2-4.5") +
    theme_void(base_size = 10,
               base_family = "Helvetica Neue") +
    geom_raster(aes(x = x, y = y, fill = preds), data = future_tibble |> filter(scenario == "ssp245") |> filter(model == "ensemble")) +
    geom_sf(aes(), data = world_shp, colour = "grey60", fill = "grey60") +
    geom_sf(aes(), data = polar_buffer$mask, fill = "white", color = "grey40", size = 0.5 / .pt) +
    geom_sf(aes(), data = fronts, colour = "black", linewidth = .1) +
    geom_text(aes(x = 2700000, y = -2100000, label = "SACCF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "black") +
    geom_text(aes(x = 3050000, y = -2550000, label = "PF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "black") +
    geom_text(aes(x = 3450000, y = -2950000, label = "SAF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "black") +
    coord_sf(xlim = c(-6400000, 6400000), ylim = c(-6400000, 6400000), expand = FALSE, crs = crs_polar, ndiscr = 1000) +
    scale_fill_discrete_gradient("Change in Habitat Suitability",
                                 colours = RColorBrewer::brewer.pal(11, "RdYlBu"),
                                 bins = 12,
                                 limits = c(-.5, .5),
                                 breaks = seq(-.5, .5, 0.25),
                                 labels = seq(-.5, .5, 0.25),
                                 oob = scales::squish,
                                 guide = guide_colourbar(
                                   nbin = 500,
                                   raster = T,
                                   frame.colour = "grey40",
                                   ticks.colour = "grey40",
                                   frame.linewidth = .1,
                                   barwidth = 20,
                                   barheight = .5,
                                   direction = "horizontal",
                                   title.position = "top",
                                   title.theme = element_text(hjust = 0.5, size = 10))) +
    xlab(NULL) + ylab(NULL) +
    theme(plot.title = element_text(hjust = 0.5, size = 10, family = "Helvetica Neue"),
          legend.text = element_text(size = 10, family = "Helvetica Neue"),
          legend.title = element_text(size = 10, family = "Helvetica Neue"),
          panel.background = element_rect(fill = "black"),
          legend.position = "bottom")
  
  p2 <- ggplot() + 
    ggtitle("SSP5-8.5") +
    theme_void(base_size = 10,
               base_family = "Helvetica Neue") +
    geom_raster(aes(x = x, y = y, fill = preds), data = future_tibble |> filter(scenario == "ssp585") |> filter(model == "ensemble")) +
    geom_sf(aes(), data = world_shp, colour = "grey60", fill = "grey60") +
    geom_sf(aes(), data = polar_buffer$mask, fill = "white", color = "grey40", size = 0.5 / .pt) +
    geom_sf(aes(), data = fronts, colour = "black", linewidth = .1) +
    geom_text(aes(x = 2700000, y = -2100000, label = "SACCF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "black") +
    geom_text(aes(x = 3050000, y = -2550000, label = "PF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "black") +
    geom_text(aes(x = 3450000, y = -2950000, label = "SAF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "black") +
    coord_sf(xlim = c(-6400000, 6400000), ylim = c(-6400000, 6400000), expand = FALSE, crs = crs_polar, ndiscr = 1000) +
    scale_fill_discrete_gradient("Change in Habitat Suitability",
                                 colours = RColorBrewer::brewer.pal(11, "RdYlBu"),
                                 bins = 12,
                                 limits = c(-.5, .5),
                                 breaks = seq(-.5, .5, 0.25),
                                 labels = seq(-.5, .5, 0.25),
                                 oob = scales::squish,
                                 guide = guide_colourbar(
                                   nbin = 500,
                                   raster = T,
                                   frame.colour = "grey40",
                                   ticks.colour = "grey40",
                                   frame.linewidth = .1,
                                   barwidth = 20,
                                   barheight = .5,
                                   direction = "horizontal",
                                   title.position = "top",
                                   title.theme = element_text(hjust = 0.5, size = 10))) +
    xlab(NULL) + ylab(NULL) +
    theme(plot.title = element_text(hjust = 0.5, size = 10, family = "Helvetica Neue"),
          legend.text = element_text(size = 10, family = "Helvetica Neue"),
          legend.title = element_text(size = 10, family = "Helvetica Neue"),
          panel.background = element_rect(fill = "black"),
          legend.position = "bottom")
  
  p3 <- ggplot() + 
    theme_void(base_size = 10,
               base_family = "Helvetica Neue") +
    geom_raster(aes(x = x, y = y, fill = preds), data = future_tibble |> filter(scenario == "ssp245") |> filter(model != "ensemble")) +
    geom_sf(aes(), data = world_shp, colour = "grey60", fill = "grey60") + 
    geom_sf(aes(), data = polar_buffer$mask, fill = "white", color = "grey40", size = 0.5 / .pt) +
    coord_sf(xlim = c(-6400000, 6400000), ylim = c(-6400000, 6400000), expand = FALSE, crs = crs_polar, ndiscr = 1000) +
    scale_fill_discrete_gradient("Change in Habitat Suitability",
                                 colours = RColorBrewer::brewer.pal(11, "RdYlBu"),
                                 bins = 12,
                                 limits = c(-.5, .5),
                                 breaks = seq(-.5, .5, 0.25),
                                 labels = seq(-.5, .5, 0.25),
                                 oob = scales::squish,
                                 guide = guide_colourbar(
                                   nbin = 500,
                                   raster = T,
                                   frame.colour = "grey40",
                                   ticks.colour = "grey40",
                                   frame.linewidth = .1,
                                   barwidth = 20,
                                   barheight = .5,
                                   direction = "horizontal",
                                   title.position = "top",
                                   title.theme = element_text(hjust = 0.5, size = 10))) +
    xlab(NULL) + ylab(NULL) +
    facet_wrap(~ model) +
    theme(plot.title = element_text(hjust = 0.5, size = 10, family = "Helvetica Neue"),
          legend.text = element_text(size = 10, family = "Helvetica Neue"),
          legend.title = element_text(size = 10, family = "Helvetica Neue"),
          panel.background = element_rect(fill = "black"),
          legend.position = "bottom")
  
  p4 <- ggplot() + 
    theme_void(base_size = 10,
               base_family = "Helvetica Neue") +
    geom_raster(aes(x = x, y = y, fill = preds), data = future_tibble |> filter(scenario == "ssp585") |> filter(model != "ensemble")) +
    geom_sf(aes(), data = world_shp, colour = "grey60", fill = "grey60") + 
    geom_sf(aes(), data = polar_buffer$mask, fill = "white", color = "grey40", size = 0.5 / .pt) +
    coord_sf(xlim = c(-6400000, 6400000), ylim = c(-6400000, 6400000), expand = FALSE, crs = crs_polar, ndiscr = 1000) +
    scale_fill_discrete_gradient("Change in Habitat Suitability",
                                 colours = RColorBrewer::brewer.pal(11, "RdYlBu"),
                                 bins = 12,
                                 limits = c(-.5, .5),
                                 breaks = seq(-.5, .5, 0.25),
                                 labels = seq(-.5, .5, 0.25),
                                 oob = scales::squish,
                                 guide = guide_colourbar(
                                   nbin = 500,
                                   raster = T,
                                   frame.colour = "grey40",
                                   ticks.colour = "grey40",
                                   frame.linewidth = .1,
                                   barwidth = 20,
                                   barheight = .5,
                                   direction = "horizontal",
                                   title.position = "top",
                                   title.theme = element_text(hjust = 0.5, size = 10))) +
    xlab(NULL) + ylab(NULL) +
    facet_wrap(~ model) +
    theme(plot.title = element_text(hjust = 0.5, size = 10, family = "Helvetica Neue"),
          legend.text = element_text(size = 10, family = "Helvetica Neue"),
          legend.title = element_text(size = 10, family = "Helvetica Neue"),
          panel.background = element_rect(fill = "black"),
          legend.position = "bottom")
  
  # pretty layout
  p <- (p1 + p2) / (p3 + p4) + 
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = 'a') & theme(legend.position = "bottom",
                                              plot.tag = element_text(size = 12,
                                                                      family = "Helvetica Neue",
                                                                      face = "bold"))
  
  ggsave(filename = output_name,
         plot = p,
         path = "figures/",
         width = 7,
         height = 8.5,
         dpi = 500)
}



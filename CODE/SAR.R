############MODELOS DE AUTOREGRESION SIMULTANEA (SAR)###########################

# Load required libraries
library(sf) # Simple features for spatial data
library(spdep) # Spatial dependence tools
library(spatialreg) # Spatial regression models
library(ggplot2) # Plotting library
library(dplyr) # Data manipulation
library(gridExtra) # For arranging multiple plots
library(stats)
library(ggspatial) # For adding north arrow and scale to maps
library(htmltools)


# Set the path for saving plots
save_path <- "G:/My Drive/INVESTIGACION/PAPERS/ELABORACION/Modelo_SAR/Figures/"

# Load shapefile data (Area of Interest)
aoi = st_read("G:/My Drive/INVESTIGACION/PAPERS/ELABORACION/Modelo_SAR/DATA/df_catchments_kmeans.gpkg", quiet = TRUE)

# Scale selected columns to standardize them for comparison
aoi2 <- aoi %>% mutate(across(c('rainfallAnnual_mean', 'elev_mean', 'rel_mean', 'area'), ~(scale(.) %>% as.vector)))

# Log transform the landslide variable to normalize the data
aoi2$logy = log(aoi$lands + 1)

# Store the geometry of the polygons for later use
aoi.geom = st_geometry(aoi)

# Calculate centroids of each polygon for neighbor identification
aoi.coords = st_centroid(aoi.geom)

# Create spatial neighbors using k-nearest neighbors (k = 5)
nb_k5 = knn2nb(knearneigh(aoi.coords, k = 5))

# Convert neighbors list to spatial weights
nb_k5_list = nb2listw(nb_k5)

# Modelo Base
col.fit1 = lm(logy ~ elev_mean + rel_mean + area + rainfallAnnual_mean, data = aoi2)
summary(col.fit1)
aoi2$fit1 <- fitted(col.fit1)
moran.mc(residuals(col.fit1),nsim = 999,listw = nb_k5_list,alternative = "greater")
lmt = lm.LMtests(col.fit1, nb_k5_list, test = c("LMerr", "LMlag"))
summary(lmt)

#Spatial Lag Model
col.fit2 = lagsarlm(logy ~ elev_mean + rel_mean + area + rainfallAnnual_mean, data = aoi2,listw = nb_k5_list)
summary(col.fit2,Nagelkerke=T)
aoi2$fit2 <- fitted(col.fit2)
moran.mc(residuals(col.fit2),nsim = 999,listw = nb_k5_list,alternative = "greater")
impacts(col.fit2, listw = nb_k5_list)
summary(impacts(col.fit2, listw = nb_k5_list, R=500), zstats = TRUE) 

#Spatial Error Model
col.fit3 = errorsarlm(logy ~ elev_mean + rel_mean + area + rainfallAnnual_mean, data = aoi2,listw = nb_k5_list)
summary(col.fit3,Nagelkerke=T)
aoi2$fit3 <- fitted(col.fit3)
moran.mc(residuals(col.fit3),nsim = 999,listw = nb_k5_list,alternative = "greater")

#Spatial Durbin Lag Model
col.fit4 = lagsarlm(logy ~ elev_mean + rel_mean + area + rainfallAnnual_mean, data = aoi2,listw = nb_k5_list,type = "mixed")
summary(col.fit4,Nagelkerke=T)
aoi2$fit4 <- fitted(col.fit4)
moran.mc(residuals(col.fit4),nsim = 999,listw = nb_k5_list,alternative = "greater")
impacts(col.fit4, listw = nb_k5_list)
summary(impacts(col.fit4, listw = nb_k5_list, R=500), zstats = TRUE) 

#SLX Spatially Lagged X
col.fit5 = lmSLX(logy ~ elev_mean + rel_mean + area + rainfallAnnual_mean, data = aoi2, listw = nb_k5_list)
summary(col.fit5,Nagelkerke=T)
aoi2$fit5 <- fitted(col.fit5)
AIC(col.fit5)
moran.mc(residuals(col.fit5),nsim = 999,listw = nb_k5_list,alternative = "greater")
impacts(col.fit5, listw = nb_k5_list)
summary(impacts(col.fit5, listw = nb_k5_list, R=500), zstats = TRUE) 

#Spatial Durbin Error
col.fit6 <- errorsarlm(logy ~ elev_mean + rel_mean + area + rainfallAnnual_mean, data = aoi2, listw = nb_k5_list, etype = "emixed")
summary(col.fit6,Nagelkerke=T)
aoi2$fit6 <- fitted(col.fit6)
moran.mc(residuals(col.fit6),nsim = 999,listw = nb_k5_list,alternative = "greater")

#Mansky (all inclusive - not recommended)
col.fit7 <- sacsarlm(logy ~ elev_mean + rel_mean + area + rainfallAnnual_mean, data = aoi2,listw = nb_k5_list, type="sacmixed") 
summary(col.fit7,Nagelkerke=T)
aoi2$fit7 <- fitted(col.fit7)
moran.mc(residuals(col.fit7),nsim = 999,listw = nb_k5_list,alternative = "greater")
impacts(col.fit7, listw = nb_k5_list)
summary(impacts(col.fit7, listw = nb_k5_list, R=500), zstats = TRUE) 

#SARAR, Kelejian-Prucha, Cliff-Ord, or SAC
col.fit8 <- sacsarlm(logy ~ elev_mean + rel_mean + area + rainfallAnnual_mean, data = aoi2,listw = nb_k5_list, type="sac")
summary(col.fit8,Nagelkerke=T)
aoi2$fit8 <- fitted(col.fit8)
moran.mc(residuals(col.fit8),nsim = 999,listw = nb_k5_list,alternative = "greater")
summary(impacts(col.fit8, listw = nb_k5_list, R=500), zstats = TRUE) 

#Comparing
aic.tbl = AIC(col.fit1 ,col.fit2, col.fit3, col.fit4, col.fit5, col.fit6, col.fit7, col.fit8)
rownames(aic.tbl) = c("Modelo Base","Lag Model", "Error Model", "Durbin Model","SLX Spatially Lagged X","Spatial Durbin Error","Mansky","SAC")
aic.tbl %>%kbl(caption = "AIC Comparison") %>%kable_classic_2(full_width = F)

########################################################################
# histograma de frecuencia

# Plot histograms to visualize the distribution of the new variables
max_y <- max(table(cut(aoi2$lands_rec, breaks = 30)), table(cut(aoi2$logy, breaks = 30)))

# Plot histograms to visualize the distribution of the new variables
p1 <- ggplot(aoi2, aes(x = lands_rec)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  xlab("Frequencia") +
  ylim(0, 250) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.y = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p2 <- ggplot(aoi2, aes(x = logy)) +
  geom_histogram(fill = "lightgreen", color = "black", bins = 30) +
  xlab("Log(frecuencia+1)") +
  ylim(0, 250) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Arrange all three plots in a 1-row, 3-column layout and save the output
combined_plot <- grid.arrange(p1, p2, nrow = 1, ncol = 2)
ggsave(filename = paste0(save_path, "Histograms.png"), plot = combined_plot, width = 15, height = 5, units = "in", dpi = 500)

########################################################################
# Matriz de vecindad

# Plot the geometry of the areas and their neighbor relationships and save the plot
neighbor_plot <- ggplot() +
  geom_sf(data = aoi.geom, fill = NA, color = "black") +
  geom_segment(data = do.call(rbind, lapply(1:length(nb_k5), function(i) {
    data.frame(
      x = st_coordinates(aoi.coords)[i, 1],
      y = st_coordinates(aoi.coords)[i, 2],
      xend = st_coordinates(aoi.coords)[nb_k5[[i]], 1],
      yend = st_coordinates(aoi.coords)[nb_k5[[i]], 2]
    )
  })), aes(x = x, y = y, xend = xend, yend = yend), color = "red", alpha = 0.5) +
  geom_sf(data = st_as_sf(aoi.coords), color = "blue", size = 1.5) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  annotation_north_arrow(location = "tl", which_north = "true", height = unit(1, "cm"), width = unit(1, "cm")) +
  annotation_scale(location = "br", style = "ticks")

ggsave(filename = paste0(save_path, "Matriz_vecindad.png"), plot = neighbor_plot, width = 10, height = 7, units = "in", dpi = 500)

########################################################################
#plot mapa y residuo

# Plot residuals and approximated fitted values in a 1-row, 2-column layout
g1 <- ggplot() +
  geom_sf(data = aoi, aes(fill = res19), color = "black") +
  annotation_scale(location = "br", style = "ticks") +
  annotation_north_arrow(location = "tr", which_north = "true", height = unit(0.7, "cm"), width = unit(0.6, "cm")) +
  scale_fill_viridis_c(name = "Residuals") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks.length = unit(-0.1, "cm"),
    axis.text.x = element_text(size = 8, margin = unit(c(t = 1, r = 0, b = 0, l = 0), "mm")),
    axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 1, b = 0, l = 0), "mm")),
    legend.text = element_text(size = 6),
    legend.title.align = 0,
    legend.position = c(0.3, 0.9),
    legend.key.size = unit(0.5, 'cm'),
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.title = element_text(size = 10, vjust = .8, hjust = .5)
  )

g2 <- ggplot() +
  geom_sf(data = aoi, aes(fill = fit19_approx), color = "black") +
  annotation_scale(location = "br", style = "ticks") +
  annotation_north_arrow(location = "tr", which_north = "true", height = unit(0.7, "cm"), width = unit(0.6, "cm")) +
  scale_fill_gradientn(colors = c("blue", "white", "red"), name = "Fitted Values (Approx)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks.length = unit(-0.1, "cm"),
    axis.text.x = element_text(size = 8, margin = unit(c(t = 1, r = 0, b = 0, l = 0), "mm")),
    axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 1, b = 0, l = 0), "mm")),
    legend.text = element_text(size = 6),
    legend.title.align = 0,
    legend.position = c(0.3, 0.9),
    legend.key.size = unit(0.5, 'cm'),
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.title = element_text(size = 10, vjust = .8, hjust = .5)
  )

combined_plot_2 <- grid.arrange(g1, g2, nrow = 1, ncol = 2)
ggsave(filename = paste0(save_path, "Residuals_and_Fitted_Values.png"), plot = combined_plot_2, width = 15, height = 7, units = "in", dpi = 500)

#######################################################################

# Load required libraries
library(spatialreg)
library(spdep)
library(sf)
library(dplyr)
library(xtable)

# Load shapefile data (Area of Interest)
aoi <- st_read("G:/My Drive/INVESTIGACION/PAPERS/ELABORACION/Modelo_SAR/DATA/df_catchments_kmeans.gpkg", quiet = TRUE)

# Ensure the continuous variables are numeric before scaling
aoi <- aoi %>%
  mutate(across(c('area', 'hypso_inte', 'slope_mean', 'rainfallAnnual_mean'), as.numeric))

# Standardize the continuous variables
aoi <- aoi %>%
  mutate(across(c('area', 'hypso_inte', 'slope_mean', 'rainfallAnnual_mean'), scale))

# Log transform the landslide variable to normalize the data
aoi$y_log <- log(aoi$lands_rec + 1)

# Store the geometry of the polygons for later use
aoi.geom <- st_geometry(aoi)

# Function to extract and format model summary
extract_model_summary <- function(model, listw) {
  summary_model <- summary(model, Nagelkerke = TRUE)
  
  # Extract coefficients, standard errors, and p-values
  coefficients <- summary_model$Coef[, "Estimate"]
  std_errors <- summary_model$Coef[, "Std. Error"]
  p_values <- summary_model$Coef[, "Pr(>|z|)"]
  lambda <- if (!is.null(summary_model$lambda)) round(summary_model$lambda, 3) else NA
  nagelkerke <- if (!is.null(summary_model$NK)) round(summary_model$NK, 3) else NA
  aic <- round(AIC(model), 3)
  
  # Moran's I test for residuals
  moran_res <- moran.mc(residuals(model), nsim = 999, listw = listw, alternative = "greater")
  moran_I <- round(moran_res$statistic, 3)
  p_value_moran <- round(moran_res$p.value, 3)
  
  return(list(coefficients = coefficients, std_errors = std_errors, p_values = p_values, 
              lambda = lambda, nagelkerke = nagelkerke, aic = aic, moran_I = moran_I, p_value_moran = p_value_moran))
}

# Format coefficients with standard errors and bold significant ones
format_coef <- function(coef, std_err, p_value) {
  if (!is.na(p_value) && p_value < 0.05) {
    return(paste0("\\textbf{", sprintf("%.2f", coef), "} (", sprintf("%.2f", std_err), ")"))
  } else {
    return(paste0(sprintf("%.2f", coef), " (", sprintf("%.2f", std_err), ")"))
  }
}

# Initialize an empty data frame to store all results
summary_results <- data.frame()

# Loop over distances from 5000 to 105000 with steps of 20000
distances <- seq(20000, 100000, by = 20000)

for (d in distances) {
  # Calculate centroids of each polygon for neighbor identification
  aoi.coords <- st_centroid(aoi.geom)
  
  # Create spatial neighbors using dnearneigh for a specified distance
  dnb <- dnearneigh(aoi.coords, d1 = 0, d2 = d)
  
  # Convert the neighbors list to a weights list using inverse distances as weights
  listw_obj <- nb2listw(dnb, glist = lapply(nbdists(dnb, aoi.coords), function(dist) 1/dist), style = "W")
  
  # Fit the Spatial Durbin Error Model (SDEM)
  SDEM_model <- errorsarlm(y_log ~ area + hypso_inte + slope_mean + rainfallAnnual_mean, data = aoi, listw = listw_obj, etype = "emixed")
  SDEM_summary <- extract_model_summary(SDEM_model, listw_obj)
  
  # Format coefficients for the table
  formatted_coefs <- sapply(1:length(SDEM_summary$coefficients), function(i) {
    format_coef(SDEM_summary$coefficients[i], SDEM_summary$std_errors[i], SDEM_summary$p_values[i])
  })
  
  # Create a data frame for this model's results
  model_results <- data.frame(
    Parameter = c("Intercepto", "Área", "Hipso", "Pendiente", "Lluvia", 
                  "Wx-Área", "Wx-Hipso", "Wx-Pendiente", "Wx-Lluvia", "Lambda", "Nagelkerke_R2", "AIC", "Moran_I", "p-value Moran"),
    Value = c(formatted_coefs, SDEM_summary$lambda, SDEM_summary$nagelkerke, SDEM_summary$aic, SDEM_summary$moran_I, SDEM_summary$p_value_moran),
    Distance = as.character(d),  # Use distance as a string for reshaping purposes
    stringsAsFactors = FALSE
  )
  
  # Append the results to the summary data frame
  summary_results <- rbind(summary_results, model_results)
}

# Reshape the data frame for LaTeX table
summary_table <- reshape(summary_results, idvar = "Parameter", timevar = "Distance", direction = "wide")

# Adjust column names to only show the distance values
colnames(summary_table) <- gsub("Value\\.", "", colnames(summary_table))

# Create a LaTeX table using xtable
latex_table <- xtable(summary_table, caption = "Model Results for SDEM at Different Distances", label = "tab:sdem_summary")

# Print the LaTeX table without escaping LaTeX commands and with custom formatting
print(latex_table, type = "latex", include.rownames = FALSE, sanitize.text.function = identity)

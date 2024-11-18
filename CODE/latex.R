# Load required libraries
library(spatialreg)
library(spdep)
library(sf)
library(dplyr)
library(xtable)

# Load shapefile data (Area of Interest)
aoi = st_read("G:/My Drive/INVESTIGACION/PAPERS/ELABORACION/Modelo_SAR/DATA/df_catchments_kmeans.gpkg", quiet = TRUE)

# Ensure the continuous variables are numeric before scaling
aoi <- aoi %>%
  mutate(across(c('area', 'hypso_inte', 'slope_mean', 'rainfallAnnual_mean'), as.numeric))

# Standardize the continuous variables
aoi <- aoi %>%
  mutate(across(c('area', 'hypso_inte', 'slope_mean', 'rainfallAnnual_mean'), scale))

# Log transform the landslide variable to normalize the data
aoi$y_log = log(aoi$lands_rec + 1)

# Store the geometry of the polygons for later use
aoi.geom = st_geometry(aoi)

# Calculate centroids of each polygon for neighbor identification
aoi.coords = st_centroid(aoi.geom)

# Create spatial 
dnb <- dnearneigh(aoi.coords, d1 = 0, d2 = 20000) 

# Convert the neighbors list to a weights list using inverse distances as weights
listw_obj <- nb2listw(dnb, glist = lapply(nbdists(dnb, aoi.coords), function(dist) 1/dist), style = "W")


# Initialize an empty data frame to store the results
results <- data.frame(
  Parameter = c("Intercepto", "Área", "Hipso", "Pendiente", "Lluvia", 
                "Wx-Área", "Wx-Hipso", "Wx-Pendiente", "Wx-Lluvia", "Rho", "Lambda", "Nagelkerke_R2", "AIC"),
  SAR = NA, SEM = NA, SAC = NA, SLX = NA, SDEM = NA, SDM = NA, GNS = NA,
  stringsAsFactors = FALSE
)

# Function to extract model summary information with standard errors and p-values
extract_model_summary <- function(model) {
  summary_model <- summary(model, Nagelkerke = TRUE)
  
  coefficients <- coef(summary_model)
  std_errors <- summary_model$Coef[, "Std. Error"]  # Extract standard errors
  p_values <- summary_model$Coef[, "Pr(>|z|)"]  # Extract p-values
  rho <- if (!is.null(summary_model$rho)) summary_model$rho else NA
  lambda <- if (!is.null(summary_model$lambda)) summary_model$lambda else NA
  nagelkerke <- if (!is.null(summary_model$NK)) summary_model$NK else NA
  aic <- AIC(model)
  
  return(list(coefficients = coefficients, std_errors = std_errors, p_values = p_values, 
              rho = rho, lambda = lambda, nagelkerke = nagelkerke, aic = aic))
}

# Format coefficients with standard errors and bold significant ones using sprintf for precise formatting
format_coef <- function(coef, std_err, p_value) {
  if (!is.na(p_value) && p_value < 0.05) {
    return(paste0("\\textbf{", sprintf("%.3f", coef), "} (", sprintf("%.2f", std_err), ")"))
  } else {
    return(paste0(sprintf("%.3f", coef), " (", sprintf("%.2f", std_err), ")"))
  }
}

# Spatial Lag Model (SAR)
SAR_model <- lagsarlm(y_log ~ area + hypso_inte + slope_mean + rainfallAnnual_mean, data = aoi, listw = listw_obj)
SAR_summary <- extract_model_summary(SAR_model)

results$SAR <- c(
  format_coef(SAR_summary$coefficients[1], SAR_summary$std_errors[1], SAR_summary$p_values[1]),  # Constant
  format_coef(SAR_summary$coefficients[3], SAR_summary$std_errors[3], SAR_summary$p_values[3]),  
  format_coef(SAR_summary$coefficients[2], SAR_summary$std_errors[2], SAR_summary$p_values[2]),  
  format_coef(SAR_summary$coefficients[4], SAR_summary$std_errors[4], SAR_summary$p_values[4]),  
  format_coef(SAR_summary$coefficients[5], SAR_summary$std_errors[5], SAR_summary$p_values[5]),  
  NA, NA, NA, NA,  # Wx_A, Wx_E, Wx_H, Wx_R (Not applicable for SAR)
  round(SAR_summary$rho, 3),  # Rho
  NA,  # Lambda (Not applicable for SAR)
  round(SAR_summary$nagelkerke, 3),  # Nagelkerke R^2
  round(SAR_summary$aic, 3)  # AIC
)

# Spatial Error Model (SEM)
SEM_model <- errorsarlm(y_log ~ area + hypso_inte + slope_mean + rainfallAnnual_mean, data = aoi, listw = listw_obj)
SEM_summary <- extract_model_summary(SEM_model)

results$SEM <- c(
  format_coef(SEM_summary$coefficients[1], SEM_summary$std_errors[1], SEM_summary$p_values[1]),  # Constant
  format_coef(SEM_summary$coefficients[3], SEM_summary$std_errors[3], SEM_summary$p_values[3]),  
  format_coef(SEM_summary$coefficients[2], SEM_summary$std_errors[2], SEM_summary$p_values[2]), 
  format_coef(SEM_summary$coefficients[4], SEM_summary$std_errors[4], SEM_summary$p_values[4]),  
  format_coef(SEM_summary$coefficients[5], SEM_summary$std_errors[5], SEM_summary$p_values[5]),  
  NA, NA, NA, NA,  # Wx_A, Wx_E, Wx_H, Wx_R (Not applicable for SEM)
  NA,  # Rho (Not applicable for SEM)
  round(SEM_summary$lambda, 3),  # Lambda
  round(SEM_summary$nagelkerke, 3),  # Nagelkerke R^2
  round(SEM_summary$aic, 3)  # AIC
)

# SAC (SARAR)
SAC_model <- sacsarlm(y_log ~ area + hypso_inte + slope_mean + rainfallAnnual_mean, data = aoi, listw = listw_obj, type = "sac")
SAC_summary <- extract_model_summary(SAC_model)

results$SAC <- c(
  format_coef(SAC_summary$coefficients[1], SAC_summary$std_errors[1], SAC_summary$p_values[1]),  
  format_coef(SAC_summary$coefficients[3], SAC_summary$std_errors[3], SAC_summary$p_values[3]),  
  format_coef(SAC_summary$coefficients[2], SAC_summary$std_errors[2], SAC_summary$p_values[2]), 
  format_coef(SAC_summary$coefficients[4], SAC_summary$std_errors[4], SAC_summary$p_values[4]),  
  format_coef(SAC_summary$coefficients[5], SAC_summary$std_errors[5], SAC_summary$p_values[5]),  
  NA, NA, NA, NA,  # Wx_A, Wx_E, Wx_H, Wx_R (Not applicable for SAC)
  round(SAC_summary$rho, 3),  # Rho
  round(SAC_summary$lambda, 3),  # Lambda
  round(SAC_summary$nagelkerke, 3),  # Nagelkerke R^2
  round(SAC_summary$aic, 3)  # AIC
)

# SLX (Spatially Lagged X Model)
SLX_model <- lmSLX(y_log ~ area + hypso_inte + slope_mean + rainfallAnnual_mean, data = aoi, listw = listw_obj)
SLX_summary <- summary(SLX_model)

results$SLX <- c(
  format_coef(coef(SLX_summary)[1], SLX_summary$coefficients[1, 2], SLX_summary$coefficients[1, 4]),  
  format_coef(coef(SLX_summary)[3], SLX_summary$coefficients[3, 2], SLX_summary$coefficients[3, 4]),  
  format_coef(coef(SLX_summary)[2], SLX_summary$coefficients[2, 2], SLX_summary$coefficients[2, 4]),  
  format_coef(coef(SLX_summary)[4], SLX_summary$coefficients[4, 2], SLX_summary$coefficients[4, 4]),  
  format_coef(coef(SLX_summary)[5], SLX_summary$coefficients[5, 2], SLX_summary$coefficients[5, 4]),  
  format_coef(coef(SLX_summary)[7], SLX_summary$coefficients[7, 2], SLX_summary$coefficients[7, 4]),  
  format_coef(coef(SLX_summary)[6], SLX_summary$coefficients[6, 2], SLX_summary$coefficients[6, 4]),  
  format_coef(coef(SLX_summary)[8], SLX_summary$coefficients[8, 2], SLX_summary$coefficients[8, 4]),  
  format_coef(coef(SLX_summary)[9], SLX_summary$coefficients[9, 2], SLX_summary$coefficients[9, 4]),  
  NA,  # Rho (Not applicable for SLX)
  NA,  # Lambda (Not applicable for SLX)
  NA,  # Nagelkerke (Not provided for SLX)
  round(AIC(SLX_model), 3)  # AIC
)

# Spatial Durbin Error Model (SDEM)
SDEM_model <- errorsarlm(y_log ~ area + hypso_inte + slope_mean + rainfallAnnual_mean, data = aoi, listw = listw_obj, etype = "emixed")
SDEM_summary <- extract_model_summary(SDEM_model)

results$SDEM <- c(
  format_coef(SDEM_summary$coefficients[1], SDEM_summary$std_errors[1], SDEM_summary$p_values[1]),  
  format_coef(SDEM_summary$coefficients[3], SDEM_summary$std_errors[3], SDEM_summary$p_values[3]),  
  format_coef(SDEM_summary$coefficients[2], SDEM_summary$std_errors[2], SDEM_summary$p_values[2]),  
  format_coef(SDEM_summary$coefficients[4], SDEM_summary$std_errors[4], SDEM_summary$p_values[4]),  
  format_coef(SDEM_summary$coefficients[5], SDEM_summary$std_errors[5], SDEM_summary$p_values[5]),  
  format_coef(SDEM_summary$coefficients[7], SDEM_summary$std_errors[7], SDEM_summary$p_values[7]),  
  format_coef(SDEM_summary$coefficients[6], SDEM_summary$std_errors[6], SDEM_summary$p_values[6]),  
  format_coef(SDEM_summary$coefficients[8], SDEM_summary$std_errors[8], SDEM_summary$p_values[8]),  
  format_coef(SDEM_summary$coefficients[9], SDEM_summary$std_errors[9], SDEM_summary$p_values[9]),  
  NA,  # Rho (Not applicable for SDEM)
  round(SDEM_summary$lambda, 3),  # Lambda
  round(SDEM_summary$nagelkerke, 3),  # Nagelkerke R^2
  round(SDEM_summary$aic, 3)  # AIC
)

# Spatial Durbin Model (SDM)
SDM_model <- lagsarlm(y_log ~ area + hypso_inte + slope_mean + rainfallAnnual_mean, data = aoi, listw = listw_obj, type = "mixed")
SDM_summary <- extract_model_summary(SDM_model)

results$SDM <- c(
  format_coef(SDM_summary$coefficients[1], SDM_summary$std_errors[1], SDM_summary$p_values[1]),  
  format_coef(SDM_summary$coefficients[3], SDM_summary$std_errors[3], SDM_summary$p_values[3]),  
  format_coef(SDM_summary$coefficients[2], SDM_summary$std_errors[2], SDM_summary$p_values[2]),  
  format_coef(SDM_summary$coefficients[4], SDM_summary$std_errors[4], SDM_summary$p_values[4]),  
  format_coef(SDM_summary$coefficients[5], SDM_summary$std_errors[5], SDM_summary$p_values[5]),  
  format_coef(SDM_summary$coefficients[7], SDM_summary$std_errors[7], SDM_summary$p_values[7]),  
  format_coef(SDM_summary$coefficients[6], SDM_summary$std_errors[6], SDM_summary$p_values[6]),  
  format_coef(SDM_summary$coefficients[8], SDM_summary$std_errors[8], SDM_summary$p_values[8]),  
  format_coef(SDM_summary$coefficients[9], SDM_summary$std_errors[9], SDM_summary$p_values[9]),  
  round(SDM_summary$rho, 3),  # Rho
  NA,  # Lambda (Not applicable for SDM)
  round(SDM_summary$nagelkerke, 3),  # Nagelkerke R^2
  round(SDM_summary$aic, 3)  # AIC
)

# Mansky Generalized Nesting Spatial (GNS) Model
GNS_model <- sacsarlm(y_log ~ area + hypso_inte + slope_mean + rainfallAnnual_mean, data = aoi, listw = listw_obj, type = "sacmixed")
GNS_summary <- extract_model_summary(GNS_model)

results$GNS <- c(
  format_coef(GNS_summary$coefficients[1], GNS_summary$std_errors[1], GNS_summary$p_values[1]), 
  format_coef(GNS_summary$coefficients[3], GNS_summary$std_errors[3], GNS_summary$p_values[3]),  
  format_coef(GNS_summary$coefficients[2], GNS_summary$std_errors[2], GNS_summary$p_values[2]),  
  format_coef(GNS_summary$coefficients[4], GNS_summary$std_errors[4], GNS_summary$p_values[4]),  
  format_coef(GNS_summary$coefficients[5], GNS_summary$std_errors[5], GNS_summary$p_values[5]),  
  format_coef(GNS_summary$coefficients[7], GNS_summary$std_errors[7], GNS_summary$p_values[7]),  
  format_coef(GNS_summary$coefficients[6], GNS_summary$std_errors[6], GNS_summary$p_values[6]),  
  format_coef(GNS_summary$coefficients[8], GNS_summary$std_errors[8], GNS_summary$p_values[8]),  
  format_coef(GNS_summary$coefficients[9], GNS_summary$std_errors[9], GNS_summary$p_values[9]),  
  round(GNS_summary$rho, 3),  # Rho
  round(GNS_summary$lambda, 3),  # Lambda
  round(GNS_summary$nagelkerke, 3),  # Nagelkerke R^2
  round(GNS_summary$aic, 3)  # AIC
)

# Create a LaTeX table
latex_table <- xtable(results, caption = "Model Results for Various Spatial Models", label = "tab:spatial_models")

# Print LaTeX table without escaping LaTeX commands and with custom hline
print(latex_table, type = "latex", include.rownames = FALSE, hline.after = c(0, 7, 9), sanitize.text.function = identity)


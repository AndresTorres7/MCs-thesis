

setwd("C:/Users/thbio/Documents/Biología/Maestría/Tesis/Avances/modogeografico")

occurrence<-read.csv("S33_distribution_maps/Micrathena_all_distributions.csv")
occurrence$Species <- gsub(" ", "_", occurrence$Species)
occurrence <- occurrence[,c(1,3,2)]

## changing col 1 name
names(occurrence) [1] <- paste("species")
names(occurrence)[2] <- paste("longitude")
names(occurrence)[3] <- paste("latitude")

setwd("E:/ModelodistribucionMicrathena")

# Cargar librerías necesarias
library(dplyr)


dir.create("Occurrences")
setwd("./Occurrences")

# Crear carpetas y guardar archivos CSV
for (especie in unique(occurrence$species)) {
  # Filtrar los datos por especie
  datos_especie <- occurrence %>% filter(occurrence$species == especie)
  

  # Crear carpeta si no existe
  if (!dir.exists(especie)) {
    dir.create(especie)
  }
  
  
  test<- datos_especie[j<-sample(nrow(datos_especie), round((length(datos_especie[,1])/2))), ]
  train <- datos_especie[-j,]
  
  
  # Guardar el archivo CSV en la carpeta correspondiente
  write.csv(datos_especie, file.path(especie, paste0(especie, "_join.csv")), row.names = F)
  write.csv(train, file.path(especie, paste0(especie, "_train.csv")), row.names = F)
  write.csv(test, file.path(especie, paste0(especie, "_test.csv")), row.names = F)
  }



#occt <- thin_data(occ, "longitude", "latitude", thin_distance = 1, save = T, #Change the Km based on environmental resolution 
                  #name = "Occ_Thinning_1km")

#tablaGenus$numOcc[k] <- nrow(occt)  ## n???mero de ourrencias por especie luego de limpieza





















##### clipping raster with america shapefile

install.packages("sf")
install.packages("RStoolbox")

# Cargar las librerías necesarias
library(sf)
library(terra)
library(RStoolbox)

# Directorio donde se encuentran los archivos raster
raster_dir <- "E:/ModelodistribucionMicrathena/env/worldclim"
newraster_dir <- "E:/ModelodistribucionMicrathena/env/worldclimclipped"

# Cargar los archivos raster
raster_files <- list.files(raster_dir, pattern = "\\.tif$", full.names = TRUE)
rasters <- lapply(raster_files, rast)

plot(rasters[[1]])

# Cargar el shapefile
shape <- st_read("E:/ModelodistribucionMicrathena/env/shape/shapeamerica.shp")
ext(shape)

plot(shape)

# Verificar y transformar la proyección si es necesario
shape_proj <- st_crs(shape)

for (i in 1:length(rasters)) {
  if (st_crs(rasters[[i]])$input != shape_proj$input) {
   cat(paste("el raster " , i ," no tiene la misma proyeccion"))
  }
}

# Crear una lista para los rasters recortados
clipped_rasters <- list()



# Recortar cada raster con el shapefile
for (i in 1:length(rasters)) {
  cropped_raster <- crop(rasters[[i]], ext(shape))
  clipped_raster <- mask(cropped_raster, vect(shape))
  
  clipped_rasters[[i]] <- clipped_raster
  
  # Guardar el raster recortado
  output_filename <- paste0("clipped_", basename(raster_files[i]))
  
  writeRaster(clipped_rasters[[i]], filename = file.path(newraster_dir, output_filename), overwrite = TRUE)
}




######  Make PCAs






########### Make M

# Cargar las librerías necesarias
library(terra)
library(sf)

# Cargar las capas raster

raster1 <- rast("E:/ModelodistribucionMicrathena/env/worldclimclipped/PCAtif/PC1.tif")
raster2 <- rast("E:/ModelodistribucionMicrathena/env/worldclimclipped/PCAtif/PC2.tif")
raster3 <- rast("E:/ModelodistribucionMicrathena/env/worldclimclipped/PCAtif/PC3.tif")
raster4 <- rast("E:/ModelodistribucionMicrathena/env/worldclimclipped/PCAtif/PC4.tif")

# Leer los datos de ocurrencia desde un archivo CSV
ocurrence_data <- read.csv("E:/ModelodistribucionMicrathena/Occurrences/Micrathena_gracilis/Micrathena_gracilis_join.csv")

# Convertir los datos de ocurrencia en un objeto spatial
ocurrence_sf <- st_as_sf(ocurrence_data, coords = c("Longitude", "Latitude"), crs = 4326)

# Crear un buffer de 30 km alrededor de cada punto de ocurrencia
ocurrence_buffer <- st_buffer(ocurrence_sf, dist = 300000)  # acomodar la distancia, está en metros

# Unir los buffers en un único polígono (si es necesario)
buffer_union <- st_union(ocurrence_buffer)

# Convertir a un objeto spatVector para usar con terra
buffer_spatvector <- vect(buffer_union)



# Recortar y enmascarar los raster layers usando el polígono del buffer
raster1_masked <- mask(crop(raster1, buffer_spatvector), buffer_spatvector)
raster2_masked <- mask(crop(raster2, buffer_spatvector), buffer_spatvector)
raster3_masked <- mask(crop(raster3, buffer_spatvector), buffer_spatvector)
raster4_masked <- mask(crop(raster4, buffer_spatvector), buffer_spatvector)

setwd("E:/ModelodistribucionMicrathena/micrathenagracilistest/Mbuffer")

# Guardar los raster layers enmascarados


writeRaster(raster1_masked, filename = "M1.asc", datatype = "FLT4S", NAflag = 999, overwrite = TRUE)
writeRaster(raster2_masked, filename = "M2.asc", datatype = "FLT4S", NAflag = 999, overwrite = TRUE)
writeRaster(raster3_masked, filename = "M3.asc", datatype = "FLT4S", NAflag = 999, overwrite = TRUE)
writeRaster(raster4_masked, filename = "M4.asc", datatype = "FLT4S", NAflag = 999, overwrite = TRUE)







##### modeling with kuenm


devtools::install_github("marlonecobos/kuenm")
library(kuenm)


kuenm_start(file.name = "gracilis_test")
setwd("E:/ModelodistribucionMicrathena/micrathenagracilistest/Micrathena_gracilis")

kuenm_varcomb("Mbuffer", "Mbu")


occ_joint <- "Micrathena_gracilis_join.csv"
occ_tra <- "Micrathena_gracilis_train.csv"
M_var_dir <- "Mbu"
batch_cal <- "Candidate_models"
out_dir <- "Candidate_Models"
reg_mult <- c(seq(0.1, 1, 0.1), 2:5)
f_clas <- "all"
args <- NULL
maxent_path <- "C:/Users/thbio/Desktop/maxent"
wait <- FALSE
run <- TRUE


kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, 
          batch = batch_cal, out.dir = out_dir, reg.mult = reg_mult, 
          f.clas = f_clas, args = args, maxent.path = maxent_path, 
          wait = wait, run = run)

occ_test <- "Micrathena_gracilis_test.csv"
out_eval <- "Calibration_results"
threshold <- 5
rand_percent <- 30
iterations <- 100
kept <- F
selection <- "OR_AICc"


cal_eval <- kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, 
                        occ.test = occ_test, batch = batch_cal, out.eval = out_eval, 
                        threshold = threshold, rand.percent = rand_percent, 
                        iterations = iterations, kept = kept, selection = selection)

?kuenm_ceval



library(sf)
library(concaveman)
library(raster)
library(readxl)
library(purrr)



# Final location of rasters
dir.create("polirdaS")


setwd("D:/OccurrenceData")
generos <- dir()
genusList <- data.frame(Genus = generos)

#Adiantum 3
#Anolis 8
# Chusquea
#nrow(genusList)
for(v in 120:nrow(genusList)) {
  k <- list()
  
  
  setwd(paste("D:/OccurrenceData", genusList$Genus[v], sep = "/"))  
  spplist <- dir()
  spplist <-  spplist[!grepl("_rasters\\.rda$",  spplist)]
  cat(genusList$Genus[v])
  
  # create general extension

  df <- data.frame()
  
  for(z in 1:length(spplist)) {
    
    setwd(file.path("D:/OccurrenceData", genusList$Genus[v], spplist[z]))
    poins <- read.csv(list.files(pattern = "_join.csv"))
    df <- rbind(df, poins)
    
  }
  
  
  pts <- st_as_sf(df, coords = c("longitude", "latitude"), crs = 4326)
  poly <- concaveman(pts, concavity = 3)
  pla <- raster::extent(poly)
  gc()
  
  #length(spplist) 154, 135 , 127, 150, 60, 95 v = 8
  for (i in 1:length(spplist)) {
    
    #setwd(file.path("D:/OccurrenceData", genusList$Genus[v], spplist[i]))
    
    #setwd(file.path("Final_Models"))  
    #models <- dir()
    
    
   # if(!file.exists("polygonModel.tif")) {
    
    setwd(file.path("D:/OccurrenceData", genusList$Genus[v], spplist[i]))

    # Leer el archivo CSV (ajustar ruta si es necesario)
    df <- read.csv(list.files(pattern = "_join.csv"))
    
    # Verifica que las columnas se llaman lon y lat (ajusta si es necesario)
    head(df)
    # Crear puntos sf
    pts <- st_as_sf(df, coords = c("longitude", "latitude"), crs = 4326)
  #  plot(pts, key.width = lcm(6))
    # Crear polígono ajustado (concave hull)
    poly <- concaveman(pts, concavity = 3)  # Puedes probar con valores de 1 a 5
   # plot(poly, key.width = lcm(6))
    
    ### crear un mega poligono de todos los puntos para calcular la extensión máxima. 
    
    # Crear raster de fondo con resolución adecuada
    r <- raster(raster::extent(poly), res = 0.05, crs = "+proj=longlat +datum=WGS84")
    
    # Rasterizar el polígono (valor 1 donde hay presencia)
    r_mask <- rasterize(as_Spatial(poly), r, field = 1)
    
   #plot(r_mask, key.width = lcm(6))
    
    r_mask <- raster::extend(r_mask, pla, NA)
    
    # Si todo es NA (polígono perdido en la resolución), forzar al menos un píxel con valor 1
    if (all(is.na(values(r_mask)))) {
      centroid <- st_centroid(st_geometry(poly))
      cell <- raster::cellFromXY(r_mask, st_coordinates(centroid))
      r_mask[cell] <- 1
    }
    
    
    # Visualizar
    plot(r_mask, key.width = lcm(9))
 
        
    
    setwd(file.path("Final_Models"))  
    
    
    raster::writeRaster(r_mask, "polygonModel.tif", overwrite = TRUE)
 # }
  
    
    cat(spplist[i])
    gc()
    dev.off()
    
  }
  
  for (l in  1:length(spplist)) {
    cat("Gathering...") 
    setwd(file.path("D:/OccurrenceData", genusList$Genus[v], spplist[l],"Final_Models"))
    modelbinar <- raster("polygonModel.tif")
    modelbinar <- readAll(modelbinar)
    gc()
    k <- c(k, list(modelbinar)) 
    gc()
    
  }
  
  
  fornames <- spplist
  
  setwd("D:/G/polirdaS")
  
  cat("Writing \\")
  names(k) <- fornames
  save(k, file = paste(genusList$Genus[v],"polirasters.rda", sep = "_"))
  
  
}





#"Anolis_kreutzi" 2       "Anolis_amplisquamosus"2

#"Anolis_purpurgularis" 2
#"Anolis_zapotecorum" no
# "Anolis_stevepoei"
#"Anolis_rubribarbaris" 
## "Anolis_wattsii"
##  "Anolis_ginaelisae"
 #"Anolis_morazani"




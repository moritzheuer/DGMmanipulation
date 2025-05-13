library("raster")
library("sp")
library("remotes")
library("whitebox")
library(plyr)
library(dplyr)
library(tmap)
library(gridExtra)
library(terra)
library(ggplot2)
library(viridis)
library(readr)
library("rdwplus")
library(tidyverse)

dem <- raster("S:/DGM_Flow_Modellierung/07_data_preprocessing/Argenschwang_DEM_50cm.tif", crs = '+init=EPSG:32632')
Auflösung <- 0.5

options(digits=20)

df <- read_csv("S:/DGM_Flow_Modellierung/07_data_preprocessing/QGIS_Projekte/Graben_Vertices_Koordinaten.csv")

df <- df[,c("id", "Ordnung", "vertex_index", "xcoord", "ycoord")]

if (Auflösung == 0.5) {
  x_min <- 401384.25
  x_max <- 402434.25
  y_min <- 5530950.75
  y_max <- 5532365.75
  
  df$xcoord <- df$xcoord - x_min
  df$xcoord <- floor(df$xcoord / 0.5) + 1
  
  df$ycoord <- df$ycoord - y_min
  df$ycoord <- floor(df$ycoord / 0.5) + 1
  df$ycoord <- 2830 - df$ycoord
  df$ycoord <- df$ycoord + 1
}

colnames(df) <- c("id", "Ordnung", "vertex_index", "column", "row")

###########################################
#Nachprägen der Gräben

zähler <- 1
grabenList <- c()
for (ordnung in sort(unique(df$Ordnung))) {
  for (id in sort(unique(df$id))) {
      #grabenZellen <- max(df[df$id == id & df$Ordnung == ordnung,]$vertex_index)
    if (nrow(df[df$id == id & df$Ordnung == ordnung,]) != 0) {
      print(c(ordnung, id))
      partDrainage <- bresenham(x = c(df[df$id == id & df$Ordnung == ordnung & df$vertex_index == 1,]$column,
                                      df[df$id == id & df$Ordnung == ordnung & df$vertex_index == 0,]$column),
                                y = c(df[df$id == id & df$Ordnung == ordnung & df$vertex_index == 1,]$row,
                                      df[df$id == id & df$Ordnung == ordnung & df$vertex_index == 0,]$row))
      for (index in (1:max(df[df$id == id & df$Ordnung == ordnung,]$vertex_index))) {
        partDrainage2 <- bresenham(x = c(df[df$id == id,]$column[index+1], df[df$id == id,]$column[index]),
                                   y = c(df[df$id == id,]$row[index+1], df[df$id == id,]$row[index]))
        
        partDrainage$x <- append(partDrainage2$x, partDrainage$x)
        partDrainage$y <- append(partDrainage2$y, partDrainage$y)
      }
      grabenList[[zähler]] <- as.data.frame(partDrainage)
      zähler <- (zähler + 1)
    }
  }
}

grabenNummer <- 1
for (grabenNummer in 1:length(grabenList)) {
  print(grabenNummer)
  for (i in length(grabenList[[grabenNummer]]$x):1) {
    print(i)
    h8 <- c(dem[grabenList[[grabenNummer]]$y[i]-1, grabenList[[grabenNummer]]$x[i]-1],
            dem[grabenList[[grabenNummer]]$y[i]-1, grabenList[[grabenNummer]]$x[i]],
            dem[grabenList[[grabenNummer]]$y[i]-1, grabenList[[grabenNummer]]$x[i]+1],
            dem[grabenList[[grabenNummer]]$y[i], grabenList[[grabenNummer]]$x[i]-1],
            dem[grabenList[[grabenNummer]]$y[i], grabenList[[grabenNummer]]$x[i]+1],
            dem[grabenList[[grabenNummer]]$y[i]+1, grabenList[[grabenNummer]]$x[i]-1],
            dem[grabenList[[grabenNummer]]$y[i]+1, grabenList[[grabenNummer]]$x[i]],
            dem[grabenList[[grabenNummer]]$y[i]+1, grabenList[[grabenNummer]]$x[i]+1])
    min(h8)
    if (min(h8) > dem[grabenList[[grabenNummer]]$y[i], grabenList[[grabenNummer]]$x[i]]) {
      dem[grabenList[[grabenNummer]]$y[i-1], grabenList[[grabenNummer]]$x[i-1]] <- 
        (dem[grabenList[[grabenNummer]]$y[i], grabenList[[grabenNummer]]$x[i]] - 0.005)
    }
    if (min(h8) < dem[grabenList[[grabenNummer]]$y[i], grabenList[[grabenNummer]]$x[i]]) {
      dem[grabenList[[grabenNummer]]$y[i-1], grabenList[[grabenNummer]]$x[i-1]] <- 
        (min(h8) - 0.01)
    }
    if (min(h8) == dem[grabenList[[grabenNummer]]$y[i], grabenList[[grabenNummer]]$x[i]]) {
      dem[grabenList[[grabenNummer]]$y[i-1], grabenList[[grabenNummer]]$x[i-1]] <- 
        (dem[grabenList[[grabenNummer]]$y[i], grabenList[[grabenNummer]]$x[i]] - 0.005)
    }
  }
}



writeRaster(dem,
            "S:/DGM_Flow_Modellierung/07_data_preprocessing/QGIS_Projekte/Argenschwang_altered_interpol.tif",
            options = c("COMPRESS=NONE"),
            overwrite = TRUE,
            datatype = "FLT8S")

temp <- raster("S:/DGM_Flow_Modellierung/07_data_preprocessing/QGIS_Projekte/Argenschwang_altered_interpol.tif", crs = '+init=EPSG:32632')


#####################################
#Gräbenverschlüsse einprägen - Alle Verschlüsse
dem <- raster("S:/DGM_Flow_Modellierung/07_data_preprocessing/QGIS_Projekte/Argenschwang_altered_interpol.tif", crs = '+init=EPSG:32632')

df <- read_csv("S:/DGM_Flow_Modellierung/07_data_preprocessing/QGIS_Projekte/Grabenverschlüsse_Vertices_Koordinaten.csv")


df <- df[,c("id", "Type", "vertex_index", "xcoord", "ycoord")]

if (Auflösung == 0.5) {
  x_min <- 401384.25
  x_max <- 402434.25
  y_min <- 5530950.75
  y_max <- 5532365.75
  
  #df$xcoord <- round_any(df$xcoord - x_min, 0.25)
  
  df$xcoord <- df$xcoord - x_min
  df$xcoord <- floor(df$xcoord / 0.5) + 1
  
  df$ycoord <- df$ycoord - y_min
  df$ycoord <- floor(df$ycoord / 0.5) + 1
  df$ycoord <- 2830 - df$ycoord
  df$ycoord <- df$ycoord + 1
}

colnames(df) <- c("id", "Type", "vertex_index", "column", "row")

df <- df[order(df$id),]
df$fid <- 1:nrow(df)
df$height = NA

for (fid in df$fid) {
  df[df$fid == fid,]$height <- dem[df[df$fid == fid,]$row, df[df$fid == fid,]$column]
}

zähler <- 1
grabenList <- c()
id <- 1
for (id in sort(unique(df$id))) {
  partDrainage <- bresenham(x = c(df[df$id == id & df$vertex_index == 1,]$column,
                                  df[df$id == id & df$vertex_index == 0,]$column),
                            y = c(df[df$id == id & df$vertex_index == 1,]$row,
                                  df[df$id == id & df$vertex_index == 0,]$row))
  grabenList[[zähler]] <- as.data.frame(partDrainage)
  zähler <- (zähler + 1)
}

for (graben in 1:length(grabenList)) {
  for (zelle in 1:nrow(grabenList[[graben]])) {
    dem[grabenList[[graben]][zelle, "y"], grabenList[[graben]][zelle, "x"]] <- max(df[df$id == graben, "height"])
  }
}

writeRaster(dem,
            "S:/DGM_Flow_Modellierung/07_data_preprocessing/QGIS_Projekte/Argenschwang_altered_Grabenverschluss.tif",
            options = c("COMPRESS=NONE"),
            overwrite = TRUE,
            datatype = "FLT8S")

#####################################
#Gräbenverschlüsse einprägen - Keyline-Verschlüsse
dem <- raster("S:/DGM_Flow_Modellierung/07_data_preprocessing/QGIS_Projekte/Argenschwang_altered_interpol.tif", crs = '+init=EPSG:32632')

df <- read_csv("S:/DGM_Flow_Modellierung/07_data_preprocessing/QGIS_Projekte/Grabenverschlüsse_Vertices_Koordinaten.csv")


df <- df[,c("id", "Type", "vertex_index", "xcoord", "ycoord")]

if (Auflösung == 0.5) {
  x_min <- 401384.25
  x_max <- 402434.25
  y_min <- 5530950.75
  y_max <- 5532365.75
  
  #df$xcoord <- round_any(df$xcoord - x_min, 0.25)
  
  df$xcoord <- df$xcoord - x_min
  df$xcoord <- floor(df$xcoord / 0.5) + 1
  
  df$ycoord <- df$ycoord - y_min
  df$ycoord <- floor(df$ycoord / 0.5) + 1
  df$ycoord <- 2830 - df$ycoord
  df$ycoord <- df$ycoord + 1
}

colnames(df) <- c("id", "Type", "vertex_index", "column", "row")

df <- df[order(df$id),]
df$fid <- 1:nrow(df)
df$height = NA

for (fid in df$fid) {
  df[df$fid == fid,]$height <- dem[df[df$fid == fid,]$row, df[df$fid == fid,]$column]
}

df <- df[df$Type == "Keyline",]

zähler <- 1
grabenList <- c()
id <- 1
for (id in sort(unique(df$id))) {
  partDrainage <- bresenham(x = c(df[df$id == id & df$vertex_index == 1,]$column,
                                  df[df$id == id & df$vertex_index == 0,]$column),
                            y = c(df[df$id == id & df$vertex_index == 1,]$row,
                                  df[df$id == id & df$vertex_index == 0,]$row))
  grabenList[[zähler]] <- as.data.frame(partDrainage)
  zähler <- (zähler + 1)
}

for (graben in 1:length(grabenList)) {
  for (zelle in 1:nrow(grabenList[[graben]])) {
    dem[grabenList[[graben]][zelle, "y"], grabenList[[graben]][zelle, "x"]] <- min(df[df$id == graben, "height"])
  }
}

writeRaster(dem,
            "S:/DGM_Flow_Modellierung/07_data_preprocessing/QGIS_Projekte/Argenschwang_altered_Grabenverschluss_für_Keylines.tif",
            options = c("COMPRESS=NONE"),
            overwrite = TRUE,
            datatype = "FLT8S")


bresenham <- function(x, y = NULL, close = TRUE)
{
  # accept any coordinate structure
  v <- xy.coords(x = x, y = y, recycle = TRUE, setLab = FALSE)
  if (!all(is.finite(v$x), is.finite(v$y)))
    stop("finite coordinates required")
  
  v[1:2] <- lapply(v[1:2], round) # Bresenham's algorithm IS for integers
  nx <- length(v$x)
  if (nx == 1) return(list(x = v$x, y = v$y)) # just one point
  if (nx > 2 && close == TRUE) { # close polygon by replicating 1st point
    v$x <- c(v$x, v$x[1])
    v$y <- c(v$y, v$y[1])
    nx <- nx + 1
  }
  # collect result in 'ans, staring with 1st point
  ans <- lapply(v[1:2], "[", 1)
  
  # process all vertices in pairs
  for (i in seq.int(nx - 1)) {
    x <- v$x[i] # coordinates updated in x, y
    y <- v$y[i]
    x.end <- v$x[i + 1]
    y.end <- v$y[i + 1]
    
    dx <- abs(x.end - x); dy <- -abs(y.end - y)
    sx <- ifelse(x < x.end, 1, -1)
    sy <- ifelse(y < y.end, 1, -1)
    err <- dx + dy
    
    # process one segment
    while(!(isTRUE(all.equal(x, x.end)) && isTRUE(all.equal(y, y.end)))) {
      e2 <- 2 * err
      if (e2 >= dy) { # increment x
        err <- err + dy
        x <- x + sx
      }
      if (e2 <= dx) { # increment y
        err <- err + dx
        y <- y + sy
      }
      ans$x <- c(ans$x, x)
      ans$y <- c(ans$y, y)
    }
  }
  # remove duplicated points (typically 1st and last)
  dups <- duplicated(do.call(cbind, ans), MARGIN = 1) 
  return(lapply(ans, "[", !dups))
}


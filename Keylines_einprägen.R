library("raster")
library("sp")
library(readr)
library(plyr)
library(dplyr)
library(tidyverse)
#library("versus")


#####################################
#Gräbenverschlüsse einprägen - Keyline-Verschlüsse
dem <- raster("S:/DGM_Flow_Modellierung/07_data_preprocessing/QGIS_Projekte/Argenschwang_altered_interpol.tif", crs = '+init=EPSG:32632')

dx <- read_csv("S:/DGM_Flow_Modellierung/07_data_preprocessing/QGIS_Projekte/Grabenverschlüsse_Vertices_Koordinaten.csv")


dx <- dx[,c("id", "Type", "vertex_index", "xcoord", "ycoord")]

if (Auflösung == 0.5) {
  x_min <- 401384.25
  x_max <- 402434.25
  y_min <- 5530950.75
  y_max <- 5532365.75
  
  #df$xcoord <- round_any(df$xcoord - x_min, 0.25)
  
  dx$xcoord <- dx$xcoord - x_min
  dx$xcoord <- floor(dx$xcoord / 0.5) + 1
  
  dx$ycoord <- dx$ycoord - y_min
  dx$ycoord <- floor(dx$ycoord / 0.5) + 1
  dx$ycoord <- 2830 - dx$ycoord
  dx$ycoord <- dx$ycoord + 1
}

colnames(dx) <- c("id", "Type", "vertex_index", "column", "row")

dx <- dx[order(dx$id),]
dx$fid <- 1:nrow(dx)
dx$height = NA

for (fid in dx$fid) {
  dx[dx$fid == fid,]$height <- dem[dx[dx$fid == fid,]$row, dx[dx$fid == fid,]$column]
}

dx <- dx[dx$Type == "Keyline",]

zähler <- 1
verschlussList <- c()
id <- 1
for (id in sort(unique(dx$id))) {
  partDrainage <- bresenham(x = c(dx[dx$id == id & dx$vertex_index == 1,]$column,
                                  dx[dx$id == id & dx$vertex_index == 0,]$column),
                            y = c(dx[dx$id == id & dx$vertex_index == 1,]$row,
                                  dx[dx$id == id & dx$vertex_index == 0,]$row))
  verschlussList[[zähler]] <- as.data.frame(partDrainage)
  zähler <- (zähler + 1)
}

######################

df <- read_csv("S:/DGM_Flow_Modellierung/07_data_preprocessing/QGIS_Projekte/Keylines_Vertices_beidseitig_Koordinaten.csv")

df <- df[,c("id", "Type", "Isometrie", "vertex_index", "xcoord", "ycoord")]


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

colnames(df) <- c("id", "Type", "Isometrie", "vertex_index", "column", "row")

df <- df[order(df$id),]
df$fid <- 1:nrow(df)

### Liste mit Koordinaten der Akkumulationsgräben erstellen

zähler <- 1
grabenList <- c()
id <- 1

for (id in sort(unique(df$id))) {
  partDrainage <- bresenham(x = c(df[df$id == id & df$Type == "Akkumulator" & df$vertex_index == 1,]$column,
                                  df[df$id == id & df$Type == "Akkumulator" & df$vertex_index == 0,]$column),
                            y = c(df[df$id == id & df$Type == "Akkumulator" & df$vertex_index == 1,]$row,
                                  df[df$id == id & df$Type == "Akkumulator" & df$vertex_index == 0,]$row))
  for (index in (1:max(df[df$id == id & df$Type == "Akkumulator",]$vertex_index))) {
    partDrainage2 <- bresenham(x = c(df[df$id == id & df$Type == "Akkumulator" & df$vertex_index == (index + 1),]$column, df[df$id == id & df$Type == "Akkumulator" & df$vertex_index == index,]$column),
                               y = c(df[df$id == id & df$Type == "Akkumulator" & df$vertex_index == (index + 1),]$row, df[df$id == id & df$Type == "Akkumulator" & df$vertex_index == index,]$row))
    
    partDrainage$x <- append(head(partDrainage2$x, -1), partDrainage$x)
    partDrainage$y <- append(head(partDrainage2$y, -1), partDrainage$y)
  }
  grabenList[[zähler]] <- as.data.frame(partDrainage)
  zähler <- (zähler + 1)
}

### Akkumulationsgräben einprägen

for (graben in 1:length(grabenList)) {
  h1 <- dem[grabenList[[graben]][length(grabenList[[graben]]$y),2], grabenList[[graben]][length(grabenList[[graben]]$x),1]]
  dem[grabenList[[graben]][1,2], grabenList[[graben]][1,1]] <- (dem[grabenList[[graben]][1,2], grabenList[[graben]][1,1]] - 0.3)
  h2 <- dem[grabenList[[graben]][1,2], grabenList[[graben]][1,1]]
  
  
  delta_h <- h1 - h2
  
  n <- length(grabenList[[graben]]$y) - 1
  for (i in 2:length(grabenList[[graben]]$y)) {
    dem[grabenList[[graben]][i,2], grabenList[[graben]][i,1]] <- 
      (dem[grabenList[[graben]][1,2], grabenList[[graben]][1,1]] + (delta_h / n) * (i-1))
    print(dem[grabenList[[graben]][i,2], grabenList[[graben]][i,1]])
  }
}

### Grabenverschlüsse einprägen

for (graben in 1:length(verschlussList)) {
  for (zelle in 1:nrow(verschlussList[[graben]])) {
    dem[verschlussList[[graben]][zelle, "y"], verschlussList[[graben]][zelle, "x"]] <- max(dx[dx$id == graben, "height"])
  }
}

### Distributionsgräben definieren


zähler <- 1
grabenList <- c()
id <- 1

for (id in sort(unique(df$id))) {
  partDrainage <- bresenham(x = c(df[df$id == id & df$Type == "Distributor" & df$vertex_index == 1,]$column,
                                  df[df$id == id & df$Type == "Distributor" & df$vertex_index == 0,]$column),
                            y = c(df[df$id == id & df$Type == "Distributor" & df$vertex_index == 1,]$row,
                                  df[df$id == id & df$Type == "Distributor" & df$vertex_index == 0,]$row))
  for (index in (1:max(df[df$id == id & df$Type == "Distributor",]$vertex_index))) {
    partDrainage2 <- bresenham(x = c(df[df$id == id & df$Type == "Distributor" & df$vertex_index == (index + 1),]$column,
                                     df[df$id == id & df$Type == "Distributor" & df$vertex_index == index,]$column),
                               y = c(df[df$id == id & df$Type == "Distributor" & df$vertex_index == (index + 1),]$row,
                                     df[df$id == id & df$Type == "Distributor" & df$vertex_index == index,]$row))
    
    partDrainage$x <- append(head(partDrainage2$x, -1), partDrainage$x)
    partDrainage$y <- append(head(partDrainage2$y, -1), partDrainage$y)
  }
  grabenList[[zähler]] <- as.data.frame(partDrainage)
  zähler <- (zähler + 1)
}

### Distributorgräben einprägen

grabenTiefen <- c()
for (graben in 1:length(grabenList)) {
  
  h <- dem[grabenList[[graben]][length(grabenList[[graben]]$y),2], grabenList[[graben]][length(grabenList[[graben]]$x),1]]
  
  for (i in 2:length(grabenList[[graben]]$y)) {
    dem[grabenList[[graben]][i,2], grabenList[[graben]][i,1]] <- h
    print(dem[grabenList[[graben]][i,2], grabenList[[graben]][i,1]])
  }
  grabenTiefen[graben] <- h
}

### Liste mit Koordinaten der Endgräben erstellen

zähler <- 1
grabenList <- c()
id <- 1

for (id in sort(unique(df$id))) {
  partDrainage <- bresenham(x = c(df[df$id == id & df$Type == "Endgraben" & df$vertex_index == 1,]$column,
                                  df[df$id == id & df$Type == "Endgraben" & df$vertex_index == 0,]$column),
                            y = c(df[df$id == id & df$Type == "Endgraben" & df$vertex_index == 1,]$row,
                                  df[df$id == id & df$Type == "Endgraben" & df$vertex_index == 0,]$row))
  for (index in (1:max(df[df$id == id & df$Type == "Endgraben",]$vertex_index))) {
    partDrainage2 <- bresenham(x = c(df[df$id == id & df$Type == "Endgraben" & df$vertex_index == (index + 1),]$column, df[df$id == id & df$Type == "Endgraben" & df$vertex_index == index,]$column),
                               y = c(df[df$id == id & df$Type == "Endgraben" & df$vertex_index == (index + 1),]$row, df[df$id == id & df$Type == "Endgraben" & df$vertex_index == index,]$row))
    
    partDrainage$x <- append(head(partDrainage2$x, -1), partDrainage$x)
    partDrainage$y <- append(head(partDrainage2$y, -1), partDrainage$y)
  }
  grabenList[[zähler]] <- as.data.frame(partDrainage)
  zähler <- (zähler + 1)
}

### Endgräben einprägen
#graben <- 1
graben <- 1
for (graben in 1:length(grabenList)) {
  h1 <- dem[grabenList[[graben]][length(grabenList[[graben]]$y),2], grabenList[[graben]][length(grabenList[[graben]]$x),1]]
  #h2 <- dem[grabenList[[graben]][1,2], grabenList[[graben]][1,1]]
  h2 <- grabenTiefen[graben]
  
  delta_h <- h1 - h2
  
  n <- length(grabenList[[graben]]$y) - 1
  
  dem[grabenList[[graben]][1,2], grabenList[[graben]][1,1]] <- h2
  
  for (i in 2:length(grabenList[[graben]]$y)) {
    dem[grabenList[[graben]][i,2], grabenList[[graben]][i,1]] <- 
      (h2 + (delta_h / n) * (i-1))
    print(dem[grabenList[[graben]][i,2], grabenList[[graben]][i,1]])
  }
}

### Liste mit Koordinaten der Endgräben erstellen - Für die Wälle

zähler <- 1
grabenList <- c()
id <- 1

for (id in sort(unique(df[df$Type == "Endgraben" & df$Isometrie == "Links",]$id))) {
  partDrainage <- bresenham(x = c(df[df$id == id & df$Type == "Endgraben" & df$vertex_index == 1,]$column,
                                  df[df$id == id & df$Type == "Endgraben" & df$vertex_index == 0,]$column),
                            y = c(df[df$id == id & df$Type == "Endgraben" & df$vertex_index == 1,]$row,
                                  df[df$id == id & df$Type == "Endgraben" & df$vertex_index == 0,]$row))
  for (index in (1:max(df[df$id == id & df$Type == "Endgraben",]$vertex_index))) {
    partDrainage2 <- bresenham(x = c(df[df$id == id & df$Type == "Endgraben" & df$vertex_index == (index + 1),]$column, df[df$id == id & df$Type == "Endgraben" & df$vertex_index == index,]$column),
                               y = c(df[df$id == id & df$Type == "Endgraben" & df$vertex_index == (index + 1),]$row, df[df$id == id & df$Type == "Endgraben" & df$vertex_index == index,]$row))
    
    partDrainage$x <- append(head(partDrainage2$x, -1), partDrainage$x)
    partDrainage$y <- append(head(partDrainage2$y, -1), partDrainage$y)
  }
  grabenList[[zähler]] <- as.data.frame(partDrainage)
  zähler <- (zähler + 1)
}

### Wälle aufschütten

k <- 1
for (k in 1:length(grabenList)) {
  
  zellen_Wall <- grabenList[[k]][length(grabenList[[k]]$x),]
  zellen_temp <- zellen_Wall 
  
  for (i in length(grabenList[[k]]$x):2) {
    
    aktuelle_Zelle <- grabenList[[k]][i,]
    folgende_Zelle <- grabenList[[k]][i-1,]
    
    zellen_Liste <- c()
    
    zelle_Diff <- folgende_Zelle - aktuelle_Zelle
    
    if(zelle_Diff$x == 1 & zelle_Diff$y == -1) {
      zellen_Liste <- c(1, 2, 3)
    }
    if(zelle_Diff$x == 0 & zelle_Diff$y == -1) {
      zellen_Liste <- c(2, 3, 4)
    }
    if(zelle_Diff$x == -1 & zelle_Diff$y == -1) {
      zellen_Liste <- c(3, 4 , 5)
    }
    if(zelle_Diff$x == -1 & zelle_Diff$y == 0) {
      zellen_Liste <- c(4, 5, 6)
    }
    if(zelle_Diff$x == -1 & zelle_Diff$y == 1) {
      zellen_Liste <- c(5, 6, 7)
    }
    if(zelle_Diff$x == 0 & zelle_Diff$y == 1) {
      zellen_Liste <- c(6, 7, 8)
    }
    if(zelle_Diff$x == 1 & zelle_Diff$y == 1) {
      zellen_Liste <- c(7, 8, 1)
    }
    if(zelle_Diff$x == 1 & zelle_Diff$y == 0) {
      zellen_Liste <- c(8, 1, 2)
    }
    
    if(1 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x, aktuelle_Zelle$y - 1))
    }
    if(2 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x - 1, aktuelle_Zelle$y - 1))
    }
    if(3 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x - 1, aktuelle_Zelle$y))
    }
    if(4 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x - 1, aktuelle_Zelle$y + 1))
    }
    if(5 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x, aktuelle_Zelle$y + 1))
    }
    if(6 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x + 1, aktuelle_Zelle$y + 1))
    }
    if(7 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x + 1, aktuelle_Zelle$y))
    }
    if(8 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x + 1, aktuelle_Zelle$y - 1))
    }
    
  }
  zellen_Wall <- zellen_Wall[-1,]
  zellen_Wall <- unique(zellen_Wall)
  
  
  
  zellen_Wall <- anti_join(zellen_Wall, grabenList[[k]])
  
  for (i in 1:nrow(zellen_Wall)) {
    dem[zellen_Wall$y[i], zellen_Wall$x[i]] <- (dem[zellen_Wall$y[i], zellen_Wall$x[i]] + 0.3)
  }
}

### Liste mit Koordinaten der Endgräben erstellen - Für die Wälle

zähler <- 1
grabenList <- c()
id <- 1

for (id in sort(unique(df[df$Type == "Endgraben" & df$Isometrie == "Rechts",]$id))) {
  partDrainage <- bresenham(x = c(df[df$id == id & df$Type == "Endgraben" & df$vertex_index == 1,]$column,
                                  df[df$id == id & df$Type == "Endgraben" & df$vertex_index == 0,]$column),
                            y = c(df[df$id == id & df$Type == "Endgraben" & df$vertex_index == 1,]$row,
                                  df[df$id == id & df$Type == "Endgraben" & df$vertex_index == 0,]$row))
  for (index in (1:max(df[df$id == id & df$Type == "Endgraben",]$vertex_index))) {
    partDrainage2 <- bresenham(x = c(df[df$id == id & df$Type == "Endgraben" & df$vertex_index == (index + 1),]$column, df[df$id == id & df$Type == "Endgraben" & df$vertex_index == index,]$column),
                               y = c(df[df$id == id & df$Type == "Endgraben" & df$vertex_index == (index + 1),]$row, df[df$id == id & df$Type == "Endgraben" & df$vertex_index == index,]$row))
    
    partDrainage$x <- append(head(partDrainage2$x, -1), partDrainage$x)
    partDrainage$y <- append(head(partDrainage2$y, -1), partDrainage$y)
  }
  grabenList[[zähler]] <- as.data.frame(partDrainage)
  zähler <- (zähler + 1)
}
### Wälle aufschütten

k <- 1
for (k in 1:length(grabenList)) {
  
  zellen_Wall <- grabenList[[k]][length(grabenList[[k]]$x),]
  zellen_temp <- zellen_Wall 
  
  for (i in length(grabenList[[k]]$x):2) {
    
    aktuelle_Zelle <- grabenList[[k]][i,]
    folgende_Zelle <- grabenList[[k]][i-1,]
    
    zellen_Liste <- c()
    
    zelle_Diff <- folgende_Zelle - aktuelle_Zelle
    
    if(zelle_Diff$x == 1 & zelle_Diff$y == -1) {
      zellen_Liste <- c(7, 6, 5)
    }
    if(zelle_Diff$x == 0 & zelle_Diff$y == -1) {
      zellen_Liste <- c(8, 7, 6)
    }
    if(zelle_Diff$x == -1 & zelle_Diff$y == -1) {
      zellen_Liste <- c(1, 8, 7)
    }
    if(zelle_Diff$x == -1 & zelle_Diff$y == 0) {
      zellen_Liste <- c(2, 1, 8)
    }
    if(zelle_Diff$x == -1 & zelle_Diff$y == 1) {
      zellen_Liste <- c(3, 2, 1)
    }
    if(zelle_Diff$x == 0 & zelle_Diff$y == 1) {
      zellen_Liste <- c(2, 3, 4)
    }
    if(zelle_Diff$x == 1 & zelle_Diff$y == 1) {
      zellen_Liste <- c(3, 4, 5)
    }
    if(zelle_Diff$x == 1 & zelle_Diff$y == 0) {
      zellen_Liste <- c(4, 5, 6)
    }
    
    if(1 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x, aktuelle_Zelle$y - 1))
    }
    if(2 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x - 1, aktuelle_Zelle$y - 1))
    }
    if(3 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x - 1, aktuelle_Zelle$y))
    }
    if(4 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x - 1, aktuelle_Zelle$y + 1))
    }
    if(5 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x, aktuelle_Zelle$y + 1))
    }
    if(6 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x + 1, aktuelle_Zelle$y + 1))
    }
    if(7 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x + 1, aktuelle_Zelle$y))
    }
    if(8 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x + 1, aktuelle_Zelle$y - 1))
    }
    
  }
  zellen_Wall <- zellen_Wall[-1,]
  zellen_Wall <- unique(zellen_Wall)
  
  
  
  zellen_Wall <- anti_join(zellen_Wall, grabenList[[k]])
  
  for (i in 1:nrow(zellen_Wall)) {
    dem[zellen_Wall$y[i], zellen_Wall$x[i]] <- (dem[zellen_Wall$y[i], zellen_Wall$x[i]] + 0.3)
    
    #dem[zellen_Wall$y[i], zellen_Wall$x[i]] <- (dem[grabenList[[k]]$y[i], grabenList[[k]]$x[i]] + 0.30)
  }
}

### Liste mit Koordinaten der Akkumulationsgräben erstellen - FÜR DIE WÄLLE

zähler <- 1
grabenList <- c()
id <- 1

for (id in sort(unique(df[df$Isometrie == "Links" & df$Type == "Akkumulator", ]$id))) {
  partDrainage <- bresenham(x = c(df[df$id == id & df$Type == "Akkumulator" & df$vertex_index == 1,]$column,
                                  df[df$id == id & df$Type == "Akkumulator" & df$vertex_index == 0,]$column),
                            y = c(df[df$id == id & df$Type == "Akkumulator" & df$vertex_index == 1,]$row,
                                  df[df$id == id & df$Type == "Akkumulator" & df$vertex_index == 0,]$row))
  for (index in (1:max(df[df$id == id & df$Type == "Akkumulator",]$vertex_index))) {
    partDrainage2 <- bresenham(x = c(df[df$id == id & df$Type == "Akkumulator" & df$vertex_index == (index + 1),]$column, df[df$id == id & df$Type == "Akkumulator" & df$vertex_index == index,]$column),
                               y = c(df[df$id == id & df$Type == "Akkumulator" & df$vertex_index == (index + 1),]$row, df[df$id == id & df$Type == "Akkumulator" & df$vertex_index == index,]$row))
    
    partDrainage$x <- append(head(partDrainage2$x, -1), partDrainage$x)
    partDrainage$y <- append(head(partDrainage2$y, -1), partDrainage$y)
  }
  grabenList[[zähler]] <- as.data.frame(partDrainage)
  zähler <- (zähler + 1)
}

### Wälle aufschütten

k <- 1
for (k in 1:length(grabenList)) {
  
  zellen_Wall <- grabenList[[k]][length(grabenList[[k]]$x),]
  zellen_temp <- zellen_Wall 
  
  for (i in length(grabenList[[k]]$x):2) {
    
    aktuelle_Zelle <- grabenList[[k]][i,]
    folgende_Zelle <- grabenList[[k]][i-1,]
    
    zellen_Liste <- c()
    
    zelle_Diff <- folgende_Zelle - aktuelle_Zelle
    
    if(zelle_Diff$x == 1 & zelle_Diff$y == -1) {
      zellen_Liste <- c(1, 2, 3)
    }
    if(zelle_Diff$x == 0 & zelle_Diff$y == -1) {
      zellen_Liste <- c(2, 3, 4)
    }
    if(zelle_Diff$x == -1 & zelle_Diff$y == -1) {
      zellen_Liste <- c(3, 4 , 5)
    }
    if(zelle_Diff$x == -1 & zelle_Diff$y == 0) {
      zellen_Liste <- c(4, 5, 6)
    }
    if(zelle_Diff$x == -1 & zelle_Diff$y == 1) {
      zellen_Liste <- c(5, 6, 7)
    }
    if(zelle_Diff$x == 0 & zelle_Diff$y == 1) {
      zellen_Liste <- c(6, 7, 8)
    }
    if(zelle_Diff$x == 1 & zelle_Diff$y == 1) {
      zellen_Liste <- c(7, 8, 1)
    }
    if(zelle_Diff$x == 1 & zelle_Diff$y == 0) {
      zellen_Liste <- c(8, 1, 2)
    }
    
    if(1 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x, aktuelle_Zelle$y - 1))
    }
    if(2 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x - 1, aktuelle_Zelle$y - 1))
    }
    if(3 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x - 1, aktuelle_Zelle$y))
    }
    if(4 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x - 1, aktuelle_Zelle$y + 1))
    }
    if(5 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x, aktuelle_Zelle$y + 1))
    }
    if(6 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x + 1, aktuelle_Zelle$y + 1))
    }
    if(7 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x + 1, aktuelle_Zelle$y))
    }
    if(8 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x + 1, aktuelle_Zelle$y - 1))
    }
    
  }
  zellen_Wall <- zellen_Wall[-1,]
  zellen_Wall <- unique(zellen_Wall)
  
  zellen_Wall <- anti_join(zellen_Wall, grabenList[[k]])
  
  for (i in 1:nrow(zellen_Wall)) {
    dem[zellen_Wall$y[i], zellen_Wall$x[i]] <- (dem[zellen_Wall$y[i], zellen_Wall$x[i]] + 0.3)
  }
}


### Distributionsgräben definieren - FÜR WÄLLE

zähler <- 1
grabenList <- c()
id <- 1

for (id in sort(unique(df[df$Isometrie == "Links" & df$Type == "Distributor", ]$id))) {
  partDrainage <- bresenham(x = c(df[df$id == id & df$Type == "Distributor" & df$vertex_index == 1,]$column,
                                  df[df$id == id & df$Type == "Distributor" & df$vertex_index == 0,]$column),
                            y = c(df[df$id == id & df$Type == "Distributor" & df$vertex_index == 1,]$row,
                                  df[df$id == id & df$Type == "Distributor" & df$vertex_index == 0,]$row))
  for (index in (1:max(df[df$id == id & df$Type == "Distributor",]$vertex_index))) {
    partDrainage2 <- bresenham(x = c(df[df$id == id & df$Type == "Distributor" & df$vertex_index == (index + 1),]$column,
                                     df[df$id == id & df$Type == "Distributor" & df$vertex_index == index,]$column),
                               y = c(df[df$id == id & df$Type == "Distributor" & df$vertex_index == (index + 1),]$row,
                                     df[df$id == id & df$Type == "Distributor" & df$vertex_index == index,]$row))
    
    partDrainage$x <- append(head(partDrainage2$x, -1), partDrainage$x)
    partDrainage$y <- append(head(partDrainage2$y, -1), partDrainage$y)
  }
  grabenList[[zähler]] <- as.data.frame(partDrainage)
  zähler <- (zähler + 1)
}

### Wälle aufschütten
k <- 1
for (k in 1:length(grabenList)) {
  
  zellen_Wall <- grabenList[[k]][length(grabenList[[k]]$x),]
  zellen_temp <- zellen_Wall 
  
  for (i in length(grabenList[[k]]$x):2) {
    
    aktuelle_Zelle <- grabenList[[k]][i,]
    folgende_Zelle <- grabenList[[k]][i-1,]
    
    zellen_Liste <- c()
    
    zelle_Diff <- folgende_Zelle - aktuelle_Zelle
    
    if(zelle_Diff$x == 1 & zelle_Diff$y == -1) {
      zellen_Liste <- c(1, 2, 3)
    }
    if(zelle_Diff$x == 0 & zelle_Diff$y == -1) {
      zellen_Liste <- c(2, 3, 4)
    }
    if(zelle_Diff$x == -1 & zelle_Diff$y == -1) {
      zellen_Liste <- c(3, 4, 5)
    }
    if(zelle_Diff$x == -1 & zelle_Diff$y == 0) {
      zellen_Liste <- c(4, 5, 6)
    }
    if(zelle_Diff$x == -1 & zelle_Diff$y == 1) {
      zellen_Liste <- c(5, 6, 7)
    }
    if(zelle_Diff$x == 0 & zelle_Diff$y == 1) {
      zellen_Liste <- c(6, 7, 8)
    }
    if(zelle_Diff$x == 1 & zelle_Diff$y == 1) {
      zellen_Liste <- c(7, 8, 1)
    }
    if(zelle_Diff$x == 1 & zelle_Diff$y == 0) {
      zellen_Liste <- c(8, 1, 2)
    }
    
    if(1 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x, aktuelle_Zelle$y - 1))
    }
    if(2 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x - 1, aktuelle_Zelle$y - 1))
    }
    if(3 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x - 1, aktuelle_Zelle$y))
    }
    if(4 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x - 1, aktuelle_Zelle$y + 1))
    }
    if(5 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x, aktuelle_Zelle$y + 1))
    }
    if(6 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x + 1, aktuelle_Zelle$y + 1))
    }
    if(7 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x + 1, aktuelle_Zelle$y))
    }
    if(8 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x + 1, aktuelle_Zelle$y - 1))
    }
    
    #Endzellen definieren
    if (i == 2) {
      zellen_Ende <- folgende_Zelle
      zellen_Ende <- rbind(zellen_Ende, list(folgende_Zelle$x - 1, folgende_Zelle$y - 1))
      zellen_Ende <- rbind(zellen_Ende, list(folgende_Zelle$x, folgende_Zelle$y - 1))
      zellen_Ende <- rbind(zellen_Ende, list(folgende_Zelle$x + 1, folgende_Zelle$y - 1))
      zellen_Ende <- rbind(zellen_Ende, list(folgende_Zelle$x - 1, folgende_Zelle$y))
      zellen_Ende <- rbind(zellen_Ende, list(folgende_Zelle$x + 1, folgende_Zelle$y))
      zellen_Ende <- rbind(zellen_Ende, list(folgende_Zelle$x - 1, folgende_Zelle$y + 1))
      zellen_Ende <- rbind(zellen_Ende, list(folgende_Zelle$x, folgende_Zelle$y + 1))
      zellen_Ende <- rbind(zellen_Ende, list(folgende_Zelle$x + 1, folgende_Zelle$y + 1))
      zellen_Ende <- anti_join(zellen_Ende, grabenList[[k]])
    }
  }
  zellen_Wall <- zellen_Wall[-1,]
  zellen_Wall <- unique(zellen_Wall)
  
  zellen_Wall <- anti_join(zellen_Wall, grabenList[[k]])
  
  for (WallZelle in 1:nrow(zellen_Wall)) {
    dem[zellen_Wall$y[WallZelle], zellen_Wall$x[WallZelle]] <- (dem[grabenList[[k]]$y[2], grabenList[[k]]$x[2]] + 0.35)
  }
}

### Grabenerweitern links
### Wälle aufschütten
k <- 1
for (k in 1:length(grabenList)) {
  
  zellen_Wall <- grabenList[[k]][length(grabenList[[k]]$x),]
  zellen_temp <- zellen_Wall 
  
  for (i in length(grabenList[[k]]$x):2) {
    
    aktuelle_Zelle <- grabenList[[k]][i,]
    folgende_Zelle <- grabenList[[k]][i-1,]
    
    zellen_Liste <- c()
    
    zelle_Diff <- folgende_Zelle - aktuelle_Zelle
    
    if(zelle_Diff$x == 1 & zelle_Diff$y == -1) {
      zellen_Liste <- c(7, 6, 5)
    }
    if(zelle_Diff$x == 0 & zelle_Diff$y == -1) {
      zellen_Liste <- c(8, 7, 6)
    }
    if(zelle_Diff$x == -1 & zelle_Diff$y == -1) {
      zellen_Liste <- c(1, 8, 7)
    }
    if(zelle_Diff$x == -1 & zelle_Diff$y == 0) {
      zellen_Liste <- c(2, 1, 8)
    }
    if(zelle_Diff$x == -1 & zelle_Diff$y == 1) {
      zellen_Liste <- c(3, 2, 1)
    }
    if(zelle_Diff$x == 0 & zelle_Diff$y == 1) {
      zellen_Liste <- c(2, 3, 4)
    }
    if(zelle_Diff$x == 1 & zelle_Diff$y == 1) {
      zellen_Liste <- c(3, 4, 5)
    }
    if(zelle_Diff$x == 1 & zelle_Diff$y == 0) {
      zellen_Liste <- c(4, 5, 6)
    }
    
    if(1 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x, aktuelle_Zelle$y - 1))
    }
    if(2 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x - 1, aktuelle_Zelle$y - 1))
    }
    if(3 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x - 1, aktuelle_Zelle$y))
    }
    if(4 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x - 1, aktuelle_Zelle$y + 1))
    }
    if(5 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x, aktuelle_Zelle$y + 1))
    }
    if(6 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x + 1, aktuelle_Zelle$y + 1))
    }
    if(7 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x + 1, aktuelle_Zelle$y))
    }
    if(8 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x + 1, aktuelle_Zelle$y - 1))
    }
    
  }
  zellen_Wall <- zellen_Wall[-1,]
  zellen_Wall <- unique(zellen_Wall)
  
  
  
  zellen_Wall <- anti_join(zellen_Wall, grabenList[[k]])
  
  for (WallZelle in 1:nrow(zellen_Wall)) {
    dem[zellen_Wall$y[WallZelle], zellen_Wall$x[WallZelle]] <- (dem[grabenList[[k]]$y[2], grabenList[[k]]$x[2]] + 0.1)
  }
  
}

f_dem <- dem
dem <- f_dem
### Liste mit Koordinaten der Akkumulationsgräben erstellen - FÜR DIE WÄLLE - RECHTS

zähler <- 1
grabenList <- c()
id <- 1

for (id in sort(unique(df[df$Isometrie == "Rechts" & df$Type == "Akkumulator", ]$id))) {
  partDrainage <- bresenham(x = c(df[df$id == id & df$Type == "Akkumulator" & df$vertex_index == 1,]$column,
                                  df[df$id == id & df$Type == "Akkumulator" & df$vertex_index == 0,]$column),
                            y = c(df[df$id == id & df$Type == "Akkumulator" & df$vertex_index == 1,]$row,
                                  df[df$id == id & df$Type == "Akkumulator" & df$vertex_index == 0,]$row))
  for (index in (1:max(df[df$id == id & df$Type == "Akkumulator",]$vertex_index))) {
    partDrainage2 <- bresenham(x = c(df[df$id == id & df$Type == "Akkumulator" & df$vertex_index == (index + 1),]$column, df[df$id == id & df$Type == "Akkumulator" & df$vertex_index == index,]$column),
                               y = c(df[df$id == id & df$Type == "Akkumulator" & df$vertex_index == (index + 1),]$row, df[df$id == id & df$Type == "Akkumulator" & df$vertex_index == index,]$row))
    
    partDrainage$x <- append(head(partDrainage2$x, -1), partDrainage$x)
    partDrainage$y <- append(head(partDrainage2$y, -1), partDrainage$y)
  }
  grabenList[[zähler]] <- as.data.frame(partDrainage)
  zähler <- (zähler + 1)
}

### Wälle aufschütten

k <- 1
for (k in 1:length(grabenList)) {
  
  zellen_Wall <- grabenList[[k]][length(grabenList[[k]]$x),]
  zellen_temp <- zellen_Wall 
  
  for (i in length(grabenList[[k]]$x):2) {
    
    aktuelle_Zelle <- grabenList[[k]][i,]
    folgende_Zelle <- grabenList[[k]][i-1,]
    
    zellen_Liste <- c()
    
    zelle_Diff <- folgende_Zelle - aktuelle_Zelle
    
    if(zelle_Diff$x == 1 & zelle_Diff$y == -1) {
      zellen_Liste <- c(7, 6, 5)
    }
    if(zelle_Diff$x == 0 & zelle_Diff$y == -1) {
      zellen_Liste <- c(8, 7, 6)
    }
    if(zelle_Diff$x == -1 & zelle_Diff$y == -1) {
      zellen_Liste <- c(1, 8, 7)
    }
    if(zelle_Diff$x == -1 & zelle_Diff$y == 0) {
      zellen_Liste <- c(2, 1, 8)
    }
    if(zelle_Diff$x == -1 & zelle_Diff$y == 1) {
      zellen_Liste <- c(3, 2, 1)
    }
    if(zelle_Diff$x == 0 & zelle_Diff$y == 1) {
      zellen_Liste <- c(2, 3, 4)
    }
    if(zelle_Diff$x == 1 & zelle_Diff$y == 1) {
      zellen_Liste <- c(3, 4, 5)
    }
    if(zelle_Diff$x == 1 & zelle_Diff$y == 0) {
      zellen_Liste <- c(4, 5, 6)
    }
    
    if(1 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x, aktuelle_Zelle$y - 1))
    }
    if(2 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x - 1, aktuelle_Zelle$y - 1))
    }
    if(3 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x - 1, aktuelle_Zelle$y))
    }
    if(4 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x - 1, aktuelle_Zelle$y + 1))
    }
    if(5 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x, aktuelle_Zelle$y + 1))
    }
    if(6 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x + 1, aktuelle_Zelle$y + 1))
    }
    if(7 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x + 1, aktuelle_Zelle$y))
    }
    if(8 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x + 1, aktuelle_Zelle$y - 1))
    }
    
  }
  zellen_Wall <- zellen_Wall[-1,]
  zellen_Wall <- unique(zellen_Wall)
  
  
  
  zellen_Wall <- anti_join(zellen_Wall, grabenList[[k]])
  
  for (i in 1:nrow(zellen_Wall)) {
    dem[zellen_Wall$y[i], zellen_Wall$x[i]] <- (dem[zellen_Wall$y[i], zellen_Wall$x[i]] + 0.3)
    
    #dem[zellen_Wall$y[i], zellen_Wall$x[i]] <- (dem[grabenList[[k]]$y[i], grabenList[[k]]$x[i]] + 0.30)
  }
}

#######hier

### Distributionsgräben definieren - FÜR WÄLLE - Rechts

zähler <- 1
grabenList <- c()
id <- 1

for (id in sort(unique(df[df$Isometrie == "Rechts" & df$Type == "Distributor", ]$id))) {
  partDrainage <- bresenham(x = c(df[df$id == id & df$Type == "Distributor" & df$vertex_index == 1,]$column,
                                  df[df$id == id & df$Type == "Distributor" & df$vertex_index == 0,]$column),
                            y = c(df[df$id == id & df$Type == "Distributor" & df$vertex_index == 1,]$row,
                                  df[df$id == id & df$Type == "Distributor" & df$vertex_index == 0,]$row))
  for (index in (1:max(df[df$id == id & df$Type == "Distributor",]$vertex_index))) {
    partDrainage2 <- bresenham(x = c(df[df$id == id & df$Type == "Distributor" & df$vertex_index == (index + 1),]$column,
                                     df[df$id == id & df$Type == "Distributor" & df$vertex_index == index,]$column),
                               y = c(df[df$id == id & df$Type == "Distributor" & df$vertex_index == (index + 1),]$row,
                                     df[df$id == id & df$Type == "Distributor" & df$vertex_index == index,]$row))
    
    partDrainage$x <- append(head(partDrainage2$x, -1), partDrainage$x)
    partDrainage$y <- append(head(partDrainage2$y, -1), partDrainage$y)
  }
  grabenList[[zähler]] <- as.data.frame(partDrainage)
  zähler <- (zähler + 1)
}

### Wälle aufschütten - Distribution - Rechts

for (k in 1:length(grabenList)) {
  
  zellen_Wall <- grabenList[[k]][length(grabenList[[k]]$x),]
  zellen_temp <- zellen_Wall 
  
  for (i in length(grabenList[[k]]$x):2) {
    
    aktuelle_Zelle <- grabenList[[k]][i,]
    folgende_Zelle <- grabenList[[k]][i-1,]
    
    zellen_Liste <- c()
    
    zelle_Diff <- folgende_Zelle - aktuelle_Zelle
    
    if(zelle_Diff$x == 1 & zelle_Diff$y == -1) {
      zellen_Liste <- c(7, 6, 5)
    }
    if(zelle_Diff$x == 0 & zelle_Diff$y == -1) {
      zellen_Liste <- c(8, 7, 6)
    }
    if(zelle_Diff$x == -1 & zelle_Diff$y == -1) {
      zellen_Liste <- c(1, 8, 7)
    }
    if(zelle_Diff$x == -1 & zelle_Diff$y == 0) {
      zellen_Liste <- c(2, 1, 8)
    }
    if(zelle_Diff$x == -1 & zelle_Diff$y == 1) {
      zellen_Liste <- c(3, 2, 1)
    }
    if(zelle_Diff$x == 0 & zelle_Diff$y == 1) {
      zellen_Liste <- c(4, 3, 2)
    }
    if(zelle_Diff$x == 1 & zelle_Diff$y == 1) {
      zellen_Liste <- c(5, 4, 3)
    }
    if(zelle_Diff$x == 1 & zelle_Diff$y == 0) {
      zellen_Liste <- c(6, 5, 4)
    }
    
    if(1 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x, aktuelle_Zelle$y - 1))
    }
    if(2 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x - 1, aktuelle_Zelle$y - 1))
    }
    if(3 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x - 1, aktuelle_Zelle$y))
    }
    if(4 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x - 1, aktuelle_Zelle$y + 1))
    }
    if(5 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x, aktuelle_Zelle$y + 1))
    }
    if(6 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x + 1, aktuelle_Zelle$y + 1))
    }
    if(7 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x + 1, aktuelle_Zelle$y))
    }
    if(8 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x + 1, aktuelle_Zelle$y - 1))
    }
  }
  zellen_Wall <- zellen_Wall[-1,]
  zellen_Wall <- unique(zellen_Wall)
  
  
  
  zellen_Wall <- anti_join(zellen_Wall, grabenList[[k]])
  
  for (WallZelle in 1:nrow(zellen_Wall)) {
    dem[zellen_Wall$y[WallZelle], zellen_Wall$x[WallZelle]] <- (dem[grabenList[[k]]$y[2], grabenList[[k]]$x[2]] + 0.35)
  }
}

###Grabenerweiterung links
### Wälle aufschütten
for (k in 1:length(grabenList)) {
  
  zellen_Wall <- grabenList[[k]][length(grabenList[[k]]$x),]
  zellen_temp <- zellen_Wall 
  
  for (i in length(grabenList[[k]]$x):2) {
    
    aktuelle_Zelle <- grabenList[[k]][i,]
    folgende_Zelle <- grabenList[[k]][i-1,]
    
    zellen_Liste <- c()
    
    zelle_Diff <- folgende_Zelle - aktuelle_Zelle
    
    if(zelle_Diff$x == 1 & zelle_Diff$y == -1) {
      zellen_Liste <- c(1, 2, 3)
    }
    if(zelle_Diff$x == 0 & zelle_Diff$y == -1) {
      zellen_Liste <- c(2, 3, 4)
    }
    if(zelle_Diff$x == -1 & zelle_Diff$y == -1) {
      zellen_Liste <- c(3, 4, 5)
    }
    if(zelle_Diff$x == -1 & zelle_Diff$y == 0) {
      zellen_Liste <- c(4, 5, 6)
    }
    if(zelle_Diff$x == -1 & zelle_Diff$y == 1) {
      zellen_Liste <- c(5, 6, 7)
    }
    if(zelle_Diff$x == 0 & zelle_Diff$y == 1) {
      zellen_Liste <- c(6, 7, 8)
    }
    if(zelle_Diff$x == 1 & zelle_Diff$y == 1) {
      zellen_Liste <- c(7, 8, 1)
    }
    if(zelle_Diff$x == 1 & zelle_Diff$y == 0) {
      zellen_Liste <- c(8, 1, 2)
    }
    
    if(1 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x, aktuelle_Zelle$y - 1))
    }
    if(2 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x - 1, aktuelle_Zelle$y - 1))
    }
    if(3 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x - 1, aktuelle_Zelle$y))
    }
    if(4 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x - 1, aktuelle_Zelle$y + 1))
    }
    if(5 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x, aktuelle_Zelle$y + 1))
    }
    if(6 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x + 1, aktuelle_Zelle$y + 1))
    }
    if(7 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x + 1, aktuelle_Zelle$y))
    }
    if(8 %in% zellen_Liste) {
      zellen_Wall <- rbind(zellen_Wall, list(aktuelle_Zelle$x + 1, aktuelle_Zelle$y - 1))
    }
    
  }
  zellen_Wall <- zellen_Wall[-1,]
  zellen_Wall <- unique(zellen_Wall)
  
  zellen_Wall <- anti_join(zellen_Wall, grabenList[[k]])
  
  for (WallZelle in 1:nrow(zellen_Wall)) {
    dem[zellen_Wall$y[WallZelle], zellen_Wall$x[WallZelle]] <- (dem[grabenList[[k]]$y[2], grabenList[[k]]$x[2]] + 0.1)
  }
  
}

writeRaster(dem,
            "S:/DGM_Flow_Modellierung/07_data_preprocessing/QGIS_Projekte/Argenschwang_KEYLINES_beiseitig_expanded.tif",
            options = c("COMPRESS=NONE"),
            overwrite = TRUE,
            datatype = "FLT8S")


### Bresenham

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
  
  

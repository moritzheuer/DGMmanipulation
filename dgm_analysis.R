library("raster")
library("sp")
library("remotes")
#install.packages("EBImageExtra")
library("whitebox")
library(dplyr)
library(tmap)
library(gridExtra)
library(terra)
library(ggplot2)
library(viridis)
library(readr)


dem <- raster("S:/DGM_Flow_Modellierung/07_data_preprocessing/Argenschwang_Hang_1m.tif", crs = '+init=EPSG:32632')

x1 <- 536
y1 <- (1855 - 225)
x2 <- 593
y2 <- (1855 - 396)

extent_raster <- extent(c(401140.5 + x1, 401140.5 + x2, 5530456 + y2, 5530456 + y1))
raster_subset_initial <- crop(dem, extent_raster)

##############
indexRaster <- raster("S:/DGM_Flow_Modellierung/07_data_preprocessing/Graben_Stützpunkte_full_vertex_index.tif", crs = '+init=EPSG:32632')

dfIndex <- data.frame(
  index = numeric(),
  column   = numeric(),
  row   = numeric()
)

for (i in 1:max(indexRaster[], na.rm = T)) {
  position_i <- which(indexRaster[] %in% i)
  
  x_i <- position_i %% indexRaster@ncols
  y_i <- (position_i - (position_i %% indexRaster@ncols)) / indexRaster@ncols + 1
  
  dfTemp <- data.frame(
    index = rep(i, length(position_i)),
    column = x_i,
    row = y_i
  )
  dfIndex <- bind_rows(dfIndex, dfTemp)
}

dfGrabenIndex <- merge(dfGrabenOrders, dfIndex, by = c("column", "row"))
dfGrabenIndex <- dfGrabenIndex[order(dfGrabenIndex$ID, dfGrabenIndex$index, decreasing = FALSE), ]



grabenList <- c()
id <- 1
for (id in unique(dfGrabenIndex$ID)) {
  grabenZellen <- nrow(dfGrabenIndex[dfGrabenIndex$ID == id,])
  partDrainage <- bresenham(x = c(dfGrabenIndex[dfGrabenIndex$ID == id,]$column[grabenZellen-1], dfGrabenIndex[dfGrabenIndex$ID == id,]$column[grabenZellen]),
                            y = c(dfGrabenIndex[dfGrabenIndex$ID == id,]$row[grabenZellen-1], dfGrabenIndex[dfGrabenIndex$ID == id,]$row[grabenZellen]))
  for (zwischenZellen in (nrow(dfGrabenIndex[dfGrabenIndex$ID == id,])-1):2) {
    partDrainage2 <- bresenham(x = c(dfGrabenIndex[dfGrabenIndex$ID == id,]$column[zwischenZellen-1], dfGrabenIndex[dfGrabenIndex$ID == id,]$column[zwischenZellen]),
                               y = c(dfGrabenIndex[dfGrabenIndex$ID == id,]$row[zwischenZellen-1], dfGrabenIndex[dfGrabenIndex$ID == id,]$row[zwischenZellen]))
    partDrainage$x <- append(partDrainage2$x, partDrainage$x)
    partDrainage$y <- append(partDrainage2$y, partDrainage$y)
  }
  grabenList[[id]] <- as.data.frame(partDrainage)
}

### Orders, und sowas, Vertices mit Koordinaten exportieren

#df <- read_delim("S:/DGM_Flow_Modellierung/07_data_preprocessing/QGIS_Projekte/Stützpunkte_mit_Koordinaten.csv", 
#                 delim = "\t", escape_double = FALSE, 
#                 col_types = cols(xcoord = col_number(), 
#                                  ycoord = col_number()), trim_ws = TRUE)

df <- read_csv("S:/DGM_Flow_Modellierung/07_data_preprocessing/QGIS_Projekte/Stützpunkte_mit_Koordinaten.csv")

df <- df[,c("id", "Ordnung", "vertex_index", "xcoord", "ycoord")]

x_min <- 401140.5
x_max <- 402511.5
y_min <- 5530455.5
y_max <- 5532310.5

b_x <- 1371
b_y <- 1855

df$xcoord <- floor(df$xcoord - x_min)

df$ycoord <- df$ycoord - y_min
df$ycoord <- (1 - (df$ycoord / b_y)) * b_y

df$ycoord <- floor(df$ycoord)

df$xcoord <- df$xcoord + 1
df$ycoord <- df$ycoord + 1

colnames(df) <- c("id", "Ordnung", "vertex_index", "column", "row")

#dem[row, col]


zähler <- 1
grabenList <- c()
for (ordnung in sort(unique(df$Ordnung))) {
  for (id in sort(unique(df$id))) {
      #grabenZellen <- max(df[df$id == id & df$Ordnung == ordnung,]$vertex_index)
    if (nrow(df[df$id == id & df$Ordnung == ordnung,]) != 0) {
      print(c(ordnung, id))
      partDrainage <- bresenham(x = c(df[df$id == id & df$Ordnung == ordnung & df$vertex_index == 1,]$column, df[df$id == id & df$Ordnung == ordnung & df$vertex_index == 0,]$column),
                                y = c(df[df$id == id & df$Ordnung == ordnung & df$vertex_index == 1,]$row, df[df$id == id & df$Ordnung == ordnung & df$vertex_index == 0,]$row))
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
  for (i in length(grabenList[[grabenNummer]]$x):1) {
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
      dem[grabenList[[grabenNummer]]$y[i-1], grabenList[[grabenNummer]]$x[i-1]] <- (dem[grabenList[[grabenNummer]]$y[i], grabenList[[grabenNummer]]$x[i]] - 0.001)
    }
    if (min(h8) < dem[grabenList[[grabenNummer]]$y[i], grabenList[[grabenNummer]]$x[i]]) {
      dem[grabenList[[grabenNummer]]$y[i-1], grabenList[[grabenNummer]]$x[i-1]] <- (min(h8) - 0.001)
    }
    
    print(dem[grabenList[[grabenNummer]]$y[i], grabenList[[grabenNummer]]$x[i]])
  }
}
  


### Ausschnitt generieren


#plot(raster_subset_initial)
#plot(raster_subset, col = pal(7))

r_points_initial = rasterToPoints(raster_subset_initial)
r_df_initial = data.frame(r_points_initial)

p1 <- ggplot(data=r_df_initial) + 
  geom_raster(aes(x=x,y=y,fill=Argenschwang_Hang_1m)) + 
  coord_equal() +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  xlab("Longitude") + ylab("Latitude") +
  scale_fill_viridis(option = "magma")


raster_subset <- crop(dem, extent_raster)
r_points = rasterToPoints(dem)
r_df = data.frame(r_points)


p2 <- ggplot(data=r_df) + 
  geom_raster(aes(x=x,y=y,fill=Argenschwang_Hang_1m)) + 
  coord_equal() +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  xlab("Longitude") + ylab("Latitude") +
  scale_fill_viridis(option = "magma")

grid.arrange(p1, p2, ncol = 2)

ggsave(
  filename = "S:/DGM_Flow_Modellierung/12_Grafiken/plots.pdf", 
  plot = marrangeGrob(list(p1, p2), nrow=1, ncol=2), 
  width = 20, height = 30
)







###############################

writeRaster(dem,'S:/DGM_Flow_Modellierung/07_data_preprocessing/Argenschwang_altered_interpol.tif',options=c('TFW=YES'), overwrite = TRUE)

wbt_fill_single_cell_pits(
  dem = "S:/DGM_Flow_Modellierung/07_data_preprocessing/Argenschwang_altered_interpol.tif", 
  output = 'S:/DGM_Flow_Modellierung/07_data_preprocessing/Argenschwang_Hang_pitFilled.tif'
)

wbt_breach_depressions_least_cost(
  dem = 'S:/DGM_Flow_Modellierung/07_data_preprocessing/Argenschwang_Hang_pitFilled.tif', 
  output = 'S:/DGM_Flow_Modellierung/07_data_preprocessing/Argenschwang_Hang_pitFilled_breached.tif', 
  dist = 5
)

wbt_d8_flow_accumulation(
  input = 'S:/DGM_Flow_Modellierung/07_data_preprocessing/Argenschwang_Hang_pitFilled_breached.tif', 
  output = 'S:/DGM_Flow_Modellierung/07_data_preprocessing/Argenschwang_Hang_d8_fac_manipulated.tif', 
  out_type = 'cells'
)

rast('S:/DGM_Flow_Modellierung/07_data_preprocessing/Argenschwang_Hang_d8_fac_manipulated.tif') %>% 
  {. ->> my_rast_d8fa}

tm_shape(my_rast_d8fa %>% log10, raster.downsample = F) +
  tm_raster(palette = 'viridis')


#####################################################


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

for (grabenNummer in 1:length(grabenList)) {
  for (i in length(grabenList[[grabenNummer]]$x):1) {
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
      dem[grabenList[[grabenNummer]]$y[i-1], grabenList[[grabenNummer]]$x[i-1]] <- (dem[grabenList[[grabenNummer]]$y[i], grabenList[[grabenNummer]]$x[i]] - 0.001)
    }
    if (min(h8) < dem[grabenList[[grabenNummer]]$y[i], grabenList[[grabenNummer]]$x[i]]) {
      dem[grabenList[[grabenNummer]]$y[i-1], grabenList[[grabenNummer]]$x[i-1]] <- (min(h8) - 0.001)
    }
    
    print(dem[grabenList[[grabenNummer]]$y[i], grabenList[[grabenNummer]]$x[i]])
  }
}

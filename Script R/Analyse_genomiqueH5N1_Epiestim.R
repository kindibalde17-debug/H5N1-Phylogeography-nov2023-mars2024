install.packages("remotes")
install.packages("devtools"); library(devtools)
remotes::install_github("sdellicour/seraphim/windows", upgrade = "never")
install.packages("diagram")
install.packages("viridisLite")

library(seraphim)
library(diagram)
library(sf)
library(viridisLite)

ls("package:seraphim")

mcc_tre = readAnnotatedNexus("Binome8_treeanot.tree")
mostRecentSamplingDatum = 2024.1693989071039
mcc_tab = mccTreeExtractions(mcc_tre, mostRecentSamplingDatum)

write.csv(mcc_tab, "Binome_8_MCC2_extracted.csv",
          row.names = FALSE, quote = FALSE)

mcc_tab = read.csv("Binome_8_MCC2_extracted.csv", header = TRUE)
localTreesDirectory = "Extracted_trees_H5N1"
allTrees = scan(file = "Binome_8b.trees",
                what = "", sep = "\n", quiet = TRUE)

burnIn = 0
randomSampling = FALSE
nberOfTreesToSample = 1000
coordinateAttributeName = "location"
treeExtractions(localTreesDirectory,
                allTrees,
                burnIn,
                randomSampling,
                nberOfTreesToSample,
                mostRecentSamplingDatum,
                coordinateAttributeName)

nberOfExtractionFiles = nberOfTreesToSample
prob = 0.80
startDatum = 2023.7164 #age du noeud le plus ancien via tracer/Figtree 
precision = 1/12   # 1 mois
polygons = suppressWarnings(spreadGraphic2(localTreesDirectory,nberOfExtractionFiles,prob,startDatum,precision))
#echelle de couleur 
colour_scale = viridis(101, option = "plasma", direction = 1)

minYear = startDatum
maxYear = mostRecentSamplingDatum

endYears_indices = round(((mcc_tab[,"endYear"] - minYear) /
                            (maxYear - minYear)) * 100) + 1
endYears_indices  = pmin(pmax(endYears_indices, 1), 101)
endYears_colours  = colour_scale[endYears_indices]

polygons_colours = rep(NA, length(polygons))
for (i in seq_along(polygons)) {
  date          = as.numeric(names(polygons[[i]]))
  polygon_index = round(((date - minYear) / (maxYear - minYear)) * 100) + 1
  polygon_index = pmin(pmax(polygon_index, 1), 101)
  polygons_colours[i] = gsub("FF","30",colour_scale[polygon_index])
}
# ============================================================
#  RASTER & SHAPEFILE
# ============================================================
borders         = shapefile("International_borders_shapefile/International_borders_shapefile/Only_international_borders.shp")
template_raster = raster("C:/Users/kindi/OneDrive/Documents/Raster_template_Europe_area (2).tif")

borders_sf = st_as_sf(borders)
borders_sf = st_simplify(borders_sf, dTolerance = 0.01, preserveTopology = TRUE)
borders    = as(borders_sf, "Spatial")
lon_range = range(c(mcc_tab$startLon, mcc_tab$endLon))
lat_range = range(c(mcc_tab$startLat, mcc_tab$endLat))
margin = 6

zoom_extent          = extent(lon_range[1] - margin, lon_range[2] + margin,
                              lat_range[1] - margin, lat_range[2] + margin)
template_raster_zoom = crop(template_raster, zoom_extent)
borders_zoom         = crop(borders, zoom_extent)

# ============================================================
#  PLOT
# ============================================================
dev.new(width = 6, height = 6.3)
par(mar=c(0,0,0,0), oma=c(3,3.5,1,1), mgp=c(0,0.4,0),
    lwd=0.2, bty="o")

# --- Fond carte ---
plot(template_raster_zoom,
     col    = "#F5F5F0",
     box    = FALSE,
     axes   = FALSE,
     colNA  = "#DCDCDC",
     legend = FALSE)

# --- Polygones HPD ---
for (i in rev(seq_along(polygons))) {
  tryCatch(
    plot(polygons[[i]], axes=FALSE, col=polygons_colours[i], add=TRUE, border=NA),
    error = function(e) message("Polygone ", i, " ignoré : ", e$message)
  )
}

# --- Frontières ---
plot(borders_zoom, add=TRUE, lwd=0.3, border="#555555")

# --- Flèches ---
for (i in 1:nrow(mcc_tab)) {
  curvedarrow(
    cbind(mcc_tab[i,"startLon"], mcc_tab[i,"startLat"]),
    cbind(mcc_tab[i,"endLon"],   mcc_tab[i,"endLat"]),
    arr.length=0.08, arr.width=0.05, lwd=0.4, lty=1,
    lcol="#333333", arr.col="#333333", arr.pos=0.95,
    curve=0.12, dr=NA, endhead=TRUE
  )
}

# --- Noeuds ---
for (i in nrow(mcc_tab):1) {
  if (i == 1) {
    points(mcc_tab[i,"startLon"], mcc_tab[i,"startLat"],
           pch=16, col=colour_scale[1], cex=0.8)
    points(mcc_tab[i,"startLon"], mcc_tab[i,"startLat"],
           pch=1, col="gray10", cex=0.8)
  }
  points(mcc_tab[i,"endLon"], mcc_tab[i,"endLat"],
         pch=16, col=endYears_colours[i], cex=0.8)
  points(mcc_tab[i,"endLon"], mcc_tab[i,"endLat"],
         pch=1, col="gray10", cex=0.8)
}

# --- Cadre & axes ---
rect(xmin(template_raster_zoom), ymin(template_raster_zoom),
     xmax(template_raster_zoom), ymax(template_raster_zoom),
     xpd=TRUE, lwd=0.2)

x_at = pretty(c(xmin(template_raster_zoom), xmax(template_raster_zoom)), n=5)
x_at = x_at[x_at > xmin(template_raster_zoom) & x_at < xmax(template_raster_zoom)]

y_at = pretty(c(ymin(template_raster_zoom), ymax(template_raster_zoom)), n=5)
y_at = y_at[y_at > ymin(template_raster_zoom) & y_at < ymax(template_raster_zoom)]

axis(1, x_at, pos=ymin(template_raster_zoom), mgp=c(0,0.4,0),
     cex.axis=0.7, lwd=0, lwd.tick=0.2, padj=-0.8, tck=-0.01,
     col.axis="gray30", labels=paste0(x_at, "°E"))

axis(2, y_at, pos=xmin(template_raster_zoom), mgp=c(0,0.6,0),
     cex.axis=0.7, lwd=0, lwd.tick=0.2, padj=1, tck=-0.01,
     col.axis="gray30", labels=paste0(y_at, "°N"))

# --- Légende colorbar (même logique que l'ancien code) ---
rast    = raster(matrix(nrow=1, ncol=2))
rast[1] = min(mcc_tab[,"startYear"])
rast[2] = max(mcc_tab[,"endYear"])

plot(rast, legend.only=TRUE, add=TRUE, col=colour_scale,
     legend.width=0.5, legend.shrink=0.3,
     smallplot=c(0.60, 0.95, 0.06, 0.10),
     horizontal=TRUE,
     legend.args=list(text="", cex=0.6, line=1.8, col="gray30"),
     axis.args=list(cex.axis=0.5, lwd=0, lwd.tick=0.2, tck=-0.5,
                    col.axis="gray30", line=0, mgp=c(0, 0.3, 0),
                    at = seq(min(mcc_tab[,"startYear"]), max(mcc_tab[,"endYear"]), length.out=5),
                    labels = round(seq(min(mcc_tab[,"startYear"]), max(mcc_tab[,"endYear"]), length.out=5), 2)))


#Analyse preliminaire de Rt#


library(EpiEstim)
library(lubridate)

tab = read.csv("dates.csv")
tab$collection_date = as.Date(gsub('"', '', tab$collection_date)) #CAR "" dans mon fichier dates

# Construire l'incidence quotidienne
dates = interval(min(tab$collection_date), tab$collection_date) %/% days(1) + 1
total_days = interval(min(tab$collection_date), max(tab$collection_date)) %/% days(1) + 1

daily_cases = rep(0, total_days)
for(i in 1:length(daily_cases)){
  daily_cases[i] = sum(dates == i)
}

# Intervalle sériel pour H5N1 
mean_si = 2
std_si  = 1

# Estimation Rt
res =  estimate_R(incid = daily_cases,
                  method = "parametric_si",
                  config = make_config(list(mean_si = mean_si,
                                            std_si = std_si)))

# Convertir les jours en dates
start_date = min(tab$collection_date)
Rt_dates = start_date + res$R$t_end - 1


dev.new()
plot(Rt_dates, res$R$`Mean(R)`,
     type="l", lwd=2, col="black",
     xlab="Date", ylab="Rt",
     main="Évolution de Rt - H5N1")
lines(Rt_dates, res$R$`Quantile.0.025(R)`, col="gray", lty=2)
lines(Rt_dates, res$R$`Quantile.0.975(R)`, col="gray", lty=2)
abline(h=1, lty=2, col="red")

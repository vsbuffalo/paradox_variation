library(rgbif)
library(rnaturalearth)
library(sf)
library(tidyverse)
library(igraph)
library(alphahull)
library(maptools)

# WGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

LAND <- rnaturalearth::ne_countries(returnclass = "sf") 
CONTINENTS <- LAND %>% group_by(continent) %>% summarize() %>% st_union()

sanitize_species <- function(x) {
  # for the Leffler et al data
  lapply(strsplit(x, ' +', perl=TRUE), 
         function(y) {
           sps <- y[2]
           if (regexpr('/', sps, fixed=TRUE) != -1) {
             spss <- strsplit(sps, '/', fixed=TRUE)
             map(spss, function(z) paste(y[1], z))
           } else {
             paste(y[1], y[2])
           }
         })
}


# get_occurrences <- function(species) {
#   key <- name_suggest(q=species, rank='species')$key[1]
#   if (length(key) == 0) {
#   }
#   occ_search(key)
# }


# get this *once*
countriesSP <- rworldmap:::getMap(resolution='low')
coords2continent = function(lat, lon) {  
  # https://stackoverflow.com/questions/21708488/get-country-and-continent-from-longitude-and-latitude-point-in-r
  points <- data.frame(lon=lon, lat=lat)
  # converting points to a SpatialPoints object
  # setting CRS directly to that from rworldmap
  pointsSP = sp:::SpatialPoints(points, proj4string=sp:::CRS(sp:::proj4string(countriesSP)))  

  # use 'over' to get indices of the Polygons object containing each point 
  indices = sp:::over(pointsSP, countriesSP)

  #indices$REGION   # returns the continent (7 continent model)
  indices$continent
}


on_map <- function(occurs, database='world', add=FALSE, 
                   direct_data=FALSE,
                   pch=19, cex=0.2, col='firebrick', ...) {
  if (!direct_data) 
    d <- occurs$data 
  else
    d <- occurs
  x <- d[, c('decimalLongitude', 'decimalLatitude')]
  if (!add) {
    plot.new()
    maps::map(database)
  }
  points(x, col=col, pch=pch, cex=cex, ...)
}


infer_is_terrestrial = function(x) {  
  if (!all(c('decimalLongitude', 'decimalLatitude') %in% colnames(x))) {
    msg  <- 'some GBIF results lack decimalLongitude and decimalLatitude as columns'
    warning(msg)
    return(NULL)
  }
  d <- filter_occurs(x) 
  lon <- d$decimalLongitude
  lat <- d$decimalLatitude
  points <- st_as_sf(x = d, 
                     coords = c("decimalLongitude", "decimalLatitude"),
                     crs=st_crs(CONTINENTS))
  indices = sf:::st_intersects(points, CONTINENTS, sparse=FALSE)
  n <- nrow(points)
  if (sum(indices) > n/2)
    return(TRUE)
  return(FALSE)
}

on_map <- function(occurs, database='world', add=FALSE, 
                   direct_data=FALSE,
                   pch=19, cex=0.2, col='firebrick', ...) {
  if (!direct_data) 
    d <- occurs$data 
  else
    d <- occurs
  x <- d[, c('decimalLongitude', 'decimalLatitude')]
  if (!add) {
    plot.new()
    maps::map(database)
  }
  points(x, col=col, pch=pch, cex=cex, ...)
}

ashape2lines <- function(x) {
  # based on: https://casoilresource.lawr.ucdavis.edu/software/r-advanced-statistical-package/working-spatial-data/converting-alpha-shapes-sp-objects/
	if(class(x) != 'ashape')
		stop('this function only works with `ashape` class objects')
	
	xdf <- as.data.frame(x$edges)
	
	# convert each edge to a line segment
	l <- lapply(seq_len(nrow(xdf)), function(i) {
		# extract line start and end points as 1x2 matrices
		p1 <- cbind(xdf$x1[i], xdf$y1[i])
		p2 <- cbind(xdf$x2[i], xdf$y2[i])
		# row-bind into 2x3 matrix
		Line(rbind(p1, p2))
  })
		
	# promote to Lines class, then to SpatialLines class
	lns <- Lines(l, ID=1)
	
	# copy over CRS data from original point data
	lspl <- tryCatch(SpatialLines(list(lns), proj4string=CRS(as.character(NA))),
                   error=function(c) {
                     warnings("could not turn ashape lines into spatial lines")
                     return(NULL)
                   })
  if (is.null(lspl)) return(NULL)
	
	# promote to SpatialLinesDataFrame, required for export to GRASS / OGR
	lspldf <- SpatialLinesDataFrame(lspl, data=data.frame(id=1), match.ID=FALSE)
	return(lspldf)
}

filter_occurs <- function(d) {
  x <- d[, c('decimalLongitude', 'decimalLatitude')]

  # filter out strange occurrences
  keep <- !duplicated(x) & apply(!is.na(x), 1, all)
  return(d[keep, ])
}

calc_area_regions <- function(occurs, ...) {
  # deprecated
  d <- filter_occurs(occurs$data)
  # infer contintinent
  d$region <- coords2continent(d$decimalLatitude, d$decimalLongitude)
  d <- d %>% filter(!is.na(region))
  d_by_region <- split(d, d$region)
  lapply(d_by_region, function(x) {
           message(sprintf("calculating α shape area for region: %s", unique(x$region)))
           calc_area(x, ...)
  })
}

calc_area_total <- function(occurs, ...) {
  d <- occurs$data 
  calc_area(d, ...)
}

occurs2coords <- function(x) {
  if (!all(c('decimalLongitude', 'decimalLatitude') %in% colnames(x))) {
    msg  <- 'some GBIF results lack decimalLongitude and decimalLatitude as columns'
    warning(msg)
    return(NULL)
  }
  return(x[, c('decimalLongitude', 'decimalLatitude')])
}

calc_area <- function(d, alpha, constrain=FALSE, 
                      clean_coords=FALSE,
                      add_outline=FALSE, is_terrestrial=NULL) {
  d <- filter_occurs(d)
  if (clean_coords)
    d <- clean_coordinates(d, lon = "decimalLongitude", lat = "decimalLatitude")
  coords <- occurs2coords(d)
  if (nrow(coords) == 0) return(list(area = NA, polygon=NULL))
  # simple way -- but often wrong
  #is_terrestrial <- !('waterBody' %in% names(x))
  if (is.null(is_terrestrial) && constrain) {
    warning('cannot constrain: is_terrestrial not set')
    constrain <- FALSE # can't constrain
  }
  
  a <- tryCatch(alphahull:::ashape(coords, alpha=alpha),
                error=function(c) {
                   if (regexpr("collinear", c$message) != -1) {
                     warning("collinear points")
                     return(NULL)
                   }
                })
   

  if (is.null(a)) return(list(area=NA, polygon=NULL))
  # convert alpha hull to spatial polygons
  #range_pg <- ah2sp(a)
  range_pg <- ashape2lines(a) 
  if (is.null(range_pg)) {
    warning("α hull did not create valid polygons")
    return(list(area=NA, polygon=NULL))
  }
  range_pg <- range_pg %>% st_as_sf() %>% 
               st_set_crs(st_crs(CONTINENTS)) %>% st_polygonize() 
  pg <- range_pg %>% st_as_sf() %>% 
           st_set_crs(st_crs(CONTINENTS)) %>%
           st_buffer(dist = 0) 

  if (constrain) {  
    if (!is_terrestrial) {
      pg <- st_geometry(st_difference(pg, st_union(CONTINENTS)))
    } else {
      pg  <- st_intersection(pg, CONTINENTS) 
    }
  }
  # in km sq
  area <- st_area(pg) %>% units:::set_units(km^2)

  # draw polygon with segments right out of alpha hull
  if (add_outline) {
    n2 <- dim(ashape)[1]
    for (i in 1:n2) {
      segments(ashape[i, "x1"], ashape[i, "y1"], ashape[i,
               "x2"], ashape[i, "y2"], col='blue')
    }
  }
  return(list(area=area, polygon=pg))
}

# Estimate range and output diagnostic plot
infer_range <- function(x, alpha=c(100, 300), col_alpha=0.4, 
                        is_terrestrial=NULL, species='',
                        constrain=TRUE) {
  if (is.null(x$data)) return(NA)
  maps::map('world', lwd=0.3, col='gray10')

  if (!is.null(is_terrestrial)) {
    if (is_terrestrial) {
      fillcol <- 'darkgreen'
      alpha <- alpha[1]
    } else {
      fillcol <- 'blue'
      alpha <- alpha[2]
    }
  } else {
    fillcol <- 'gray'
    alpha <- alpha[1]
  }

  # non-terrestrial are *not* done by region
  area <- calc_area_total(x, alpha=alpha, constrain=constrain, 
                          is_terrestrial=is_terrestrial)
  if (is.null(area$polygon)) {
    total_area <- NA
  } else {
    plot(st_geometry(area$polygon),  add=TRUE,
         border=NA, col=scales:::alpha(fillcol, col_alpha))
    total_area <- sum(unlist(area$area), na.rm=TRUE)
  }
  on_map(x, add=TRUE)
  tot <- round(total_area)
  title(main=sprintf("%s, %d km², alpha = %d, log10 area = %0.2f", species, 
                     as.integer(tot), as.integer(alpha), round(log10(tot), 1)), 
        cex.main=0.8)
  if (length(total_area) == 0) return(NA)
  return(total_area)
}


is_genus <- function(x, genus) {
  regexpr(genus, x) != -1
}

valid_occurs <- function(x) {
  coord_cols <- c('decimalLatitude', 'decimalLongitude')
  any(sapply(x, function(y) all(coord_cols %in% colnames(y$data))))
}

consensus  <- function(x) {
  names(sort(table(unlist(x)), decreasing=TRUE))[1]
}

get_terrestriality <- function(x) {
  as.logical(consensus(sapply(x, function(y) {
      infer_is_terrestrial(y$data)
  })))
}

count_occs <- function(x) {
  sum(sapply(x, function(y) {
               if (is.null(y) || is.null(y$data)) return(0L)
               nrow(y$data)
   }))
}



# 00 - Extract the five IUCN-mapped Somniosus ranges to a local GeoPackage
#
# IUCN range polygons may not be redistributed, so they are NOT in this
# repository (data/iucn/ is git-ignored). Run this script once after obtaining
# the polygons yourself; 02_figure.R then reads data/iucn/somniosus_ranges.gpkg.
#
# Two sources are accepted, chosen by environment variable:
#
#   IUCN_SHP=/path/to/SHARKS_RAYS_CHIMAERAS.shp   (or .gpkg)
#       The IUCN Red List spatial download for Sharks, Rays & Chimaeras
#       (https://www.iucnredlist.org/resources/spatial-data-download), v2024-2.
#
#   SPECIES_RANGES_DIR=/path/to/chondrichthyes.species.ranges/data
#       A local clone of the lab's chondrichthyes.species.ranges, which stores
#       the same IUCN v2024-2 polygons as one .RData sf object per species.
#       This is the default, at ~/Documents/GitHub/chondrichthyes.species.ranges/data.
#
# Filter (as in the lab's mpa.chondrichthyes): extant (presence 1), native
# (origin 1), dissolved to one feature per species.
#
# Somniosus has six valid species; IUCN maps five. S. cheni Hsu, Lin & Joung
# 2022 (Taiwan) has no IUCN range polygon and is therefore not drawn.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(sf)
})

SPECIES <- c("Somniosus microcephalus", "Somniosus rostratus", "Somniosus pacificus",
             "Somniosus antarcticus", "Somniosus longus")
OUT <- here("data", "iucn", "somniosus_ranges.gpkg")
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)

shp <- Sys.getenv("IUCN_SHP", "")
if (nzchar(shp)) {
  raw <- st_read(shp, quiet = TRUE) |> filter(sci_name %in% SPECIES)
} else {
  dir <- Sys.getenv("SPECIES_RANGES_DIR",
                    path.expand("~/Documents/GitHub/chondrichthyes.species.ranges/data"))
  raw <- map(SPECIES, function(sp) {
    f <- file.path(dir, paste0(gsub(" ", "_", sp), ".RData"))
    if (!file.exists(f)) stop("No IUCN range file for ", sp, ": ", f,
                              "\nSet IUCN_SHP or SPECIES_RANGES_DIR.")
    e <- new.env(); obj <- load(f, envir = e); e[[obj]]
  }) |> bind_rows()
}

sf_use_s2(FALSE)
ranges <- raw |>
  filter(presence == 1, origin == 1) |>
  st_make_valid() |>
  group_by(sci_name) |>
  summarise(.groups = "drop")

stopifnot(setequal(ranges$sci_name, SPECIES))
st_write(ranges, OUT, delete_dsn = TRUE, quiet = TRUE)
message("wrote ", OUT, " (", nrow(ranges), " species)")

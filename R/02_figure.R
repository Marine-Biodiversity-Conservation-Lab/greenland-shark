# 02 - Three-panel figure: sleeper-shark distribution, metabolism, r_max
#
#   A  World map of the five IUCN-mapped Somniosus ranges, plus published
#      Greenland shark records that fall outside the IUCN polygon.
#   B  Resting metabolic power corrected to a 1 kg body with the fitted mass
#      exponent, against MEASUREMENT temperature, 26 species, with the model's
#      temperature term and 95% credible band (01_fit_rmr_model.R).
#   C  Maximum intrinsic rate of population increase (r_max) against MAXIMUM
#      body mass, 72 species, coloured by ENVIRONMENTAL temperature. No model.
#
# Panels B and C carry two different thermal quantities and must not share a
# legend: B is the respirometer water (Greenland shark 4.9 C), C is the annual
# mean temperature over each species' range at its median depth (Greenland
# shark 1.72 C). That is why the Greenland shark marker is mid blue in B and
# dark blue in C.
#
# Inputs
#   data/iucn/somniosus_ranges.gpkg   from 00_extract_iucn_ranges.R (not in repo)
#   data/gls_point_records.csv
#   data/rmr_interspecific_26spp.csv
#   outputs/m2_posterior_draws.csv    from 01_fit_rmr_model.R (committed)
#   data/rmax_72spp.csv
#   data/phylopic_somniosus_microcephalus.rds
#
# Writes outputs/figures/figure.{pdf,png,svg}. The SVG keeps text live for
# editing in Illustrator.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(sf)
  library(rnaturalearth)   # also needs rnaturalearthdata for scale = 50
  library(ggrepel)
  library(patchwork)
  library(rphylopic)
})

FIG_DIR <- here("outputs", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

#### Shared constants ####

K_BOLTZ <- 8.617333e-5
GLS   <- "Somniosus microcephalus"
ENDO  <- "Isurus oxyrinchus"       # shortfin mako, regional endotherm
WHALE <- "Rhincodon typus"

GLS_BLUE <- "#2a78d6"   # GLS range colour on the map; every GLS label uses it

# Warm-to-cold temperature ramp (checked for colour-vision deficiency)
ISO_COLS <- c(`25` = "#c1272d", `15` = "#eda100", `5` = "#35a3d9", `0` = "#2a4d9b")
GLS_FILL_B <- ISO_COLS[["5"]]   # measurement temperature 4.9 C
GLS_FILL_C <- ISO_COLS[["0"]]   # environmental temperature 1.72 C

INK <- "#0b0b0b"; INK_SECOND <- "#52514e"; MIDGREY <- "#808080"; GRID <- "#e6e5e2"
FONT <- "Arial"
SZ_TICK <- 9; SZ_TITLE <- 12; SZ_TAG <- 12; SZ_LEGEND <- 9
LAB <- 9 / .pt                  # 9 pt text, in mm for geom_text

theme_pub <- function() {
  theme_classic(base_size = SZ_TICK, base_family = FONT) +
    theme(axis.text    = element_text(size = SZ_TICK, colour = "grey20"),
          axis.title   = element_text(size = SZ_TITLE, colour = "grey20"),
          axis.line    = element_blank(),
          axis.ticks   = element_line(linewidth = 0.3, colour = "grey60"),
          panel.border = element_rect(fill = NA, colour = "grey60", linewidth = 0.4),
          panel.grid   = element_blank(),
          plot.margin  = margin(4, 12, 4, 4))
}
TAG <- theme(plot.tag = element_text(family = FONT, face = "bold", size = SZ_TAG,
                                     colour = INK))

# PhyloPic silhouette of Somniosus microcephalus, Margot Michaud 2018, CC0 1.0
# (uuid 50f67ea6-a922-4ede-a054-af9632cff470), cached so no API call is needed.
gls_img <- readRDS(here("data", "phylopic_somniosus_microcephalus.rds"))


#### Panel A - Somniosus ranges ####

RANGES_GPKG <- here("data", "iucn", "somniosus_ranges.gpkg")
if (!file.exists(RANGES_GPKG))
  stop("IUCN ranges not found. Obtain the IUCN polygons and run R/00_extract_iucn_ranges.R first.")

SPECIES <- tribble(
  ~sci_name,                  ~colour,
  "Somniosus microcephalus",  GLS_BLUE,
  "Somniosus rostratus",      "#eb6834",
  "Somniosus pacificus",      "#1baf7a",
  "Somniosus antarcticus",    "#eda100",
  "Somniosus longus",         "#e87ba4"
)
pal <- setNames(SPECIES$colour, SPECIES$sci_name)
MAP_INK <- "#2b2b2b"
short_name <- function(x) sub("Somniosus ", "S. ", x)

sf_use_s2(FALSE)
ranges <- st_read(RANGES_GPKG, quiet = TRUE) |>
  rename(geometry = any_of(c("geom", "geometry"))) |>
  mutate(sci_name = factor(sci_name, levels = SPECIES$sci_name)) |>
  arrange(sci_name)
stopifnot(nrow(ranges) == 5)

land <- ne_countries(scale = 50, returnclass = "sf") |> st_make_valid() |> summarise()

gls_range <- filter(ranges, sci_name == GLS)
LON0 <- round(mean(st_bbox(gls_range)[c("xmin", "xmax")]))   # map centre longitude

# Point records: drawn only if OUTSIDE the IUCN polygon (tested on the ellipsoid)
gls_records <- read.csv(here("data", "gls_point_records.csv")) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
local({
  old <- sf_use_s2(TRUE); on.exit(suppressMessages(sf_use_s2(old)))
  poly <- st_make_valid(st_union(st_transform(gls_range, 4326)))
  inside <- as.vector(st_intersects(gls_records, poly, sparse = FALSE))
  if (any(inside)) message("GLS records inside the IUCN range, not drawn: ",
                           paste(gls_records$record[inside], collapse = ", "))
  gls_records <<- gls_records[!inside, ]
})

lonlat_box <- function(xmin, xmax, ymin, ymax, step = 1) {
  xs <- seq(xmin, xmax, length.out = ceiling((xmax - xmin) / step) + 1)
  ys <- seq(ymin, ymax, length.out = ceiling((ymax - ymin) / step) + 1)
  ring <- rbind(cbind(xs, ymin), cbind(xmax, ys), cbind(rev(xs), ymax), cbind(xmin, rev(ys)))
  st_sfc(st_polygon(list(ring)), crs = 4326)
}

CRS_MAP <- sprintf("+proj=robin +lon_0=%d +datum=WGS84 +units=m", LON0)
break_map <- function(x) suppressWarnings(st_break_antimeridian(x, lon_0 = LON0)) |>
  st_transform(CRS_MAP)
ranges_m <- break_map(ranges)
land_m   <- break_map(land)
globe    <- lonlat_box(LON0 - 180 + 1e-6, LON0 + 180 - 1e-6, -90, 90) |> st_transform(CRS_MAP)

# Larger ranges underneath, GLS on top
draw_order <- function(x) x |> mutate(.a = as.numeric(st_area(x))) |>
  arrange(sci_name == GLS, desc(.a)) |> select(-.a)

# Default label anchor: a point on each species' largest polygon piece
lab <- ranges_m |>
  st_collection_extract("POLYGON", warn = FALSE) |>
  st_cast("MULTIPOLYGON", warn = FALSE) |> st_cast("POLYGON", warn = FALSE) |>
  mutate(.a = as.numeric(st_area(geometry))) |>
  group_by(sci_name) |> slice_max(.a, n = 1, with_ties = FALSE) |> ungroup() |>
  (\(d) suppressWarnings(st_point_on_surface(d)))() |>
  mutate(lab = short_name(sci_name),
         x = st_coordinates(geometry)[, 1], y = st_coordinates(geometry)[, 2]) |>
  st_drop_geometry()

# Hand-placed labels (lon/lat of label; hj/vj its justification)
LAB_POS <- tribble(
  ~lab,               ~lon, ~lat, ~hj, ~vj,
  "S. microcephalus",  -50,    3, 0.5, 0.5,
  "S. rostratus",       -3,   12, 0,   1,
  "S. longus",          95,   38, 0.5, 0.5
)
pos <- st_coordinates(st_transform(st_as_sf(LAB_POS, coords = c("lon", "lat"), crs = 4326), CRS_MAP))
lab <- lab |>
  left_join(tibble(lab = LAB_POS$lab, lx = pos[, 1], ly = pos[, 2],
                   hj = LAB_POS$hj, vj = LAB_POS$vj), by = "lab") |>
  mutate(sx = NA_real_, sy = NA_real_)

# Vertical leader from (sx, sy) up to where it first enters the range
range_of <- function(sp) st_union(ranges_m[ranges_m$sci_name == sp, ])
vertical_anchor <- function(lab, label, sp, sx, sy, inset = 1.5e5) {
  i <- which(lab$lab == label)
  up  <- st_sfc(st_linestring(rbind(c(sx, sy), c(sx, sy + 2e7))), crs = st_crs(ranges_m))
  hit <- st_coordinates(suppressWarnings(st_intersection(up, range_of(sp))))
  if (nrow(hit) == 0) stop("vertical leader for ", label, " misses its range")
  lab$sx[i] <- sx; lab$sy[i] <- sy; lab$x[i] <- sx; lab$y[i] <- min(hit[, "Y"]) + inset
  lab
}
i <- which(lab$lab == "S. microcephalus")
lab <- vertical_anchor(lab, "S. microcephalus", GLS, lab$lx[i], lab$ly[i])
i <- which(lab$lab == "S. rostratus")
lab <- vertical_anchor(lab, "S. rostratus", "Somniosus rostratus", lab$lx[i] + 4e5, lab$ly[i])

# S. pacificus: label directly below the southernmost point of its range
pac_pts <- st_coordinates(range_of("Somniosus pacificus"))
pac_s   <- pac_pts[which.min(pac_pts[, "Y"]), c("X", "Y")]
i <- which(lab$lab == "S. pacificus")
lab[i, c("lx", "ly", "sx", "sy", "x", "y", "hj", "vj")] <-
  list(pac_s[1], pac_s[2] - 1.1e6, pac_s[1], pac_s[2] - 1.1e6, pac_s[1], pac_s[2], 0.5, 0.5)

lab$col <- ifelse(lab$lab == "S. microcephalus", GLS_BLUE, MAP_INK)
lab$sx  <- ifelse(is.na(lab$sx), lab$lx, lab$sx)
lab$sy  <- ifelse(is.na(lab$sy), lab$ly, lab$sy)
rep_lab <- lab[is.na(lab$lx), ]    # S. antarcticus: placed by ggrepel
fix_lab <- lab[!is.na(lab$lx), ]

pA <- ggplot() +
  geom_sf(data = globe, fill = "white", colour = "#9a9992", linewidth = 0.3) +
  geom_sf(data = land_m, fill = "#e4e3df", colour = NA) +
  geom_sf(data = draw_order(ranges_m), aes(fill = sci_name, colour = sci_name),
          alpha = 0.38, linewidth = 0.35) +
  geom_sf(data = st_transform(gls_records, CRS_MAP), inherit.aes = FALSE,
          shape = 21, size = 1.5, stroke = 0.35, fill = GLS_BLUE, colour = MAP_INK) +
  geom_label_repel(data = rep_lab, aes(x, y, label = lab), fontface = "italic",
                   size = 2.6, colour = rep_lab$col, fill = alpha("white", 0.8),
                   label.size = NA, label.padding = unit(0.6, "mm"),
                   segment.colour = MAP_INK, segment.size = 0.25,
                   min.segment.length = 0, box.padding = 0.4, seed = 15) +
  geom_segment(data = fix_lab, aes(x = sx, y = sy, xend = x, yend = y),
               colour = MAP_INK, linewidth = 0.25) +
  geom_label(data = fix_lab, aes(lx, ly, label = lab, hjust = hj, vjust = vj),
             fontface = "italic", size = 2.6, colour = fix_lab$col,
             fill = alpha("white", 0.8), linewidth = 0.25,
             label.padding = unit(0.6, "mm")) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_colour_manual(values = pal, guide = "none") +
  coord_sf(crs = CRS_MAP, expand = FALSE, datum = NA) +
  theme_void(base_size = 8) +
  theme(plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA))


#### Panel B - mass-corrected metabolic power vs measurement temperature ####

rmr <- read.csv(here("data", "rmr_interspecific_26spp.csv")) |>
  mutate(ln_watts = log(metabolic_power_w), ln_mass = log(mass_kg))
T_REF <- median(rmr$temp_c + 273.15)
x_T_at <- function(temp_c) 1 / (K_BOLTZ * T_REF) - 1 / (K_BOLTZ * (temp_c + 273.15))

draws <- read.csv(here("outputs", "m2_posterior_draws.csv"))
E_hat <- mean(draws$b_x_T)        # posterior means, as brms::fixef reports
b_m   <- mean(draws$b_ln_mass)

# Each species' power corrected to a 1 kg body: W * (1 / M)^b
rmr <- rmr |>
  mutate(mc_watts = exp(ln_watts - b_m * ln_mass),
         tband = cut(temp_c, c(-Inf, 10, 20, Inf), labels = c("5", "15", "25")))

t_seq <- seq(min(rmr$temp_c), max(rmr$temp_c), length.out = 100)
mc_line <- map_dfr(t_seq, function(tc) {
  v <- draws$b_Intercept + draws$b_x_T * x_T_at(tc)
  tibble(temp_c = tc, fit = exp(median(v)),
         lo = exp(quantile(v, 0.025)), hi = exp(quantile(v, 0.975)))
})

gls_b  <- rmr |> filter(species == GLS)  |> transmute(x = temp_c, y = mc_watts)
endo_b <- rmr |> filter(species == ENDO) |> transmute(x = temp_c, y = mc_watts)

pB <- ggplot() +
  geom_ribbon(data = mc_line, aes(temp_c, ymin = lo, ymax = hi), fill = GRID, colour = NA) +
  geom_line(data = mc_line, aes(temp_c, fit), linewidth = 0.6, colour = INK_SECOND) +
  geom_segment(data = gls_b, aes(x = x, y = 0.031, xend = x, yend = y),
               linewidth = 0.3, colour = GLS_BLUE) +
  geom_point(data = rmr, aes(temp_c, mc_watts, fill = tband),
             shape = 21, size = 2.1, stroke = 0.5, colour = "white") +
  geom_point(data = gls_b, aes(x, y), shape = 21, size = 2.8, stroke = 0.7,
             fill = GLS_FILL_B, colour = INK) +
  annotate("text", x = gls_b$x, y = 0.026, label = "Greenland shark",
           hjust = 0, size = LAB, family = FONT, colour = GLS_BLUE) +
  add_phylopic(img = gls_img, x = gls_b$x + 3.2, y = 0.0185, height = 0.22,
               fill = GLS_BLUE, alpha = 1) +
  geom_segment(data = endo_b, aes(x = x, y = 1.05, xend = x, yend = y * 1.12),
               linewidth = 0.3, colour = MIDGREY) +
  annotate("text", x = endo_b$x, y = 1.25, label = "Shortfin mako shark",
           hjust = 0.5, size = LAB, family = FONT, colour = MIDGREY) +
  annotate("text", x = 30.5, y = 0.017, hjust = 1, vjust = 0, size = LAB, family = FONT,
           colour = INK_SECOND,
           label = sprintf("E = %.2f eV, mass exponent %.2f", E_hat, b_m)) +
  scale_fill_manual(values = ISO_COLS[c("5", "15", "25")], guide = "none") +
  scale_x_continuous(limits = c(3.5, 31), breaks = seq(5, 30, 5)) +
  scale_y_log10(labels = scales::label_number(drop0trailing = TRUE), limits = c(0.014, 1.7)) +
  annotation_logticks(sides = "l", linewidth = 0.2, short = unit(0.4, "mm"),
                      mid = unit(0.7, "mm"), long = unit(1, "mm"), colour = INK_SECOND) +
  labs(x = "Measurement temperature (°C)", y = "Metabolic power at 1 kg (W)") +
  theme_pub()


#### Panel C - r_max vs maximum body mass, 72 species ####

rmax <- read.csv(here("data", "rmax_72spp.csv")) |>
  mutate(tband = cut(temp_env_c, c(-Inf, 5, 10, 20, Inf), labels = c("0", "5", "15", "25")))
gls_c   <- rmax |> filter(species == GLS)
whale_c <- rmax |> filter(species == WHALE)

# y-axis label: quoted exponent keeps a hyphen-minus that Arial has
YLAB_RMAX <- expression(paste(italic(r)[plain(max)], " (year"^"-1", ")"))

pC <- ggplot() +
  geom_segment(data = gls_c, aes(x = 1800, y = rmax, xend = max_mass_kg, yend = rmax),
               linewidth = 0.3, colour = GLS_BLUE) +
  geom_segment(data = whale_c, aes(x = 13000, y = rmax, xend = max_mass_kg, yend = rmax),
               linewidth = 0.3, colour = MIDGREY) +
  geom_point(data = rmax, aes(max_mass_kg, rmax, fill = tband),
             shape = 21, size = 2.1, stroke = 0.5, colour = "white") +
  geom_point(data = gls_c, aes(max_mass_kg, rmax), shape = 21, size = 2.8, stroke = 0.7,
             fill = GLS_FILL_C, colour = INK) +
  annotate("text", x = 1650, y = gls_c$rmax, label = "Greenland shark",
           hjust = 1, size = LAB, family = FONT, colour = GLS_BLUE) +
  annotate("text", x = 12000, y = whale_c$rmax, label = "Whale shark",
           hjust = 1, size = LAB, family = FONT, colour = MIDGREY) +
  scale_fill_manual(values = ISO_COLS[c("0", "5", "15", "25")],
                    labels = c("< 5", "5–10", "10–20", "> 20"),
                    name = "Environmental\ntemperature (°C)",
                    guide = guide_legend(reverse = TRUE)) +
  scale_x_log10(limits = range(rmax$max_mass_kg),
                labels = scales::label_number(drop0trailing = TRUE),
                expand = expansion(mult = c(0.06, 0.08))) +
  scale_y_log10(labels = scales::label_number(drop0trailing = TRUE),
                expand = expansion(mult = c(0.30, 0.06))) +
  annotation_logticks(sides = "bl", linewidth = 0.2, short = unit(0.4, "mm"),
                      mid = unit(0.7, "mm"), long = unit(1, "mm"), colour = INK_SECOND) +
  labs(x = "Maximum body mass (kg)", y = YLAB_RMAX) +
  theme_pub() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.02, 0.02),
        legend.justification = c(0, 0),
        legend.title = element_text(size = SZ_LEGEND, colour = "grey20"),
        legend.text  = element_text(size = SZ_LEGEND, colour = "grey20"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.size = unit(3.2, "mm"))


#### Assemble and export ####

fig <- (pA + labs(tag = "A") + TAG) /
  (pB + labs(tag = "B") + TAG) /
  (pC + labs(tag = "C") + TAG) +
  plot_layout(heights = c(0.72, 1, 1))

out <- file.path(FIG_DIR, "figure")
ggsave(paste0(out, ".pdf"), fig, width = 120, height = 230, units = "mm", device = cairo_pdf)
ggsave(paste0(out, ".png"), fig, width = 120, height = 230, units = "mm", dpi = 400,
       device = ragg::agg_png)
ggsave(paste0(out, ".svg"), fig, width = 120, height = 230, units = "mm",
       device = svglite::svglite, fix_text_size = FALSE)
message("wrote ", out, ".{pdf,png,svg}")

# Greenland shark: range, metabolism and population growth

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24208838.svg)](https://doi.org/10.5281/zenodo.24208838)

This repository holds the minimal data and code needed to reproduce the
three-panel figure on the Greenland shark (*Somniosus microcephalus*) in:

Dulvy, N. K. & VanderWright, W. J. (2027). Quick Guide Greenland Shark. *Current Biology,* 37**,** R1–R3. 

The figure puts the Greenland shark in context of (A) its genus's
distribution, (B) the temperature dependence of shark resting metabolic rate,
and (C) the maximum intrinsic rate of population increase across sharks.

## Quick start

```sh
# 1. IUCN range maps are not distributed here (see "IUCN range maps" below).
#    Point the extractor at your own copy of the IUCN polygons:
IUCN_SHP=/path/to/SHARKS_RAYS_CHIMAERAS.shp Rscript R/00_extract_iucn_ranges.R

# 2. Refit the metabolic-rate model (about 30 s; optional, because its
#    posterior draws are committed in outputs/)
Rscript R/01_fit_rmr_model.R

# 3. Draw the figure -> outputs/figures/figure.{pdf,png,svg}
Rscript R/02_figure.R
```

Run from the repository root (or open `greenland-shark.Rproj`);
paths are resolved with `here`.

## Repository layout

| Path | Contents |
|---|---|
| `R/00_extract_iucn_ranges.R` | Extracts the five mapped *Somniosus* ranges from the IUCN polygons to `data/iucn/somniosus_ranges.gpkg` (git-ignored). |
| `R/01_fit_rmr_model.R` | Fits the phylogenetic metabolic-rate model (panel B); writes the posterior draws and a summary. |
| `R/02_figure.R` | Builds all three panels and writes the figure. |
| `data/rmr_interspecific_26spp.csv` | Species-level resting metabolic rate, body mass and measurement temperature (26 species). |
| `data/phylogeny_26spp.nex` | Phylogeny pruned to those 26 species. |
| `data/rmax_72spp.csv` | *r*<sub>max</sub>, maximum body mass and environmental temperature (72 species). |
| `data/gls_point_records.csv` | Published Greenland shark records, with coordinates and sources (16 records). |
| `data/phylopic_somniosus_microcephalus.rds` | Cached Greenland shark silhouette (PhyloPic, CC0). |
| `outputs/m2_posterior_draws.csv` | Posterior draws of the panel B model (12 000 draws). |
| `outputs/m2_summary.csv` | Posterior medians and 95% credible intervals. |

Git-ignored and regenerated locally: `data/iucn/` (IUCN polygons),
`outputs/fits/` (the `brms` fit object) and `outputs/figures/` (the rendered
figure, whose PDF and SVG embed the IUCN polygons as vector geometry).

## Analyses

### Panel A — distribution of the sleeper sharks

Range polygons for the five *Somniosus* species mapped by the IUCN Red List
(*S. microcephalus*, *S. rostratus*, *S. pacificus*, *S. antarcticus*,
*S. longus*; spatial data version 2024-2) are filtered to extant (presence = 1)
and native (origin = 1) parts and dissolved to one feature per species. The
sixth valid species, *S. cheni* Hsu, Lin & Joung 2022, has no IUCN range map
and is not shown. The map uses a Robinson projection centred on the central
longitude of the Greenland shark range, with Natural Earth land (1:50 m).

Points are published Greenland shark records **that fall outside the IUCN
polygon**;  records inside the polygon (here records 11–13) are not drawn. Records 1–7 are given with
coordinates in the species account of Mecklenburg et al. (2018). Records 9–16
were read off that account's range map after georeferencing it on its known
points (azimuthal equidistant; positions good to about ±0.3°). Record 8 is a
capture off Glover's Reef, Belize (Kasana et al. 2022) identified only to genus
(*Somniosus* sp.). Records 14–16 lie outside the polygon by less than the
reading error and should not be read as range extensions. Primary citations
for each record are in `data/gls_point_records.csv`.

### Panel B — temperature dependence of resting metabolic rate

**Data.** Resting metabolic rates for 26 shark species compiled from published
respirometry studies (Gravel et al. 2024), converted to whole-animal metabolic
power (W) using an oxycalorific coefficient of 13.6 J mg⁻¹ O₂ (Elliott &
Davison 1975), after converting rates reported per unit mass back to whole-animal
rates using each study's stated scaling exponent. Each species contributes one
value (a published mean, or a mean over individuals) at its measured body mass
and measurement temperature. The Greenland shark value is a single individual
of 40.8 kg measured at 4.9 °C.

**Model.** A phylogenetic Bayesian regression in `brms` (Bürkner 2017):

ln *W*<sub>i</sub> = β₀ + *E* · *x*<sub>T,i</sub> + *b* · ln *M*<sub>i</sub> + ε,
  ε ~ MVN(0, σ² **V**(λ))

where *x*<sub>T</sub> = 1/(*k T*<sub>ref</sub>) − 1/(*k T*<sub>i</sub>) is the
Boltzmann–Arrhenius temperature term (Gillooly et al. 2001; *k* =
8.617 × 10⁻⁵ eV K⁻¹, *T*<sub>ref</sub> = median measurement temperature, 294.69 K),
so *E* is an apparent activation energy (eV) and *b* the mass-scaling exponent.
**V**(λ) = λ**A** + (1 − λ)**I**, where **A** is the phylogenetic correlation
matrix from a tree pruned from the chondrichthyan phylogeny of Stein et al.
(2018), and λ is Pagel's (1999) λ. With one observation per species λ is not
separately identifiable from the residual variance, so it was fixed at its
profile-likelihood maximum, λ = 0.816 (support interval 0.38–0.98; λ = 0
rejected, ΔAIC = 9.2). Priors: *b* ~ N(0.8, 0.2), *E* ~ N(0.65, 0.3),
σ ~ Student-t(3, 0, 1); 4 chains × 4000 iterations (1000 warm-up).

| Term | Posterior median | 95% CrI |
|---|---|---|
| Activation energy *E* (eV) | 0.42 | 0.19 – 0.66 |
| Mass exponent *b* | 0.94 | 0.84 – 1.05 |
| Intercept (ln W at 1 kg, 21.5 °C) | −1.46 | −1.92 – −0.99 |

**Plotted.** Each species' metabolic power corrected to a 1 kg body,
*W* · *M*<sup>−*b*</sup>, using the posterior mean of *b*, against measurement
temperature. The line and band are the posterior median and 95% credible
interval of exp(β₀ + *E x*<sub>T</sub>). The model is phylogenetic, so the line
need not pass through the centre of the raw points. Point colour bins
measurement temperature (< 10, 10–20, > 20 °C). The shortfin mako (*Isurus
oxyrinchus*) is labelled because, as a regional endotherm, its measurement
temperature understates its tissue temperature.

### Panel C — maximum intrinsic rate of population increase

**Data.** For 72 shark species: *r*<sub>max</sub>, maximum body mass and
environmental temperature, taken from the lab's *r*<sub>max</sub> compilation
(`RMAX_GSH`).

- ***r*<sub>max</sub>** is estimated from age at maturity, maximum age and
  annual fecundity, with natural mortality derived from lifespan (Pardo et al.
  2016), using life-history data compiled by Finucci et al. (2024) and Gravel
  et al. (2024).
- **Maximum body mass** is *a* · *L*<sub>max</sub><sup>*b*</sup>, with the
  length–weight coefficients *a*, *b* the FishBase Bayesian estimates (Froese,
  Thorson & Reyes 2014) obtained with `rfishbase` (Boettiger, Lang &
  Wainwright 2012) and *L*<sub>max</sub> the maximum total length from the
  lab's habitat compilation. Two exceptions: *L*<sub>max</sub> for
  *S. microcephalus* (756 cm TL) and *S. pacificus* (700 cm TL) is from Ebert,
  Dando & Fowler (2021); and the whale shark (*Rhincodon typus*) is set to
  34 t, because extrapolating *a* · *L*<sup>*b*</sup> to its 20 m length gives
  an implausible 77.5 t. Columns `max_mass_source` and `lmax_source` record
  which applies.
- **Environmental temperature** is the World Ocean Atlas 2018 annual mean
  temperature (Locarnini et al. 2018) at each species' median depth, averaged
  over its core range (AquaMaps probability ≥ 0.9; Kaschner et al.). It is a
  different quantity from the measurement temperature of panel B: 1.72 °C
  against 4.9 °C for the Greenland shark, hence the different marker colours.

**Plotted.** *r*<sub>max</sub> against maximum body mass on log–log axes,
coloured by environmental temperature (< 5, 5–10, 10–20, > 20 °C). No model is
fitted in this panel.

## IUCN range maps

The IUCN Red List range polygons are **not** included, in line with the IUCN
Red List terms of use, which do not permit redistribution of the spatial data.
To reproduce panel A:

1. Download the *Sharks, Rays and Chimaeras* range data (version 2024-2) from
   <https://www.iucnredlist.org/resources/spatial-data-download>.
2. Run `IUCN_SHP=/path/to/the.shp Rscript R/00_extract_iucn_ranges.R`.
   Lab members can instead use a clone of
   `Marine-Biodiversity-Conservation-Lab/chondrichthyes.species.ranges`, which
   holds the same polygons: set `SPECIES_RANGES_DIR` to its `data/` folder (the
   default is `~/Documents/GitHub/chondrichthyes.species.ranges/data`).

A later IUCN version may draw slightly different ranges.

## Software

Run with R 4.5.2 and: `brms` 2.23.0 (`rstan` 2.32.7 backend), `posterior`
1.6.1, `ape` 5.8.1, `sf` 1.1.3, `rnaturalearth` 1.2.0 + `rnaturalearthdata`
1.0.0, `ggrepel` 0.9.8, `patchwork` 1.3.2, `rphylopic` 1.7.0, `ragg` 1.5.2,
`svglite` 2.2.2, `tidyverse` 2.0.0, `here` 1.0.2. Text is set in Arial, which
must be installed for the PDF and PNG to match.

Refitting reproduces the reported estimates to two decimals; draws are not
bit-identical to `outputs/m2_posterior_draws.csv` across machines or Stan
versions.

## Data sources to cite

| Data | Reference |
|---|---|
| Resting metabolic rate | Gravel et al. (2024), and the primary respirometry studies it compiles |
| Oxycalorific coefficient | Elliott & Davison (1975) |
| *r*<sub>max</sub> method | Pardo et al. (2016) |
| *r*<sub>max</sub> life-history inputs | Finucci et al. (2024) ; Gravel et al. (2024) |
| Length–weight coefficients | Froese, Thorson & Reyes (2014); Froese & Pauly, FishBase; Boettiger, Lang & Wainwright (2012) |
| *L*<sub>max</sub> for two *Somniosus* | Ebert, Dando & Fowler (2021) |
| Maximum length (other species) | Finucci et al. (2024) ; Gravel et al. (2024) |
| Whale shark maximum mass (34 t) | Calculated by Wade VanderWright based on FishBase L-W coefficients |
| Environmental temperature | Matsushiba et al. (submitted) |
| Phylogeny | Stein, Mull et al. (2018) |
| Range maps | IUCN (2024) |
| Greenland shark point records | Mecklenburg et al. (2018); Kasana et al. (2022); Fossen et al. (2008) |
| Land outlines | Natural Earth |
| Silhouette | Margot Michaud (2018), PhyloPic, CC0 1.0; Gearty & Jones (2023) |

## References

- Boettiger, C., Lang, D.T. & Wainwright, P.C. (2012) rfishbase: exploring, manipulating and visualizing FishBase data from R. *Journal of Fish Biology* 81: 2030–2039.
- Bürkner, P.-C. (2017) brms: an R package for Bayesian multilevel models using Stan. *Journal of Statistical Software* 80(1): 1–28.
- Ebert, D.A., Dando, M. & Fowler, S. (2021) *Sharks of the World: A Complete Guide*. Princeton University Press.
- Elliott, J.M. & Davison, W. (1975) Energy equivalents of oxygen consumption in animal energetics. *Oecologia* 19: 195–201.
- Finucci, B. et al. (2024) Fishing for oil and meat drives irreversible defaunation of deepwater sharks and rays. *Science* 383: 1135–1141. 
- Fossen, I., Cotton, C.F., Bergstad, O.A. & Dyb, J.E. (2008) Species composition and distribution patterns of fishes captured by longlines on the Mid-Atlantic Ridge. *Deep-Sea Research II* 55: 203–217. 
- Froese, R., Thorson, J.T. & Reyes, R.B. (2014) A Bayesian approach for estimating length-weight relationships in fishes. *Journal of Applied Ichthyology* 30: 78–85.
- Froese, R. & Pauly, D. (eds) FishBase. www.fishbase.org, 9th November 2026.
- Gearty, W. & Jones, L.A. (2023) rphylopic: an R package for fetching, transforming, and visualising PhyloPic silhouettes. *Methods in Ecology and Evolution* 14: 2700–2708.
- Gillooly, J.F., Brown, J.H., West, G.B., Savage, V.M. & Charnov, E.L. (2001) Effects of size and temperature on metabolic rate. *Science* 293: 2248–2251.
- Gravel, S., Bigman, J.S., Pardo, S.A., Wong, S. & Dulvy, N.K. (2024) Metabolism, population growth and the fast-slow continuum of marine fishes *Fish and Fisheries* 25: 349–361. 
- IUCN (2024) *The IUCN Red List of Threatened Species*. Version 2024-2. https://www.iucnredlist.org
- Kasana, D., Martinez, H. D., Faux, O., Monzon, N., Guerra, E. & Chapman, D. D. (2022). First report of a sleeper shark (*Somniosus* sp.) in the western Caribbean, off the insular slope of a coral atoll. *Marine Biology,* 169: 101-104.
- Locarnini, R.A. et al. (2018) *World Ocean Atlas 2018, Volume 1: Temperature*. NOAA Atlas NESDIS 81. 
- Matsushiba, J. H., Ruiz-Garcia, D., Aitchison, R. M., VanderWright, W. J., Arnold, A. E., Sherman, C. S., Haque, A. B. & Dulvy, N. K. (Submitted). ocean3d: An R package for three-dimensional spatial analyses of marine habitat. *ROpenScience*.  https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d
- Mecklenburg, C.W., Lynghammar, A., Johannesen, E., Byrkjedal, I., Christiansen, J.S., Dolgov, A.V., Karamushko, O.V., Mecklenburg, T.A., Møller, P.R., Steinke, D. & Wienerroither, R.M. (2018) *Marine Fishes of the Arctic Region*. Conservation of Arctic Flora and Fauna (CAFF), Akureyri, Iceland. 
- Pagel, M. (1999) Inferring the historical patterns of biological evolution. *Nature* 401: 877–884.
- Pardo, S.A., Kindsvater, H.K., Reynolds, J.D. & Dulvy, N.K. (2016) Maximum intrinsic rate of population increase in sharks, rays, and chimaeras: the importance of survival to maturity. *Canadian Journal of Fisheries and Aquatic Sciences* 73: 1159–1163.
- Stein, R.W., Mull, C.G., Kuhn, T.S., Aschliman, N.C., Davidson, L.N.K., Joy, J.B., Smith, G.J., Dulvy, N.K. & Mooers, A.O. (2018) Global priorities for conserving the evolutionary history of sharks, rays and chimaeras. *Nature Ecology & Evolution* 2: 288–298.

## Citation

Please cite the paper above. To cite this repository (data and code), use:

Dulvy, N. K. & VanderWright, W. J. (2026). Greenland shark: range, metabolism
and population growth (data and code). Zenodo.
https://doi.org/10.5281/zenodo.24208838

This concept DOI always resolves to the latest version; v1.0.0 is
https://doi.org/10.5281/zenodo.24208839.

## Licence

Code (`R/`) is released under the MIT licence (`LICENSE`). Data in `data/` and
`outputs/` are released under CC BY 4.0 (`LICENSE-data`). The IUCN range data
are not distributed here, are not covered by either licence, and remain
subject to the IUCN Red List terms of use. The PhyloPic silhouette is CC0 1.0.

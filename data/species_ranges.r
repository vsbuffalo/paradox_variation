library(rgbif)
library(rnaturalearth)
library(sf)
library(tidyverse)
library(readxl)
library(igraph)
library(units)

geomean <- function(x) prod(x)^(1/length(x))
source('../R/range_funcs.r')


#### Load the data in 

# we don't load in Corbett-Detig data, which has ranges though some are off.
# This is merged in combined_dataset.r. I make one correction here though, to
# Equus ferus przewalskii.

## Leffler data
message('loading Leffler et al data...')
dl <- read_tsv('./leffler_et_al_2012_updated.tsv')
lsp <- sanitize_species(dl$species)
names(lsp) <- dl$species

# get Leffler species names
uniq_lsp_names <- unique(names(lsp))
uniq_lsp <- lsp[uniq_lsp_names]

## Stapley data 
message('loading Stapley et al data...')
ds <- read_tsv('./stapley_et_al_2017_updated.tsv')
# Stapley is all unique species, double check
stopifnot(length(unique(ds$species)) == nrow(ds))

## Romiguier data 
# there are some multiple observations, but not combined like Leffler data
# (hence no need for key)
message('loading Romiguier et al data...')
drs <- read_tsv('./romiguier_et_al_2014_updated.tsv')
drs_sps <- unique(drs$species)

## Combine all data
# flatten the Leffler data
species_df <- bind_rows(map2(lsp, names(lsp), function(x, n) {
                 tibble(key=rep(n, length(x)), species=unlist(x),
                        dataset='Leffler')
}))

species_df <- bind_rows(species_df, tibble(key=ds$species, 
                        species=ds$species, dataset='Stapley'),
                        tibble(key=drs_sps, species=drs_sps,
                               dataset='Romiguier'))

uniq_species_df <- species_df %>% distinct(.keep_all=TRUE)
nrow(uniq_species_df)

#### Species keys and occurence data from GBIF
safe_occ_search <- safely(occ_search)

OCC_FILE <- 'gbif_species_occ.Rdata'
if (!file.exists(OCC_FILE)) {
  message('cached GBIF occurence data et al not found -- downloading from GBIF')
  message('getting suggested GBIF species names...')
  uniq_species_df <- uniq_species_df %>% 
                       mutate(taxon_key = map(species, 
                                        ~ name_suggest(., rank='species')))

  message('downloading GBIF species occurrences...')
  n <- nrow(uniq_species_df)
  uniq_species_df <- uniq_species_df %>% 
                       mutate(occs = map2(taxon_key, 1:nrow(uniq_species_df),
                         function(x, i) {
                           keep_keys <- x$key[!is.na(x$canonicalName)]
                           map(keep_keys, function(k) {
                             msg <-"    species '%s' %d/%d\n"
                             sps <- species_df$species[i]
                             message(sprintf(msg, sps, i, n))
                             out <- safe_occ_search(k, 
                                                    fields='minimal', 
                                                    limit = 30000)
                             return(out$result)
                         })
                        }))

  save(uniq_species_df, file=OCC_FILE)
} else {
  message('loading cached GBIF occurence data...   ', appendLF=FALSE)
  load(OCC_FILE)
}
message('done.')


## for citation (grr, GBIF makes this hard)
# keys <- unlist(map(uniq_species_df$taxon_key, 'key'))
# dwld <- occ_download(
#   pred_in("taxonKey", keys),
#   format = "SIMPLE_CSV",
#   user="vsbuffalo", pwd="hidden", email="vsbuffalo@gmail.com")


# combine in alphas
DEFAULT_ALPHA <- c(8, 70)
RANGE_DIR <- 'total_range_plots'
RAW_RANGES <- 'gbif_inferred_ranges_raw.Rdata'
if (!file.exists(RAW_RANGES)) {
  ### Infer terrestriality
  # filter entries without lat/long and infer terrestriality 
  ranges <- uniq_species_df %>% 
              #filter(species == 'Anoplopoma fimbria') %>%
              #filter(species == 'Aptenodytes patagonicus') %>%
              mutate(valid_geo = map_lgl(occs, valid_occurs)) %>%
              filter(valid_geo) %>%
              distinct(species, .keep_all=TRUE) %>%
              mutate(is_terrestrial = map_lgl(occs, get_terrestriality))

  ## Manually fix things
  # merge in the manually fixed data and fix the terrestriality flags
  source('./species_range_fixes.r')
  ranges <- ranges %>% 
    mutate(genus = extract_genus(species)) %>%
    mutate(is_terrestrial_fixed =
           pmap_lgl(list(genus, species, is_terrestrial),
                    function(g, s, it) {
                     case_when((g %in% marine_gnr) || 
                               (s %in% marine_sps) ~ FALSE,
                               (g %in% terrestrial_gnr) ||
                               (s %in% terrestrial_sps) ~ TRUE,
                               TRUE ~ it)
                    })) %>%
    mutate(is_terrestrial_orig = is_terrestrial, 
           is_terrestrial = is_terrestrial_fixed) %>%
    select(-is_terrestrial_fixed)

  # how many terrestriality fixes?
  table(ranges$is_terrestrial != ranges$is_terrestrial_orig)

  ### Infer Range
  ranges <- ranges %>%  
   mutate(range_area=pmap(list(occs, species, genus, is_terrestrial),
                    function(x, sp, gn, t) {
                      clean_name <- gsub(' +', '_', gsub('/', '_', sp))
                      pdf(sprintf('%s/%s.pdf', RANGE_DIR, clean_name))
                      if (length(x) == 0) {
                        # another degenerate case: just plot an empty plot 
                        maps::map('world', lwd=0.3, col='gray10')
                        dev.off()
                        return(NA)
                      }

                      alpha <- DEFAULT_ALPHA
                      if (sp %in% names(custom_alphas_sps))
                        alpha <- custom_alphas_sps[[sp]]
                      if (gn %in% names(custom_alphas_gnr))
                        alpha <- custom_alphas_gnr[[gn]]
                      # if any alphas are NA, set to default
                      alpha[is.na(alpha)] <- DEFAULT_ALPHA[is.na(alpha)]
                      if (gn %in% no_contrain_gnr || sp %in% no_contrain_sps)
                        constrain <- FALSE
                      else
                        constrain <- TRUE
                      areas <- map_dbl(x, infer_range, alpha=alpha, 
                                      is_terrestrial=t, species=sp, 
                                      constrain=constrain)
                      dev.off()
                      return(areas)
                    }
                    )) %>%
  mutate(n_occ = map_int(occs, count_occs))

  save(ranges, file=RAW_RANGES)
} else {
  message('loading raw range data...   ', appendLF=FALSE)
  load(RAW_RANGES)
}
message('done.')

ranges_df <- ranges %>% select(-occs, -taxon_key, -valid_geo) %>%
               mutate(range = map_dbl(range_area, max, na.rm=TRUE)) %>%
               select(-range_area)

## Manual fixes
# to ensure we are assining *over*, or correcting a value, 
# not creating a new entry, we use this
is_entry <- function(species) stopifnot(any(ranges_df$species == species))

# Human range size is just wrong. Occurrences weren't logged everywhere.
earth_land_area <- 510.1e6 # km^2
antartica_land_area <- 13.21e6 # km^2
human_range <- earth_land_area - antartica_land_area
is_entry('Homo sapiens')
ranges_df[ranges_df$species == 'Homo sapiens', ]$range <- human_range

# Nasonia vitripennis 
# see: https://www.nature.com/articles/hdy2009160/figures/1
us_land_area <- 9.834e6 # km^2
europe_land_area <- 10.18e6  # km^2
nv_range <- us_land_area + europe_land_area
is_entry('Nasonia vitripennis')
ranges_df[ranges_df$species == 'Nasonia vitripennis', ]$range <- nv_range

# Caenorhabditis
# https://www.wolframalpha.com/input/?i=area+of+california+%2B+area+of+washington+state+%2B+area+of+oregon+
# based on Frezel and Felix (2015)
western_us <- 863400 # km^2 
celegans_range <- western_us + europe_land_area
is_entry('Caenorhabditis elegans')
ranges_df[ranges_df$species == 'Caenorhabditis elegans', ]$range <- celegans_range

# Pinus balfouriana
# this was done manually in photoshop from the image 
# https://en.wikipedia.org/wiki/Pinus_balfouriana#/media/File:Pinus_balfouriana_range_map_1.png
ca_area <- 423971 # km^2
pinus_bal_sel <- (174016 - 172459) / 174016 
is_entry('Pinus balfouriana')
ranges_df[ranges_df$species == 'Pinus balfouriana', ]$range <- ca_area* pinus_bal_sel

# Pinus_massoniana
# This is an old map -- I base the area off korea
korea_rel_area <- 9758
pinus_mas_rel_area <- 62666
korea_area <- 100210 + 120540
pinus_mas_area <- korea_area * pinus_mas_rel_area / korea_rel_area
is_entry('Pinus massoniana')
ranges_df[ranges_df$species == 'Pinus massoniana', ]$range <- pinus_mas_area

# Pinus pinaster
pinus_pins_rel_area <- 20880
germany_rel_area <- 101984
germany_area <- 357000 # km^2
pinus_pins_area <- pinus_pins_rel_area / germany_rel_area * germany_area
is_entry('Pinus pinaster')
ranges_df[ranges_df$species == 'Pinus pinaster', ]$range <- pinus_pins_area

# Crassostrea gigas's distribution is coastal, which alpha hull 
# measures have trouble with. 
# 
crass_rel_area <- 15909
australia_rel_area <- 3429
australia_area <- 7.741e6
crass_range <- crass_rel_area / australia_rel_area * australia_area
is_entry('Crassostrea gigas')
ranges_df[ranges_df$species == 'Crassostrea gigas', ]$range <- crass_range

# Cystodytes dellechiajei's automated range is way off. 
# This is because it is a coastal species and the occurence data
# is being overfit. 
# I have used the map here: https://www.sealifebase.ca/summary/Cystodytes-dellechiajei.html#
# and used photoshop to estimate the red areas
cyst_rel_area <- 9538 
australia_rel_area <- 4094
cyst_range <- cyst_rel_area / australia_rel_area * australia_area
is_entry('Cystodytes dellechiajei')
ranges_df[ranges_df$species == 'Cystodytes dellechiajei', ]$range <- cyst_range


# Ciona intestinalis's automated range is way off. 
# https://www.aquamaps.org/receive.php?type_of_map=regular
# note: this isn't that far off from the inferred range from occurrence data, which is reassuring
ciona_i_rel_area <- 2442
australia_rel_area <- 4139
ciona_i_range <- ciona_i_rel_area / australia_rel_area * australia_area
is_entry('Ciona intestinalis')
ranges_df[ranges_df$species == 'Ciona intestinalis', ]$range <- ciona_i_range

# Hippocampus Kuda is off beacuse it's coastal
# https://www.aquamaps.org/receive.php?type_of_map=regular
# note: this isn't that far off from the inferred range from occurrence data, which is reassuring
# note, this leads to N ~ 10^9
# this seems high but searching around, I've found
# reports of 1e7 dead seahorses being poached
# so it's not unreasonable.
# https://abcnews.go.com/Technology/wireStory/peruvian-authorities-123-million-dried-seahorses-seized-66066372
hk_rel_area <- 3644
australia_rel_area <- 4147  # roughly same map source as above
hk_range <- hk_rel_area / australia_rel_area * australia_area
is_entry('Hippocampus kuda')
ranges_df[ranges_df$species == 'Hippocampus kuda', ]$range <- hk_range


# from https://www.iucnredlist.org/species/2815/123789863
is_entry('Bison bison')
ranges_df[ranges_df$species == 'Bison bison', ]$range <- 143253

# Some ranges were very poorly inferred, due mostly do low 
# occurrence data
ranges_df[ranges_df$genus == "Chlamydomonas", ]$range <- NA
ranges_df[ranges_df$genus == "Plasmodium", ]$range <- NA
ranges_df[ranges_df$genus == "Neurospora", ]$range <- NA
ranges_df[ranges_df$genus == "Phytophthora", ]$range <- NA
ranges_df[ranges_df$genus == "Paracoccidioides", ]$range <- NA
ranges_df[ranges_df$genus == "Histoplasma", ]$range <- NA
ranges_df[ranges_df$genus == "Saccharomyces", ]$range <- NA

# Let's talk about Chlorocebus aethiops.
# Looking at the inferred range, it's clear to me that 
# these are generally Chlorocebus sps. labeled as Chlorocebus.
# see: https://www.researchgate.net/figure/Distribution-of-African-green-monkeys-Chlorocebus-species-in-Africa-Chlorocebus_fig5_304030985
# and https://en.wikipedia.org/wiki/Grivet#/media/File:Grivet_area.png
# this is manually fixed using the Wikipedia range.
grivet_rel_area <- 1787
saudia_arabia_rel_area <- 3136
saudia_arabia_area <- 1.961e6 # km^2 -- wolfram alpha
grivet_range <- grivet_rel_area / saudia_arabia_rel_area * saudia_arabia_area
is_entry('Chlorocebus aethiops')
ranges_df[ranges_df$species == 'Chlorocebus aethiops', ]$range <- grivet_range 


# Like Chlorocebus, Callithrix jacchus occurrences seem to 
# generally include all Callithrix species.
# Compare: https://commons.wikimedia.org/wiki/File:Distribution_Callithrix.jpg
# with https://commons.wikimedia.org/wiki/File:Callithrix_jacchus_distribution.svg
common_marmoset_rel_area <- 1630
brazil_rel_area <- 12730 + common_marmoset_rel_area
brazil_area <- 8.515e6 # km^2 wolfram alpha
common_marmoset_area <- common_marmoset_rel_area / brazil_rel_area * brazil_area 
is_entry('Callithrix jacchus')
ranges_df[ranges_df$species == 'Callithrix jacchus', ]$range <- common_marmoset_area


# Pongo species -- Pongo abelii and P. pygmaeus
# (in the Leffler data they are combined) 
# Pongo species are way off -- here alpha hull approach is too liberal for
# these endagered species. Their ranges are (sadly) much more fragmented. IUCN
# redlist data does a better job:
# https://www.iucnredlist.org/species/121097935/123797627#population
# 
# I use the data on that page, from Wich et al. (2016)
is_entry('Pongo abelii')
ranges_df[ranges_df$species == 'Pongo abelii', ]$range <- 16775 # km^2
# https://www.iucnredlist.org/species/17975/123809220#text-fields
is_entry('Pongo pygmaeus')
ranges_df[ranges_df$species == 'Pongo pygmaeus', ]$range <- 97716 #km^2


## Termites
# the occurrence data is bad for these two species. Works well for flavipes
#https://www.researchgate.net/figure/Geographic-distribution-of-Reticulitermes-species-Isoptera-Rhinotermitidae-in-western_fig1_226108638
is_entry('Reticulitermes grassei')
is_entry('Reticulitermes lucifugus')
ret_grassei_rel <- 24005
ret_lucif_rel <- 10064
austria_area <- 32386
austria_rel_area <- 4385
ranges_df[ranges_df$species == 'Reticulitermes grassei', ]$range <- ret_grassei_rel / austria_rel_area * austria_area
ranges_df[ranges_df$species == 'Reticulitermes lucifugus', ]$range <- ret_lucif_rel / austria_rel_area * austria_area


## Daphnia
# Daphnia is hard -- clearly, though, terrestrial inference is way off since
# these live in freshwater. As a *very* rough approximation we say that 4% of
# land is lakes (this is true globally:
# https://blog.nationalgeographic.org/2014/09/15/117-million-lakes-found-in-latest-world-count/) 
PERCENT_LAKE <- 0.037 
daphnia_ranges <- ranges_df[ranges_df$genus == 'Daphnia', ]$range 
ranges_df[ranges_df$genus == 'Daphnia', ]$range  <- PERCENT_LAKE * daphnia_ranges


## Shrimp
# There were poorly inferred from GBIF data for a known and expected 
# reason: they're primarily coastal
# The manual fixes are based off a map (downloaded
# from http://en.aquaculture.ifremer.fr/Sectors/Crustacean-sector/Discoveries/Shrimps
# and in range_maps)
# 
# first image in the range map pair
australia_rel_area <- 5307
mars_rel_area <- 3102
mars_area <- mars_rel_area  / australia_rel_area * australia_area
is_entry('Marsupenaeus japonicus')
ranges_df[ranges_df$species == 'Marsupenaeus japonicus', ]$range <- mars_area

# second image in the range map pair
lito_rel_area <- 6638
australia_rel_area2 <- 5203 # 2, for the second map in this range map pair
lito_area <- lito_rel_area / australia_rel_area2 * australia_area
is_entry('Litopenaeus vannamei')
ranges_df[ranges_df$species == 'Litopenaeus vannamei', ]$range <- lito_area

panaeus_rel_area <- 7398
australia_rel_area2 <- 5203 # 2, for the second map in this range map pair
panaeus_area <- panaeus_rel_area / australia_rel_area2 * australia_area
is_entry('Penaeus monodon')
ranges_df[ranges_df$species == 'Penaeus monodon', ]$range <- panaeus_area

## Canis lupus
# the occurrence dataset has a bunch of odd entries. I instead 
# use the IUCN redlist range maps
australia_rel_area <- 1632
canisl_rel_area <- 31281
canisl_area <- canisl_rel_area  / australia_rel_area * australia_area
is_entry('Canis lupus')
ranges_df[ranges_df$species == 'Canis lupus', ]$range <- canisl_area


## Argopecten irradians -- bay scallop
# like shrimp, these are coastal and their range is poorly inferred by 
# GBIF occurrence data using an alpha hull
# Despite a few occurrences across the globe, this is an Atlantic coast species
# as far as I can tell.
# I use the range map from here: 
# https://www.sciencedirect.com/topics/agricultural-and-biological-sciences/argopecten-irradians
# The range map was too roughly colored, so I've very rougly traced it
# in bright green (again, on log scale this noise is fine)
argo_rel_area <- 6864
bahamas_rel_area <- 369
bahamas_area <- 13880 # km2, wolfram alpha
argo_area <- argo_rel_area / bahamas_rel_area * bahamas_area
is_entry('Argopecten irradians')
ranges_df[ranges_df$species == 'Argopecten irradians', ]$range <- argo_area

## bird range fixes
# Aquila clanga/pomarina -- from http://datazone.birdlife.org/species/factsheet/lesser-spotted-eagle-clanga-pomarina -- geometric mean of both species, since their diversity data is averaged
is_entry('Aquila clanga')
ranges_df[ranges_df$species == 'Aquila clanga', ]$range <- 5340000
is_entry('Aquila pomarina')
ranges_df[ranges_df$species == 'Aquila pomarina', ]$range <- 18100000

# not a huge difference, but I'll use the birdlife value
is_entry('Ficedula albicollis')
ranges_df[ranges_df$species == 'Ficedula albicollis', ]$range <- 4430000
# this range is muuuch smaller according to birdlife.org
is_entry('Taeniopygia guttata')
ranges_df[ranges_df$species == 'Taeniopygia guttata', ]$range <- 308000
# again, not a huge difference but I'll use the birdlife.org values
is_entry('Zonotrichia albicollis')
ranges_df[ranges_df$species == 'Zonotrichia albicollis', ]$range <- 8040000


## Equus ferus przewalskii
# from: https://ielc.libguides.com/sdzg/factsheets/przewalskishorse/behavior
# Hustai National Park: 120-2,400 ha (King and Gurnell 2005)
# Great Gobi B Strictly Protected Area: 150-825 km² (Kaczensky et al. 2008).
# note that the Corbett-Detig dataset lists this range as 10^7.13, which is around
# what regular horsies have... and this just seems wrong.
# 
# since this is a correction to Corbett Detig's data, which hasn't been 
# merged in yet, I add this row first.
eqfp <- tibble(key = "Equus ferus przewalskii", 
               species = "Equus ferus przewalskii",
               dataset = "Corbett-Detig", is_terrestrial = TRUE, 
               genus="Equus", 
               is_terrestrial_orig = TRUE, n_occ = NA, range = NA)
ranges_df <- rbind(ranges_df, eqfp)
is_entry('Equus ferus przewalskii')
ranges_df[ranges_df$species == 'Equus ferus przewalskii', ]$range <- geomean(0.01 * c(120, 2400)) + geomean(c(150, 825))


## Metriaclima zebra, an endemic cichlid
# Issue here is it's endemc in Lake Malawi, with few occurrences logged.
# I just put the whole Lake Malawi area
lake_malawi_area <- 29600 # km2
is_entry('Metriaclima zebra')
ranges_df[ranges_df$species == 'Metriaclima zebra', ]$range <- lake_malawi_area

## Fenneropenaeus chinensis, Chinese white shrimp
# this is an obvoius outlier in the data, and has too few occurrences to
# accurately infer the range of
is_entry('Fenneropenaeus chinensis')
ranges_df[ranges_df$species == 'Fenneropenaeus chinensis', ]$range <- NA

# Let's fix the Drosophilids. 
# The data is from 
# https://web.archive.org/web/20170710053101/http://www.drosophila-speciation-patterns.com/
# where I've found a spreadsheet of the absolute range sizes
drosoph <- read_xlsx('./RAW_dataset_Jan_06_2012.xlsx')[, c(1, 27)] %>%
             rename(species=`Sp 1`, 
                    range=`absolute range (sp.1) * in square km`) %>%
             mutate(species = paste('Drosophila', species)) %>%
             mutate(range = as.numeric(range)) %>%
             filter(!is.na(range)) %>% distinct()


ranges_df <- ranges_df %>% left_join(drosoph, by='species') %>%
                 mutate(range = pmax(range.x, range.y, na.rm=TRUE))  %>%
                 select(-range.y, -range.x)

write_tsv(ranges_df, "gbif_ranges.tsv")

## some diagnostic images

# leffd %>%
#   group_by(species) %>% summarize_at(c('size', 'range', 'diversity'), mean, na.rm=TRUE) %>%
#   ggplot(aes(size, range, color=genus)) + scale_x_log10() +  geom_point() +
#     scale_y_log10() + #geom_smooth(method='lm') + 
#     geom_text_repel(aes(size, range, label=species), size=1.8) + 
#     theme(legend.position="none")

# leffd %>%
#   ggplot(aes(range / size, diversity)) + scale_x_log10() + 
#     scale_y_log10() + geom_smooth(method='lm') + 
#     geom_text_repel(aes(range/size, diversity, label=species), size=1.8) 


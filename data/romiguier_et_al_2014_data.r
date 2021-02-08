library(tidyverse)
library(readxl)
library(ritis)

dr_raw <- read_xls('romiguier_et_al_2014_supptable_2.xls')

colnames(dr_raw) <- c("species", "common_name", "family", "order",
                      "class", "phylum", "number_of_inds", "contigs_number",
                      "average_length_of_contigs", "SNP_number", "piS",
                      "piN", "piNpiS", "fit", "invasive_status", 
                      "maximal_geoographic_distance",
                      "mean_geographic_distance", "mean_distance_from_equator")

# remove quotes around names
cleanup_text <- function(x) gsub('"([^"]+)"', '\\1', x)

clean_species <- function(x) {
  sub('([^_]+)_([^_]+).*', '\\1 \\2', x)
}

dr <- dr_raw %>% mutate_at(vars(species:phylum), cleanup_text) %>%
         mutate(species=map_chr(species, clean_species)) %>%
         mutate(piS=as.numeric(piS), 
                piN=as.numeric(piN))


# Supplementary table 3 data, extracted from PDF using Tabula 
# and then cleaned up by hand
s3_files <- list.files('.', pattern='romiguier_et_al_supptable_3_.*\\.csv', 
              full.names=TRUE)

s3_type <- sub('.*_([^\\.]+)\\.csv', '\\1', s3_files)
names(s3_files)  <- s3_type

s3d_raw <- Map(function(x, n) {
                 d <- read_csv(x)
                 colnames(d)[3:4] <- paste(n, colnames(d)[3:4], sep='_')
                 d
              }, s3_files, s3_type)

# join everything
drs <- reduce(s3d_raw, inner_join, by='species') %>% 
           # convert to kg
           mutate(mass = bodymass_g / 1000)  %>%
           # convert to m
           mutate(size = size_cm / 100) %>%
           # convert to m
           mutate(propagule = propagule_size_cm / 100) %>%
           mutate(log10_size=log10(size),
                  log10_mass=log10(mass),
                  log10_propagule=log10(propagule)) %>%
           mutate(species=map_chr(species, clean_species)) %>%
           left_join(dr) %>% 
        mutate(genus = map_chr(species, ~ strsplit(., ' ')[[1]][1]))

# There is a very odd size entry for Lineus longissimus -- 10m.
# They can get this long, but Wikipedia says they are usually "5 to 10mm"
# (https://en.wikipedia.org/wiki/Lineus_longissimus)
# Given they were an outlier on range / size plots, I've made this change
drs$size[drs$species == 'Lineus longissimus'] <- 0.0075


# there is another error in the Romiguier et al size data, as 
# far as I can tell. This page says Reticulitermes flavipes

# gets to be 4-5mm, not cm
# https://animaldiversity.org/accounts/Reticulitermes_flavipes/#:~:text=flavipes%20soldiers%20are%20slightly%20bigger,of%20up%20to%209cm%20long.
drs$size[drs$species == 'Reticulitermes flavipes'] <- 0.0075
# same with R. lucifigus: https://mapadetermitas.org/uploads/library/LIB_80D061-8E6B05-CC523F-E0DFFD-97A0E5-AF2292.pdf
drs$size[drs$species == 'Reticulitermes lucifigus'] <- 0.012
# and R. grassei too.
# source: https://www.researchgate.net/figure/The-species-Reticulitermes-grassei-Clement-A-Alate-caste-characteristically-dark_fig5_237044135#:~:text=The%20alates%20are%20dark%20in,5).
drs$size[drs$species == 'Reticulitermes grassei'] <- 0.005


write_tsv(drs, 'romiguier_et_al_2014_updated.tsv')

# # there's a tight relation between mass and size
# ggplot(drs) + geom_point(aes(log10_mass, log10_size))

# ggplot(drs) + geom_point(aes(log10(size), log10(fecundity), color=phylum))

# # propagule size and fecundity
# ggplot(drs) + geom_point(aes(log10(fecundity), log10_propagule, color=phylum))

# ggplot(drs) + geom_point(aes(log10(fecundity), log10(piS), color=phylum))

# drs %>% group_by(family) %>%
#   summarize(log10_size = mean(log10_size, na.rm=TRUE),
#             log10_fecundity = mean(log10(fecundity), na.rm=TRUE)) %>%
# ggplot() + geom_point(aes(log10_size, log10_fecundity, color=family)) +
#    geom_smooth(aes(log10_size, log10_fecundity), method='lm')  

# ggplot(drs) + geom_point(aes(log10_size, log10(fecundity), color=phylum)) +
#    geom_smooth(aes(log10_size, log10(fecundity)), method='lm')  


# ggplot(drs) + geom_point(aes(log10_size, log10_propagule, color=phylum)) + 
#   geom_text_repel(aes(log10_size, log10_propagule, label=species))


# ggplot(drs) + geom_point(aes(log10_size, log10(piS), color=phylum))


# ggplot(drs) + geom_point(aes(log10(fecundity), log10(piS), color=phylum))
# ggplot(drs) + geom_point(aes(log10(longevity_years), log10(piS), color=phylum))

# ggplot(drs) + geom_point(aes(log10_size, piN/piS, color=phylum))

# ggplot(drs) + geom_point(aes(log10_propagule, log10(piS), color=phylum)) + 
#   geom_text_repel(aes(log10_propagule, log10(piS), label=species))

# ggplot(drs) + 
#   geom_point(aes(x=log10_size, y=piN, color=phylum)) +
#   geom_segment(aes(x=log10_size, y=piN, xend=log10_size, yend=piS, color=phylum), 
#                            arrow=arrow(length = unit(0.009, "npc"))) + scale_y_log10()


# table(dr$species %in% dc$species)


# dc %>% inner_join(drs, by='species') %>% 
#   ggplot() + geom_point(aes(map_length, log10_propagule)) + 
#   geom_text_repel(aes(map_length, log10_propagule, label=species))


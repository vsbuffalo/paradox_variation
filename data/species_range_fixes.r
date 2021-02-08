library(tidyverse)

terrestriality_fixes <- tribble(~ name, ~ is_terrestrial,
        'Anopheles scanloni', TRUE,
        'Arabidopsis', TRUE,
        'Aphanarthrum subglabrum/glabrum', TRUE,
        'Ciona', FALSE, 
        'Caenorhabditis', TRUE,
        'Crassostrea gigas', FALSE,
        'Cryptomeria japonica', TRUE,
        # note: daphnia is freshwater, so should be found on continents
        'Daphnia', TRUE,
        'Drosophila', TRUE,
        # occurrences of gray whales are from coasts, so this is 
        # inferred as a terrestrial, yes, land gray whales... *facepalm* 
        'Eschrichtius robustus', FALSE,
        'Ficedula albicollis', TRUE,
        # Gasterosteus aculeatus / stickelback is classified as terrestrial here --
        # mostly inland (no changes)
        'Hoolock leuconedys', TRUE,
        # chinook is hard-- it's both. We assign as terretrial here
        # as that gives a better estimate
        'Oncorhynchus tshawytscha', TRUE,
        # all urchins are not terrestrial
        'Strongylocentrotu', FALSE,
        'Saccharomyces', TRUE,
        'Bombyx', TRUE,
        'Agaricus', TRUE,
        'Agrostis', TRUE,
        'Agrostis', TRUE,
        'Amaranthus', TRUE,
        'Ananas', TRUE,
        'Bracon', TRUE,
        'Brassica', TRUE,
        'Camelina', TRUE,
        'Capsicum', TRUE,
        'Carica', TRUE,
        'Cocos', TRUE,
        'Coffea', TRUE,
        'Cucurbita', TRUE,
        'Cynodon', TRUE,
        'Daucus', TRUE,
        # penguin
        'Eudyptes', TRUE,
        'Festuca', TRUE,
        'Gossypium', TRUE,
        'Hibiscus', TRUE,
        'Hippocampus', FALSE,
        'Lactuca', TRUE,
        'Lineus', FALSE,
        'Litopenaeus', FALSE,
        'Lolium', TRUE,
        'Lupinus', TRUE,
        'Macropus', TRUE,
        'Malus', TRUE,
        'Mangifera', TRUE,
        'Medicago', TRUE,
        'Mimulus', TRUE,
        'Mytilus', FALSE,
        'Oncorhynchus mykiss', TRUE,
        'Oncorhynchus kisutch', FALSE,
        'Oncorhynchus nerka', FALSE,
        'Oryza', TRUE,
        'Pinus', TRUE,
        'Populus', TRUE,
        'Rubus', TRUE,
        'Salmo salar', TRUE,
        'Sepia officinalis', FALSE,
        'Seriola quinqueradiata', FALSE,
        'Sorghum', TRUE,
        'Strongylocentrotus', FALSE,
        'Vaccinium', TRUE,
        'Vicia', TRUE,
        'Vitis', TRUE, 
        'Quercus', TRUE,
        'Crassostrea', FALSE
)

alpha_fixes <- tribble(~name, ~alpha,
       'Crassostrea', c(NA, 10),
       'Eunicella', c(NA, 10),
       'Haliotis', c(NA, 10),
       # seahorse -- coastal
       'Hippocampus', c(NA, 10),
       'Hippoglossus', c(NA, 40),
       'Leptogorgia', c(NA, 40),
       'Mytilus', c(NA, 10),
       'Drosophila', c(80, NA),
       'Oncorhynchus', c(10, 40),
       'Daphnia magna', c(10, NA),
       'Ciona', c(NA, 10),
       'Penaeus monodon', c(NA, 40),
       'Strongylocentrotus droebachiensis', c(NA, 40),
       # this was an outlier in the range/size diversity plot
       # and I compared to 
       'Populus tremula', c(30, NA),
       'Nasonia', c(40, NA),
       # cuttlefish -- coastal
       'Sepia officinalis', c(NA, 10),
       # domesticated maize, cattle, and cucumber need a smaller alpha
       'Zea mays', c(4, NA),
       'Bos taurus', c(4, NA),
       'Cucumis sativus', c(4, NA)
)


extract_genus <- function(x) {
  map_chr(strsplit(x, ' +'), 1)
}

whole_genus <- function(x) {
  map_int(strsplit(x, ' +'), length) == 1
}

terrestriality_fixes <- terrestriality_fixes %>% 
                          mutate(genus = extract_genus(name),
                                 whole_genus = whole_genus(name),
                                 species = ifelse(whole_genus, NA, name))

# create sets
terrestrial_sps <- terrestriality_fixes %>% filter(!whole_genus, is_terrestrial) %>% 
                       filter(!is.na(species)) %>% pull(species)
terrestrial_gnr <- terrestriality_fixes %>% filter(whole_genus, is_terrestrial) %>% 
                       pull(genus)
marine_sps <- terrestriality_fixes %>% filter(!whole_genus, !is_terrestrial) %>% 
                       filter(!is.na(species)) %>% pull(species)
marine_gnr <- terrestriality_fixes %>% filter(whole_genus, !is_terrestrial) %>% 
                       pull(genus)

# create alpha sets
alpha_fixes <- alpha_fixes %>% mutate(genus = extract_genus(name),
                                 whole_genus = whole_genus(name),
                                 species = ifelse(whole_genus, NA, name))

custom_alphas_gnr <- 
  structure(
  alpha_fixes %>% filter(whole_genus) %>% pull(alpha),
  .Names=alpha_fixes %>% filter(whole_genus) %>% pull(genus))

custom_alphas_sps <- structure(
  alpha_fixes %>% filter(!whole_genus) %>% pull(alpha),
  .Names=alpha_fixes %>% filter(!whole_genus) %>% pull(species))

# a few species to not constrain -- anadramous fishes mostly
no_contrain_gnr <- c('Oncorhynchus')
no_contrain_sps <- c()


